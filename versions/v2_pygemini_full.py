#!/usr/bin/env python3
"""PY-GEMINI - core-genome and recombination-aware phylogeny from protein FASTA files.

Command-line rewrite of the original PY-GEMINI.py. The scientific method, the analysis
parameters and the output layout are unchanged; see AI-ASSISTANCE-RECORD.md for a full
account of what differs from the original script and why.

Pipeline
    A  prepare a working copy of the input          prepare_working_copy
    B  prepare the output directory                 prepare_output_directory
       normalise the input file names               normalise_input_filenames
    C  all-vs-all BLASTP                            run_all_vs_all_blast
    D  homolog extraction                           extract_homologs
    E  choose the reciprocal reference genome       select_reference_genome
    F  write homolog sequences                      write_homolog_sequences
    G  reciprocal hit pairs                         find_reciprocal_pairs
       core gene selection                          select_core_genes
       core gene tables and sequences               write_core_gene_tables
    H  concatenated supermatrix                     build_supermatrix
    I  multiple sequence alignment                  run_alignment
    J  Newick tree                                  build_tree
    K  one-line FASTA / ortholog group transpose    write_single_line_fasta
                                                    transpose_gene_sets
    L  align each ortholog group                    align_gene_sets
    M  PhiPack recombination test                   run_recombination_tests
    N  parse PHI values / split recombinant genes   parse_phi_values
                                                    split_recombinant_genes
                                                    write_nonrecombinant_sequences
    O  alignment of the non-recombinant supermatrix run_alignment
    P  non-recombinant Newick tree                  build_tree

Every stage is a single function so that each can be checked, or reused, on its own.
"""

import argparse
import logging
import os
import shutil
import subprocess
import sys
from typing import NamedTuple

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

__version__ = '2.0.0'

LOG = logging.getLogger('pygemini')

# The original spelled this folder 'homologousPID=>85%Genome'; '>' is not a legal character in a
# Windows path, so the name is sanitised. Contents are unchanged.
HOMO_PID_DIR = 'homologousPID_ge_85pc_Genome'


# --------------------------------------------------------------------------------------
# external programs
# --------------------------------------------------------------------------------------

def find_executable(*names):
    """Return the first of `names` found on PATH, or None."""
    for name in names:
        path = shutil.which(name)
        if path:
            return path
    return None


def resolve_tools():
    """Locate every external program the pipeline needs.

    Returns a dict of tool -> absolute path. Raises SystemExit naming everything missing,
    rather than failing part-way through a long run.
    """
    tools = {
        'makeblastdb': find_executable('makeblastdb'),
        'blastp': find_executable('blastp'),
        'clustalw': find_executable('clustalw', 'clustalw2'),
        'fasttree': find_executable('fasttree', 'FastTree', 'fasttreeMP', 'FastTreeMP'),
        'phi': find_executable('Phi', 'phipack-phi', 'phi', 'PhiPack'),
    }
    missing = sorted(name for name, path in tools.items() if path is None)
    if missing:
        raise SystemExit(
            'Required programs not found on PATH: ' + ', '.join(missing) + '\n'
            'Install them into the active environment, for example:\n'
            '    micromamba install -c conda-forge -c bioconda blast clustalw fasttree phipack\n'
            'or create the bundled environment:\n'
            '    micromamba env create -f environment.yml')
    for name, path in sorted(tools.items()):
        LOG.debug('using %s: %s', name, path)
    return tools


def run_tool(cmd, stdout_path=None, stdin_path=None, mode='w', check=True):
    """Run an external program.

    Arguments are passed as a list so paths containing spaces need no shell quoting.
    A non-zero exit aborts the run unless check=False -- carrying on past a failed tool is
    how an empty alignment ends up presented as a finished tree.
    """
    LOG.debug('run: %s', ' '.join(cmd))
    out = inp = None
    if stdout_path:
        # create the destination folder here so no caller has to remember to
        os.makedirs(os.path.dirname(os.path.abspath(stdout_path)), exist_ok=True)
    try:
        out = open(stdout_path, mode) if stdout_path else None
        inp = open(stdin_path, 'r') if stdin_path else None
        code = subprocess.call(cmd, stdout=out, stdin=inp)
    finally:
        if out:
            out.close()
        if inp:
            inp.close()
    if check and code != 0:
        raise SystemExit(
            '%s failed (exit status %d):\n    %s\n'
            'PY-GEMINI stopped rather than continue with incomplete data.'
            % (os.path.basename(cmd[0]), code, ' '.join(cmd)))
    return code


def require_nonempty(path, description):
    """Abort if a tool exited 0 but produced nothing."""
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        raise SystemExit(
            '%s is empty or was never written: %s\n'
            'Check the matching log file before trusting any earlier output.' % (description, path))


# --------------------------------------------------------------------------------------
# small shared helpers
# --------------------------------------------------------------------------------------

def stem_of(name):
    """Filename without its extension. Never use str.rstrip for this."""
    return os.path.splitext(name)[0]


def path_contains(parent, child):
    """True if `child` is `parent` itself or lies anywhere inside it.

    Comparison is case-insensitive where the platform is, and `rstrip` on the separator keeps
    a drive root ('C:\\') or the POSIX root ('/') working as a parent.
    """
    parent = os.path.normcase(os.path.normpath(parent))
    child = os.path.normcase(os.path.normpath(child))
    return child == parent or child.startswith(parent.rstrip(os.sep) + os.sep)


def blast_row_passes(fields, pident_min, gapopen_max, evalue_max):
    """Apply the identity / gap-open / e-value filter to one BLAST -outfmt 6 row."""
    try:
        return (float(fields[2]) >= pident_min
                and int(fields[5]) <= gapopen_max
                and float(fields[-2]) < evalue_max)
    except (IndexError, ValueError):
        return False


def read_id_list(path):
    """Read a one-identifier-per-line file, ignoring blank lines."""
    with open(path, 'r') as handle:
        return [line.strip() for line in handle if line.strip()]


def write_id_list(path, identifiers):
    """Write one identifier per line."""
    with open(path, 'w') as handle:
        for identifier in identifiers:
            handle.write(identifier + '\n')


def index_fasta(path):
    """Map sequence id -> SeqRecord for one FASTA file.

    Indexing once avoids the original's re-parse of the whole proteome per header, which
    made the sequence-collection stages quadratic.
    """
    return {record.id: record for record in SeqIO.parse(path, 'fasta')}


def copy_records(index, identifiers, source):
    """Fetch `identifiers` from a FASTA index, preserving order. Aborts on a missing id."""
    missing = [i for i in identifiers if i not in index]
    if missing:
        raise SystemExit('%s: %d sequence(s) named in the tables are absent from the FASTA, '
                         'first missing: %s' % (source, len(missing), missing[0]))
    return [SeqRecord(Seq(str(index[i].seq)), id=index[i].id, description=index[i].description)
            for i in identifiers]


# --------------------------------------------------------------------------------------
# stage A/B - directories
# --------------------------------------------------------------------------------------

def prepare_working_copy(input_dir, work_dir):
    """Copy the input proteomes aside.

    BLAST writes database files next to the sequences and the files are renamed to .fasta,
    so every mutation happens here and the user's input directory is never touched.
    """
    if os.path.exists(work_dir):
        shutil.rmtree(work_dir)
    shutil.copytree(input_dir, work_dir)
    LOG.info('working copy: %s', work_dir)
    return work_dir


def prepare_output_directory(output_dir, force):
    """Create the results directory, refusing to clobber existing data unless forced."""
    if os.path.exists(output_dir):
        if not os.path.isdir(output_dir):
            raise SystemExit('Output path exists and is not a directory: ' + output_dir)
        if os.listdir(output_dir):
            if not force:
                raise SystemExit(
                    'Output directory is not empty: %s\n'
                    'Use --force to delete and recreate it, or choose another --output.'
                    % output_dir)
            LOG.warning('--force given, deleting existing %s', output_dir)
            shutil.rmtree(output_dir)
    os.makedirs(output_dir, exist_ok=True)
    LOG.info('results: %s', output_dir)
    return output_dir


def normalise_input_filenames(work_dir):
    """Rename every proteome to <stem>.fasta and return the sorted file list.

    The duplicate-stem check runs BEFORE any rename: 'x.faa' and 'x.fa' both become
    'x.fasta', which would silently destroy one of them.
    """
    entries = sorted(os.listdir(work_dir))
    if not entries:
        raise SystemExit('No files found in the input directory.')
    stems = [stem_of(name) for name in entries]
    duplicates = sorted({s for s in stems if stems.count(s) > 1})
    if duplicates:
        raise SystemExit('These input names collide once the extension is removed: '
                         + ', '.join(duplicates))
    for name in entries:
        target = stem_of(name) + '.fasta'
        if target != name:
            os.rename(os.path.join(work_dir, name), os.path.join(work_dir, target))
    files = sorted(os.listdir(work_dir))
    if len(files) < 2:
        raise SystemExit('At least two protein FASTA files are needed, found %d.' % len(files))
    return files


# --------------------------------------------------------------------------------------
# stage C - BLAST
# --------------------------------------------------------------------------------------

def run_all_vs_all_blast(tools, work_dir, files, blast_dir, log_path, threads):
    """Build a BLAST database per genome and query every other genome against it.

    Produces <blast_dir>/<db stem>/<db stem>_<query stem>.tsv, so column 0 of a row is the
    QUERY genome's gene and column 1 is the DATABASE genome's gene.
    """
    stems = [stem_of(f) for f in files]
    for i, db_file in enumerate(files):
        os.makedirs(os.path.join(blast_dir, stems[i]), exist_ok=True)
        db_path = os.path.join(work_dir, db_file)
        run_tool([tools['makeblastdb'], '-in', db_path, '-dbtype', 'prot',
                  '-title', stems[i], '-out', db_path + 'db'],
                 stdout_path=log_path, mode='a')
        for j, query_file in enumerate(files):
            if i == j:
                continue
            LOG.info('blastp: %s vs %s', stems[j], stems[i])
            run_tool([tools['blastp'],
                      '-query', os.path.join(work_dir, query_file),
                      '-db', db_path + 'db',
                      '-out', os.path.join(blast_dir, stems[i],
                                           stems[i] + '_' + stems[j] + '.tsv'),
                      '-num_threads', str(threads), '-outfmt', '6'])
    return stems


# --------------------------------------------------------------------------------------
# stage D - homologs
# --------------------------------------------------------------------------------------

def extract_homologs(blast_dir, homolog_dir, stems, pident_min, gapopen_max, evalue_max):
    """Keep the genes of each genome that are hit from EVERY other genome.

    The per-file hit sets are intersected. The original summed occurrences across files,
    so a gene hit twice by one genome and never by another still reached the threshold.
    """
    os.makedirs(homolog_dir, exist_ok=True)
    counts = {}
    for stem in stems:
        genome_dir = os.path.join(blast_dir, stem)
        per_file = []
        for name in sorted(os.listdir(genome_dir)):
            hits = set()
            with open(os.path.join(genome_dir, name), 'r') as handle:
                for line in handle:
                    fields = line.rstrip('\n').split('\t')
                    if blast_row_passes(fields, pident_min, gapopen_max, evalue_max):
                        hits.add(fields[1])
            per_file.append(hits)
        if len(per_file) == len(stems) - 1 and per_file:
            common = sorted(set.intersection(*per_file))
        else:
            common = []
        write_id_list(os.path.join(homolog_dir, stem + '.txt'), common)
        counts[stem] = len(common)
        LOG.info('%s: %d homologs', stem, len(common))
    return counts


# --------------------------------------------------------------------------------------
# stage E - reference genome
# --------------------------------------------------------------------------------------

def select_reference_genome(homolog_counts):
    """The genome with the fewest homologs becomes the reciprocal reference."""
    if not homolog_counts:
        raise SystemExit('No homolog sets were produced.')
    if min(homolog_counts.values()) == 0:
        empty = sorted(s for s, n in homolog_counts.items() if n == 0)
        raise SystemExit(
            'These genomes share no homolog with all the others, so there is no core genome: '
            + ', '.join(empty) + '\nRelax --pident / --evalue, or check the inputs.')
    reference = min(sorted(homolog_counts), key=lambda s: homolog_counts[s])
    LOG.info('reference genome: %s (%d homologs)', reference, homolog_counts[reference])
    return reference


# --------------------------------------------------------------------------------------
# stage F - homolog sequences
# --------------------------------------------------------------------------------------

def write_homolog_sequences(work_dir, homolog_dir, homolog_seq_dir, stems):
    """Write the FASTA of each genome's homologs."""
    os.makedirs(homolog_seq_dir, exist_ok=True)
    for stem in stems:
        identifiers = read_id_list(os.path.join(homolog_dir, stem + '.txt'))
        index = index_fasta(os.path.join(work_dir, stem + '.fasta'))
        records = copy_records(index, identifiers, stem)
        SeqIO.write(records, os.path.join(homolog_seq_dir, stem + '.fasta'), 'fasta')
    return homolog_seq_dir


# --------------------------------------------------------------------------------------
# stage G - reciprocal pairs
# --------------------------------------------------------------------------------------

def find_reciprocal_pairs(blast_dir, homolog_dir, reciprocal_dir, stems, reference,
                          pident_min, gapopen_max, evalue_max):
    """Pair the reference genome with each other genome using hits seen in BOTH directions.

    The match is forced ONE-TO-ONE in both directions: a reference gene backs at most one
    partner and a partner backs at most one reference gene. Without the second constraint two
    near-identical paralogs in the reference could both claim the same partner, putting one
    sequence into the supermatrix twice. Rows are consumed in descending bitscore order so the
    best reciprocal hit wins and the outcome does not depend on the order BLAST printed rows.

    Returns {other stem: {reference gene: partner gene}}.
    """
    full_dir = os.path.join(reciprocal_dir, 'reciprocal_blast_fullGenome')
    homo_rows_dir = os.path.join(reciprocal_dir, 'reciprocal_blast_homologousGenome')
    pid_dir = os.path.join(reciprocal_dir, HOMO_PID_DIR)
    renamed_dir = os.path.join(reciprocal_dir, 'reciprocal_renamed')
    genome_dir = os.path.join(reciprocal_dir, 'reciprocal_genome')
    for path in (full_dir, homo_rows_dir, pid_dir, renamed_dir, genome_dir):
        os.makedirs(path, exist_ok=True)

    homologs = {s: set(read_id_list(os.path.join(homolog_dir, s + '.txt'))) for s in stems}

    def homologous_rows(db_stem, query_stem):
        """Copy one BLAST file aside and keep rows whose subject is a homolog of the db genome."""
        name = db_stem + '_' + query_stem + '.tsv'
        source = os.path.join(blast_dir, db_stem, name)
        shutil.copyfile(source, os.path.join(full_dir, name))
        kept = []
        with open(source, 'r') as handle:
            for line in handle:
                fields = line.rstrip('\n').split('\t')
                if len(fields) > 1 and fields[1] in homologs[db_stem]:
                    kept.append(line)
        with open(os.path.join(homo_rows_dir, name), 'w') as handle:
            handle.writelines(kept)
        return kept, name

    def dump_pairs(path, pairs):
        with open(path, 'w') as handle:
            handle.write(''.join(p[0] + ':' + p[1] + '\n' for p in pairs))

    pair_maps = {}
    for other in [s for s in stems if s != reference]:
        # forward: database = reference, query = other -> col 0 is other's gene, col 1 reference's
        forward_rows, forward_name = homologous_rows(reference, other)
        # reverse: database = other, query = reference -> col 0 is reference's gene, col 1 other's
        reverse_rows, reverse_name = homologous_rows(other, reference)

        # pairs are carried as tuples and only joined with ':' when written, so an identifier
        # containing a colon cannot be mis-parsed on the way back in
        forward_pairs = []
        for line in forward_rows:
            fields = line.rstrip('\n').split('\t')
            if blast_row_passes(fields, pident_min, gapopen_max, evalue_max):
                forward_pairs.append((fields[1], fields[0]))
        reverse_pairs = []
        for line in reverse_rows:
            fields = line.rstrip('\n').split('\t')
            if blast_row_passes(fields, pident_min, gapopen_max, evalue_max):
                try:
                    bitscore = float(fields[-1])
                except (IndexError, ValueError):
                    bitscore = 0.0
                reverse_pairs.append((fields[0], fields[1], bitscore))

        dump_pairs(os.path.join(pid_dir, forward_name), forward_pairs)
        dump_pairs(os.path.join(pid_dir, reverse_name), reverse_pairs)
        # inspection copies. NOTE: both are written reference:partner. The original wrote the
        # '1.' file the other way round; the orientation is unified here.
        dump_pairs(os.path.join(renamed_dir, '1.' + reference + '_' + other + '.tsv'),
                   forward_pairs)
        dump_pairs(os.path.join(renamed_dir, '2.' + reference + '_' + other + '.tsv'),
                   reverse_pairs)

        partners = set(p[1] for p in forward_pairs)
        used_reference = set()
        used_partner = set()
        mapping = {}
        out_path = os.path.join(genome_dir, reference + '_vs_' + other + '.txt')
        with open(out_path, 'w') as handle:
            for ref_gene, partner, _score in sorted(reverse_pairs, key=lambda t: -t[2]):
                if (partner in partners
                        and ref_gene not in used_reference
                        and partner not in used_partner):
                    used_reference.add(ref_gene)
                    used_partner.add(partner)
                    mapping[ref_gene] = partner
                    handle.write(ref_gene + ':' + partner + '\n')
        # kept in memory rather than parsed back out of the file just written
        pair_maps[other] = mapping
        LOG.info('%s vs %s: %d reciprocal pairs', reference, other, len(mapping))
    return pair_maps


# --------------------------------------------------------------------------------------
# core gene selection
# --------------------------------------------------------------------------------------

def select_core_genes(pair_maps):
    """Reference genes present in every pairwise reciprocal set, sorted for reproducibility."""
    if not pair_maps:
        raise SystemExit('No reciprocal comparisons were produced.')
    core = sorted(set.intersection(*[set(m) for m in pair_maps.values()]))
    if not core:
        raise SystemExit('No core genes are shared by all the genomes.')
    LOG.info('core genes: %d', len(core))
    return core


def write_core_gene_tables(core_dir, homolog_seq_dir, pair_maps, core_genes, reference):
    """Write per-genome core gene name tables and their FASTA.

    Every table is written in the same core_genes order, so row i is the same ortholog group
    in every organism -- the index alignment the original README promises.
    """
    names_dir = os.path.join(core_dir, 'gene_names')
    seq_dir = os.path.join(core_dir, 'gene_sequences')
    for path in (names_dir, seq_dir, os.path.join(core_dir, 'msa_seq')):
        os.makedirs(path, exist_ok=True)

    write_id_list(os.path.join(names_dir, reference), core_genes)
    for other, mapping in sorted(pair_maps.items()):
        write_id_list(os.path.join(names_dir, other), [mapping[g] for g in core_genes])

    for stem in sorted(os.listdir(names_dir)):
        identifiers = read_id_list(os.path.join(names_dir, stem))
        index = index_fasta(os.path.join(homolog_seq_dir, stem + '.fasta'))
        records = copy_records(index, identifiers, stem)
        SeqIO.write(records, os.path.join(seq_dir, stem + '.fasta'), 'fasta')
    return names_dir, seq_dir


# --------------------------------------------------------------------------------------
# stages H / I / J - supermatrix, alignment, tree (reused for both trees)
# --------------------------------------------------------------------------------------

def build_supermatrix(seq_dir, per_genome_dir, combined_path):
    """Concatenate each genome's genes into one sequence and collect them into one FASTA.

    Gene order inside a genome comes from its core gene table, which is identical across
    genomes, so column i of the supermatrix is the same ortholog group for every organism.
    """
    os.makedirs(per_genome_dir, exist_ok=True)
    labels = []
    with open(combined_path, 'w') as combined:
        for name in sorted(os.listdir(seq_dir)):
            joined = ''.join(str(r.seq) for r in SeqIO.parse(os.path.join(seq_dir, name), 'fasta'))
            label = stem_of(name)
            labels.append(label)
            with open(os.path.join(per_genome_dir, name), 'w') as handle:
                handle.write('>' + label + '\n' + joined + '\n')
            combined.write('>' + label + '\n' + joined + '\n')
    LOG.info('supermatrix: %d genomes -> %s', len(labels), combined_path)
    return combined_path


def run_alignment(tools, in_path, out_path, log_path, check=True):
    """Align a FASTA with ClustalW using the original's parameter set."""
    run_tool([tools['clustalw'],
              '-INFILE=' + in_path, '-QUICKTREE', '-OUTFILE=' + out_path, '-OUTPUT=FASTA',
              '-KTUPLE=2', '-TOPDIAGS=5', '-WINDOW=5', '-PAIRGAP=3',
              '-SCORE=PERCENT', '-GAPOPEN=8.0'],
             stdout_path=log_path, check=check)
    return out_path


def build_tree(tools, alignment_path, tree_path, description):
    """Infer a Newick tree from an alignment with FastTree."""
    require_nonempty(alignment_path, description + ' alignment')
    run_tool([tools['fasttree']], stdin_path=alignment_path, stdout_path=tree_path)
    require_nonempty(tree_path, description + ' tree')
    LOG.info('%s tree: %s', description, tree_path)
    return tree_path


# --------------------------------------------------------------------------------------
# stage K - one-line FASTA and ortholog group transpose
# --------------------------------------------------------------------------------------

def write_single_line_fasta(seq_dir, out_dir):
    """Rewrite each genome's core genes as strict two-line records."""
    os.makedirs(out_dir, exist_ok=True)
    for name in sorted(os.listdir(seq_dir)):
        with open(os.path.join(out_dir, name), 'w') as handle:
            for record in SeqIO.parse(os.path.join(seq_dir, name), 'fasta'):
                handle.write('>' + record.id + '\n' + str(record.seq) + '\n')
    return out_dir


def transpose_gene_sets(single_line_dir, out_dir, gene_count):
    """Turn per-genome files into per-ortholog-group files.

    rec<N>.fasta holds ortholog group N taken from every organism, which is what PhiPack
    needs. Relies on every genome's file listing the groups in the same order.
    """
    os.makedirs(out_dir, exist_ok=True)
    for name in sorted(os.listdir(single_line_dir)):
        with open(os.path.join(single_line_dir, name), 'r') as handle:
            lines = handle.readlines()
        if len(lines) != 2 * gene_count:
            raise SystemExit('%s holds %d records, expected %d; the gene sets would not line up.'
                             % (name, len(lines) // 2, gene_count))
        for position in range(0, len(lines), 2):
            group = os.path.join(out_dir, 'rec%d.fasta' % (position // 2))
            with open(group, 'a') as handle:
                handle.write(lines[position] + lines[position + 1])
    return out_dir


# --------------------------------------------------------------------------------------
# stages L / M - per-group alignment and recombination test
# --------------------------------------------------------------------------------------

def align_gene_sets(tools, gene_set_dir, aligned_dir, log_dir):
    """Align every ortholog group. A single failed group is tolerated and reported later."""
    for path in (aligned_dir, log_dir):
        os.makedirs(path, exist_ok=True)
    names = sorted(os.listdir(gene_set_dir))
    for name in names:
        base = stem_of(name)
        run_alignment(tools,
                      os.path.join(gene_set_dir, name),
                      os.path.join(aligned_dir, base + '.aln.fasta'),
                      os.path.join(log_dir, base + '.log'),
                      check=False)
    LOG.info('aligned %d ortholog groups', len(names))
    return aligned_dir


def run_recombination_tests(tools, aligned_dir, phi_dir):
    """Run PhiPack on every aligned ortholog group.

    A non-zero exit is tolerated: PhiPack legitimately declines groups with too few
    informative sites. Those groups are reported by parse_phi_values.
    """
    os.makedirs(phi_dir, exist_ok=True)
    for name in sorted(os.listdir(aligned_dir)):
        run_tool([tools['phi'], '-f', os.path.join(aligned_dir, name),
                  '-t', 'A', '-w', '10', '-v', '-g', 'i'],
                 stdout_path=os.path.join(phi_dir, name.split('.aln.fasta')[0]),
                 check=False)
    return phi_dir


# --------------------------------------------------------------------------------------
# stage N - PHI values and the recombinant / non-recombinant split
# --------------------------------------------------------------------------------------

def parse_phi_values(phi_dir, gene_count):
    """Read the 'PHI (Normal):' p-value of every ortholog group.

    Every group is seeded as untested ('--') so a group whose PhiPack or ClustalW run produced
    no PHI line still appears in the report instead of vanishing from it. Untested groups are
    treated as non-recombinant -- absence of evidence -- and their number is returned so the
    caller can say so out loud.
    """
    values = {'rec%d' % i: '--' for i in range(gene_count)}
    tested = set()
    for name in sorted(os.listdir(phi_dir)):
        if name not in values:
            continue
        with open(os.path.join(phi_dir, name), 'r') as handle:
            for line in handle:
                if 'PHI (Normal):' in line:
                    raw = line.rstrip('\n').split(':')[1].strip()
                    try:
                        values[name] = float(raw)
                    except ValueError:
                        values[name] = '--'
                    tested.add(name)
                    break
    untested = sorted(set(values) - tested, key=lambda k: int(k[3:]))
    if untested:
        LOG.warning('%d of %d ortholog groups produced no PHI value (too few informative sites, '
                    'or the PhiPack/ClustalW run failed); they are reported as "--" and kept in '
                    'the non-recombinant set', len(untested), gene_count)
    return values, untested


def split_recombinant_genes(values, alpha, core_dir, nonrec_dir, rec_data_dir):
    """Partition each genome's core gene table into recombinant and non-recombinant halves.

    Indexing is by the group numbers PhiPack actually flagged. The original deleted list
    positions 0,1,2... from a shrinking list, so it removed the wrong genes entirely.
    """
    for path in (nonrec_dir, rec_data_dir):
        os.makedirs(path, exist_ok=True)

    order = sorted(values, key=lambda k: int(k[3:]))
    with open(os.path.join(rec_data_dir, 'sorted_List_Phi_Values.txt'), 'w') as handle:
        for key in order:
            handle.write('%s:%s\n' % (key, values[key]))

    flagged = []
    with open(os.path.join(rec_data_dir, 'all_core_genes_having_recombination.txt'), 'w') as handle:
        for key in order:
            value = values[key]
            if value != '--' and float(value) < alpha:
                flagged.append(int(key[3:]))
                handle.write('%s:%s\n' % (key, value))
    flagged_set = set(flagged)
    LOG.info('recombinant core genes: %d of %d', len(flagged), len(values))

    names_dir = os.path.join(core_dir, 'gene_names')
    rec_tables = os.path.join(rec_data_dir, 'All-organisms-Reco-data')
    nonrec_tables = os.path.join(nonrec_dir, 'All-organisms-nReco-data')
    for path in (rec_tables, nonrec_tables):
        os.makedirs(path, exist_ok=True)
    for stem in sorted(os.listdir(names_dir)):
        identifiers = read_id_list(os.path.join(names_dir, stem))
        write_id_list(os.path.join(rec_tables, stem),
                      [g for i, g in enumerate(identifiers) if i in flagged_set])
        write_id_list(os.path.join(nonrec_tables, stem),
                      [g for i, g in enumerate(identifiers) if i not in flagged_set])
    if len(flagged) == len(values):
        raise SystemExit('Every core gene was flagged as recombinant, so no non-recombinant tree '
                         'can be built. Consider a stricter --phi-alpha.')
    return nonrec_tables, flagged


def write_nonrecombinant_sequences(single_line_dir, table_dir, out_dir):
    """Collect each genome's non-recombinant core gene sequences."""
    os.makedirs(out_dir, exist_ok=True)
    for name in sorted(os.listdir(single_line_dir)):
        table = os.path.join(table_dir, stem_of(name))
        if not os.path.isfile(table):
            continue
        identifiers = read_id_list(table)
        index = index_fasta(os.path.join(single_line_dir, name))
        records = copy_records(index, identifiers, stem_of(name))
        SeqIO.write(records, os.path.join(out_dir, name), 'fasta')
    return out_dir


# --------------------------------------------------------------------------------------
# orchestration
# --------------------------------------------------------------------------------------

def run_pipeline(settings, tools):
    """Drive every stage in order and return the two tree paths.

    `settings` is a Settings; `tools` is the dict returned by resolve_tools().
    """
    out = settings.output_dir
    prepare_output_directory(out, settings.force)
    work_dir = os.path.join(out, 'working_sequences')

    LOG.info('[A] preparing working copy')
    prepare_working_copy(settings.input_dir, work_dir)
    LOG.info('[B] normalising input file names')
    files = normalise_input_filenames(work_dir)

    LOG.info('[C] all-vs-all BLASTP (%d genomes, %d threads)', len(files), settings.threads)
    blast_dir = os.path.join(out, 'first_blast')
    stems = run_all_vs_all_blast(tools, work_dir, files, blast_dir,
                                 os.path.join(out, 'Blast_log', 'BLAST-logfile'), settings.threads)

    LOG.info('[D] extracting homologs')
    homolog_dir = os.path.join(out, 'homogenome')
    counts = extract_homologs(blast_dir, homolog_dir, stems,
                              settings.pident, settings.gapopen, settings.evalue)

    LOG.info('[E] selecting the reference genome')
    reference = select_reference_genome(counts)

    LOG.info('[F] writing homolog sequences')
    homolog_seq_dir = write_homolog_sequences(work_dir, homolog_dir,
                                              os.path.join(out, 'homogenome_seq'), stems)

    LOG.info('[G] finding reciprocal pairs')
    pair_maps = find_reciprocal_pairs(blast_dir, homolog_dir, os.path.join(out, 'Reciprocal'),
                                      stems, reference, settings.pident, settings.gapopen, settings.evalue)
    core_genes = select_core_genes(pair_maps)
    core_dir = os.path.join(out, 'CORE_Genes')
    _, core_seq_dir = write_core_gene_tables(core_dir, homolog_seq_dir, pair_maps,
                                             core_genes, reference)

    LOG.info('[H] building the core-gene supermatrix')
    clustal_dir = os.path.join(out, 'CLUSTALW')
    os.makedirs(clustal_dir, exist_ok=True)
    core_matrix = build_supermatrix(core_seq_dir, os.path.join(core_dir, 'msa_seq'),
                                    os.path.join(clustal_dir, 'all_seq_for_msa.fasta'))
    LOG.info('[I] aligning the core-gene supermatrix')
    core_aln = run_alignment(tools, core_matrix,
                             os.path.join(clustal_dir, 'all_seq_for_msa.aln.fasta'),
                             os.path.join(clustal_dir, 'CLUSTALW-logfile'))
    LOG.info('[J] building the core-gene tree')
    core_tree = build_tree(tools, core_aln,
                           os.path.join(out, 'FAST-TREE', 'core_genes.nwk'), 'core-gene')

    LOG.info('[K] transposing ortholog groups')
    nonrec_dir = os.path.join(out, 'Non-Recombinant')
    single_line_dir = write_single_line_fasta(core_seq_dir,
                                              os.path.join(nonrec_dir, 'Single-Liners_FASTA'))
    gene_sets = transpose_gene_sets(single_line_dir,
                                    os.path.join(nonrec_dir, 'Gene_set_FASTA'), len(core_genes))

    LOG.info('[L] aligning each ortholog group')
    aligned = align_gene_sets(tools, gene_sets,
                              os.path.join(nonrec_dir, 'Gene_set_FASTA_aligned'),
                              os.path.join(nonrec_dir, 'Log-Files'))
    LOG.info('[M] testing each ortholog group for recombination')
    phi_dir = run_recombination_tests(tools, aligned, os.path.join(nonrec_dir, 'Phipack'))

    LOG.info('[N] splitting recombinant from non-recombinant genes')
    values, _untested = parse_phi_values(phi_dir, len(core_genes))
    table_dir, _flagged = split_recombinant_genes(
        values, settings.phi_alpha, core_dir, nonrec_dir,
        os.path.join(nonrec_dir, 'Recombinant-data'))
    nonrec_seq_dir = write_nonrecombinant_sequences(
        single_line_dir, table_dir, os.path.join(nonrec_dir, 'All-organisms-nReco-Sequences'))

    nr_clustal_dir = os.path.join(nonrec_dir, 'nR-CLUSTALW')
    os.makedirs(nr_clustal_dir, exist_ok=True)
    nonrec_matrix = build_supermatrix(nonrec_seq_dir,
                                      os.path.join(nonrec_dir, 'nReco-msa_seq'),
                                      os.path.join(nr_clustal_dir, 'nR_all_seq_for_msa.fasta'))
    LOG.info('[O] aligning the non-recombinant supermatrix')
    nonrec_aln = run_alignment(tools, nonrec_matrix,
                               os.path.join(nr_clustal_dir, 'nR_all_seq_for_msa.aln.fasta'),
                               os.path.join(nr_clustal_dir, 'CLUSTALW-logfile'))
    LOG.info('[P] building the non-recombinant tree')
    nonrec_tree = build_tree(tools, nonrec_aln,
                             os.path.join(nonrec_dir, 'FAST-TREE', 'nReco_genes.nwk'),
                             'non-recombinant')

    if settings.clean:
        LOG.info('removing the working copy (--clean)')
        shutil.rmtree(work_dir, ignore_errors=True)
    return core_tree, nonrec_tree


# --------------------------------------------------------------------------------------
# command line
# --------------------------------------------------------------------------------------

def build_parser():
    """Define the command-line interface."""
    parser = argparse.ArgumentParser(
        prog='pygemini',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='Core-genome and recombination-aware phylogeny from protein FASTA files.',
        epilog='example:\n'
               '  pygemini -i proteomes/ -o results/ -t 8\n\n'
               'Requires makeblastdb, blastp, clustalw, fasttree and Phi on PATH:\n'
               '  micromamba env create -f environment.yml && micromamba activate pygemini\n')
    parser.add_argument('-i', '--input', required=True, metavar='DIR',
                        help='directory of protein FASTA files, one per genome (never modified)')
    parser.add_argument('-o', '--output', metavar='DIR',
                        help='directory for results (default: <input>_Result)')
    parser.add_argument('-t', '--threads', type=int, default=1, metavar='N',
                        help='threads passed to blastp (default: 1, this machine has %d)'
                             % (os.cpu_count() or 1))
    group = parser.add_argument_group('analysis thresholds')
    group.add_argument('--pident', type=float, default=85.0, metavar='PCT',
                       help='minimum percentage identity for a BLAST hit (default: 85)')
    group.add_argument('--gapopen', type=int, default=3, metavar='N',
                       help='maximum gap opens in a BLAST hit (default: 3)')
    group.add_argument('--evalue', type=float, default=1e-4, metavar='E',
                       help='maximum BLAST e-value (default: 1e-4). The original script wrote '
                            '"10e-5", which is 1e-4 rather than the 1e-5 it resembles; the '
                            'original value is kept as the default so results stay comparable')
    group.add_argument('--phi-alpha', type=float, default=0.05, metavar='P',
                       help='PHI p-value below which a gene counts as recombinant (default: 0.05)')
    behaviour = parser.add_argument_group('behaviour')
    behaviour.add_argument('--force', action='store_true',
                           help='delete and recreate the output directory if it is not empty')
    behaviour.add_argument('--clean', action='store_true',
                           help='delete the working copy of the sequences when finished')
    behaviour.add_argument('--log-file', metavar='FILE', help='also write the log to FILE')
    verbosity = behaviour.add_mutually_exclusive_group()
    verbosity.add_argument('-v', '--verbose', action='store_true', help='show every command run')
    verbosity.add_argument('-q', '--quiet', action='store_true', help='warnings and errors only')
    parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)
    return parser


def configure_logging(verbose, quiet):
    """Send progress to stderr. Handlers are reset so a second call cannot double every line."""
    # the logger itself stays wide open and each handler filters, so -q quietens the console
    # without also emptying a --log-file
    LOG.setLevel(logging.DEBUG)
    LOG.handlers.clear()
    stream = logging.StreamHandler()
    stream.setLevel(logging.DEBUG if verbose else logging.WARNING if quiet else logging.INFO)
    stream.setFormatter(logging.Formatter('%(levelname)-8s %(message)s'))
    LOG.addHandler(stream)


def add_log_file(log_file):
    """Tee the log to a file.

    Called only after the arguments validate, so that creating the log cannot make the output
    directory look non-empty, nor hold a handle open across the --force rmtree.
    """
    os.makedirs(os.path.dirname(os.path.abspath(log_file)), exist_ok=True)
    handler = logging.FileHandler(log_file, mode='w')
    handler.setLevel(logging.DEBUG)   # the file always gets the full run, whatever -q/-v say
    handler.setFormatter(logging.Formatter('%(asctime)s %(levelname)-8s %(message)s'))
    LOG.addHandler(handler)


class Settings(NamedTuple):
    """Everything the pipeline needs, fixed once at start-up and never mutated afterwards.

    This is the whole contract between the command line and run_pipeline. Presentation-only
    options (--verbose, --quiet, --log-file) are deliberately absent: they belong to main(),
    not to the analysis.
    """
    input_dir: str        # directory of protein FASTA files; read-only
    output_dir: str       # directory every result is written under
    threads: int          # passed to blastp -num_threads
    pident: float         # minimum percentage identity for a BLAST hit
    gapopen: int          # maximum gap opens in a BLAST hit
    evalue: float         # maximum BLAST e-value
    phi_alpha: float      # PHI p-value below which a gene counts as recombinant
    force: bool           # delete a non-empty output directory instead of refusing
    clean: bool           # remove the working copy when finished


def settings_from_args(args):
    """Turn the parsed command line into validated Settings.

    Rejects any combination that would destroy data, and resolves both directories with
    realpath so that a symlink cannot slip past the string containment tests below.
    """
    input_dir = os.path.realpath(os.path.expanduser(args.input))
    if not os.path.isdir(input_dir):
        raise SystemExit('Not a directory: ' + input_dir)
    output = args.output if args.output is not None else input_dir.rstrip(os.sep) + '_Result'
    output_dir = os.path.realpath(os.path.expanduser(output))

    if os.path.normcase(output_dir) == os.path.normcase(input_dir):
        raise SystemExit('--output must not be the input directory: ' + output_dir)
    if path_contains(input_dir, output_dir):
        raise SystemExit('--output must not sit inside the input directory: ' + output_dir)
    # the reverse matters just as much: --force deletes the whole output tree, so an output
    # directory that CONTAINS the input would destroy the proteomes it is about to read
    if path_contains(output_dir, input_dir):
        raise SystemExit('--output must not contain the input directory: ' + output_dir +
                         '\nIt would be deleted with --force, taking your input with it.')

    if args.threads < 1:
        raise SystemExit('--threads must be 1 or more.')
    if not 0 < args.pident <= 100:
        raise SystemExit('--pident must be between 0 and 100.')
    if args.gapopen < 0:
        raise SystemExit('--gapopen must be 0 or more.')
    if args.evalue <= 0:
        raise SystemExit('--evalue must be greater than 0.')
    if not 0 < args.phi_alpha < 1:
        raise SystemExit('--phi-alpha must be between 0 and 1.')

    return Settings(input_dir=input_dir, output_dir=output_dir, threads=args.threads,
                    pident=args.pident, gapopen=args.gapopen, evalue=args.evalue,
                    phi_alpha=args.phi_alpha, force=args.force, clean=args.clean)


def validated_log_path(log_file, output_dir):
    """Resolve --log-file and keep it out of the output tree.

    A log written inside --output would make the directory look non-empty, and under --force
    would hold an open handle across the rmtree.
    """
    path = os.path.realpath(os.path.expanduser(log_file))
    if path_contains(output_dir, path):
        raise SystemExit('--log-file must not sit inside the output directory: ' + path)
    return path


def main(argv=None):
    """Entry point. Returns a process exit status."""
    args = build_parser().parse_args(argv)
    configure_logging(args.verbose, args.quiet)
    try:
        settings = settings_from_args(args)
        if args.log_file:
            add_log_file(validated_log_path(args.log_file, settings.output_dir))
        tools = resolve_tools()
        core_tree, nonrec_tree = run_pipeline(settings, tools)
    except SystemExit as exc:
        if isinstance(exc.code, str):
            LOG.error('%s', exc.code)
            return 1
        raise
    except KeyboardInterrupt:
        LOG.error('interrupted')
        return 130
    print('\nCore gene tree      : ' + core_tree)
    print('Non-recombinant tree: ' + nonrec_tree)
    return 0


if __name__ == '__main__':
    sys.exit(main())
