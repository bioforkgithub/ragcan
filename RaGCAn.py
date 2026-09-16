#!/usr/bin/env python3
"""RaGCAn (Rapid Genome Coherence Analyzer) - rapid core-genome screen for taxonomic coherence.

Answers one question quickly: do a set of genomes share a core genome at all, and if not,
which genomes are responsible for the collapse?

A taxon whose members genuinely belong together retains a substantial core. When a strain does
not belong - a novel organism, or one classified into the wrong taxon - the shared core falls
away. That collapse is the signal this program measures, and it needs no alignment and no tree,
so it runs in minutes rather than days.

What it does NOT do, by design: multiple sequence alignment, recombination testing, and tree
inference. Those are left to tools built for them (MAFFT, IQ-TREE, RAxML-NG). The core gene
sequences are written out per genome and per ortholog group so they can be taken straight into
such a pipeline. The earlier full pipeline is preserved as pygemini_full.py.

Method
    1  every proteome is pooled into ONE DIAMOND database, searched all-against-all in a
       single pass                                          run_all_vs_all_search
    2  best cross-genome hits are recorded                  parse_best_hits
    3  reciprocal best hits define orthology                reciprocal_best_hits
    4  the genome with the FEWEST homologs anchors the core select_reference_genome
    5  core = reference genes with a partner everywhere     core_proteins
    6  genomes are binned into candidate sub-taxa           bin_genomes
    7  leave-one-out and greedy removal re-apply the anchor rule inside each subset

The anchoring rule is the author's, and it is the idea the program is built on. The core genome
is by definition contained in every genome's gene set, so the genome sharing the fewest genes
gives the TIGHTEST superset of the core: the most efficient pivot for the intersection, and the
most conservative one. Requiring a direct reciprocal best hit from that genome to every other is
also stricter than joining orthologs transitively through a graph, which can chain through
errors. Because the rule is re-applied within every subset, excluding the anchor genome itself
is handled naturally rather than as a special case.

Every stage is a single function so that each can be checked, or reused, on its own.
"""

import argparse
import logging
import os
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from typing import NamedTuple

import numpy as np

__version__ = '3.0.0'

LOG = logging.getLogger('RaGCAn')

FASTA_SUFFIXES = ('.faa', '.fa', '.fasta', '.fas', '.pep', '.protein')


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
    """Locate DIAMOND. It is the only external program this version needs."""
    diamond = find_executable('diamond')
    if diamond is None:
        raise SystemExit(
            'DIAMOND was not found on PATH.\n'
            'Install it into the active environment:\n'
            '    micromamba install -c conda-forge -c bioconda diamond\n'
            'or create the bundled environment:\n'
            '    micromamba env create -f environment.yml')
    LOG.debug('using diamond: %s', diamond)
    return {'diamond': diamond}


def run_tool(cmd, stdout_path=None, check=True):
    """Run an external program, aborting on a non-zero exit unless check=False."""
    LOG.debug('run: %s', ' '.join(cmd))
    handle = None
    if stdout_path:
        os.makedirs(os.path.dirname(os.path.abspath(stdout_path)), exist_ok=True)
    try:
        handle = open(stdout_path, 'w') if stdout_path else None
        code = subprocess.call(cmd, stdout=handle)
    finally:
        if handle:
            handle.close()
    if check and code != 0:
        raise SystemExit('%s failed (exit status %d):\n    %s'
                         % (os.path.basename(cmd[0]), code, ' '.join(cmd)))
    return code


# --------------------------------------------------------------------------------------
# small shared helpers
# --------------------------------------------------------------------------------------

def stem_of(name):
    """Filename without its extension. Never use str.rstrip for this."""
    return os.path.splitext(name)[0]


def path_contains(parent, child):
    """True if `child` is `parent` itself or lies anywhere inside it."""
    parent = os.path.normcase(os.path.normpath(parent))
    child = os.path.normcase(os.path.normpath(child))
    return child == parent or child.startswith(parent.rstrip(os.sep) + os.sep)


def iter_fasta(path):
    """Yield (identifier, description, sequence) from a FASTA file.

    Deliberately not Biopython: this is the only parsing the program does, it is trivial,
    and dropping the dependency keeps start-up instant.
    """
    identifier = description = None
    chunks = []
    with open(path, 'r') as handle:
        for line in handle:
            line = line.rstrip('\n').rstrip('\r')
            if line.startswith('>'):
                if identifier is not None:
                    yield identifier, description, ''.join(chunks)
                header = line[1:].strip()
                parts = header.split(None, 1)
                identifier = parts[0] if parts else ''
                description = parts[1] if len(parts) > 1 else ''
                chunks = []
            elif line:
                chunks.append(line.strip())
    if identifier is not None:
        yield identifier, description, ''.join(chunks)


def human_count(value):
    """Thousands-separated integer, for log lines."""
    return format(int(value), ',d')


# --------------------------------------------------------------------------------------
# proteome loading
# --------------------------------------------------------------------------------------

class Proteomes(NamedTuple):
    """Every protein from every genome, flattened and integer-indexed.

    Integer indexing is what makes the rest of the program fast: the DIAMOND database is
    written with the serial number as each sequence's identifier, so parsing its output is
    two int() calls per line with no string lookups at all.
    """
    genomes: list        # genome names, in sorted order; index = genome id
    names: list          # names[i]       = original identifier of protein i
    descriptions: list   # descriptions[i]= original description of protein i
    sequences: list      # sequences[i]   = amino acid sequence of protein i
    genome_of: object    # int32 array, genome_of[i] = genome id of protein i
    offsets: list        # offsets[g] = (first protein index, one past last) for genome g


def load_proteomes(input_dir):
    """Read every FASTA in `input_dir` and flatten it into integer-indexed arrays."""
    files = sorted(f for f in os.listdir(input_dir)
                   if os.path.isfile(os.path.join(input_dir, f))
                   and f.lower().endswith(FASTA_SUFFIXES))
    if len(files) < 2:
        raise SystemExit('At least two protein FASTA files are needed, found %d in %s\n'
                         'Recognised extensions: %s'
                         % (len(files), input_dir, ', '.join(FASTA_SUFFIXES)))
    stems = [stem_of(f) for f in files]
    duplicates = sorted({s for s in stems if stems.count(s) > 1})
    if duplicates:
        raise SystemExit('These input names collide once the extension is removed: '
                         + ', '.join(duplicates))

    names, descriptions, sequences, genome_of, offsets = [], [], [], [], []
    for genome_id, filename in enumerate(files):
        start = len(names)
        for identifier, description, sequence in iter_fasta(os.path.join(input_dir, filename)):
            if not sequence:
                continue
            names.append(identifier)
            descriptions.append(description)
            sequences.append(sequence)
            genome_of.append(genome_id)
        offsets.append((start, len(names)))
        if len(names) == start:
            raise SystemExit('No sequences found in ' + filename)
        LOG.info('%s: %s proteins', stems[genome_id], human_count(len(names) - start))
    LOG.info('%d genomes, %s proteins total', len(files), human_count(len(names)))
    return Proteomes(genomes=stems, names=names, descriptions=descriptions, sequences=sequences,
                     genome_of=np.array(genome_of, dtype=np.int32), offsets=offsets)


def write_pooled_fasta(proteomes, path):
    """Write every protein into one FASTA whose identifiers are the serial indices."""
    with open(path, 'w') as handle:
        for index, sequence in enumerate(proteomes.sequences):
            handle.write('>%d\n%s\n' % (index, sequence))
    return path


# --------------------------------------------------------------------------------------
# the single all-against-all search
# --------------------------------------------------------------------------------------

def run_all_vs_all_search(tools, pooled_fasta, work_dir, settings):
    """Build one DIAMOND database and search it against itself in a single pass.

    The original ran makeblastdb plus blastp once per ordered genome pair - N(N-1) searches,
    9,900 of them at 100 genomes. One pooled search replaces all of it.

    The search dominates the runtime and everything after it takes seconds, so a finished
    search is recorded with a marker file and reused on a later run. DIAMOND writes its output
    progressively and cannot resume, so a search that was interrupted has no marker and is
    redone from the start.
    """
    database = os.path.join(work_dir, 'pooled')
    hits_path = os.path.join(work_dir, 'all_vs_all.tsv')
    marker = hits_path + '.complete'
    if os.path.exists(marker) and os.path.getsize(hits_path) > 0 and not settings.rerun_search:
        with open(marker) as handle:
            recorded = handle.read().strip()
        LOG.info('reusing the finished search at %s (%.2f GB)',
                 hits_path, os.path.getsize(hits_path) / 1e9)
        LOG.info('  it was run as: %s', recorded)
        LOG.info('  pass --rerun-search to discard it and search again')
        return hits_path
    for stale in (marker, hits_path):
        if os.path.exists(stale):
            os.remove(stale)
    run_tool([tools['diamond'], 'makedb', '--in', pooled_fasta, '-d', database,
              '--threads', str(settings.threads), '--quiet'])
    command = [tools['diamond'], 'blastp',
               '-q', pooled_fasta, '-d', database, '-o', hits_path,
               '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'length',
               'qlen', 'slen', 'gapopen', 'evalue', 'bitscore',
               '--threads', str(settings.threads),
               '--evalue', repr(settings.evalue),
               '--max-target-seqs', str(settings.max_target_seqs),
               # DIAMOND's peak memory is roughly 6 x block-size GB. The default of 2.0 wants
               # about 12 GB, which will thrash or be killed on an ordinary machine.
               '--block-size', repr(settings.block_size),
               '--index-chunks', str(settings.index_chunks),
               '--quiet']
    if settings.sensitivity:
        command.append('--' + settings.sensitivity)
    LOG.info('running DIAMOND all-against-all (%s mode, %d threads)',
             settings.sensitivity or 'default', settings.threads)
    LOG.info('this is the slow step; everything after it takes seconds')
    run_tool(command)
    if os.path.getsize(hits_path) == 0:
        raise SystemExit('The search produced no output. It was probably interrupted; nothing '
                         'downstream can run.')
    # only now is the file trustworthy
    with open(marker, 'w') as handle:
        handle.write('%s mode, evalue %g, max-target-seqs %d'
                     % (settings.sensitivity or 'default', settings.evalue,
                        settings.max_target_seqs))
    LOG.info('search finished: %.2f GB of hits', os.path.getsize(hits_path) / 1e9)
    return hits_path


# --------------------------------------------------------------------------------------
# hit parsing
# --------------------------------------------------------------------------------------

class HitTables(NamedTuple):
    """Cross-genome hit summaries, one row per protein and one column per genome."""
    best_hit: object     # int32  [protein, genome] -> best-scoring subject there, or -1
    best_pident: object  # float32[protein, genome] -> that subject's percentage identity
    hit_from: object     # bool   [protein, genome] -> some protein there hit this one


def parse_best_hits(hits_path, proteomes, settings):
    """Stream the DIAMOND output once and record, for every protein, its best hit in each
    other genome and which genomes hit it.

    The identity threshold is deliberately NOT applied here. AAI is computed downstream as the
    mean identity over reciprocal best hits, and pre-filtering at, say, 85% would discard every
    pair below that and report an AAI biased upwards towards the cutoff. Identity filtering
    belongs to core-gene selection alone, and is applied there.

    Same-genome hits are skipped: only cross-genome relationships matter, and that also removes
    every self-hit without a special case.
    """
    protein_count = len(proteomes.names)
    genome_count = len(proteomes.genomes)
    best_hit = np.full((protein_count, genome_count), -1, dtype=np.int32)
    best_score = np.zeros((protein_count, genome_count), dtype=np.float32)
    best_pident = np.zeros((protein_count, genome_count), dtype=np.float32)
    hit_from = np.zeros((protein_count, genome_count), dtype=bool)
    genome_of = proteomes.genome_of

    kept = examined = 0
    with open(hits_path, 'r') as handle:
        for line in handle:
            fields = line.split('\t')
            if len(fields) < 9:
                continue
            examined += 1
            query = int(fields[0])
            subject = int(fields[1])
            query_genome = genome_of[query]
            subject_genome = genome_of[subject]
            if query_genome == subject_genome:
                continue
            if int(fields[6]) > settings.gapopen:
                continue
            if settings.min_coverage > 0.0:
                length = int(fields[3])
                if (length / int(fields[4]) < settings.min_coverage
                        or length / int(fields[5]) < settings.min_coverage):
                    continue
            if float(fields[7]) >= settings.evalue:
                continue
            kept += 1
            hit_from[subject, query_genome] = True
            score = float(fields[8])
            if score > best_score[query, subject_genome]:
                best_score[query, subject_genome] = score
                best_pident[query, subject_genome] = float(fields[2])
                best_hit[query, subject_genome] = subject
    LOG.info('hits: %s examined, %s passed the filters', human_count(examined), human_count(kept))
    LOG.debug('hit tables: %.2f GB', (best_hit.nbytes + best_pident.nbytes
                                      + hit_from.nbytes) / 1e9)
    if kept == 0:
        raise SystemExit('No hits passed the filters. Relax --evalue, raise --sensitivity, '
                         'or check that the inputs are protein FASTA.')
    del best_score          # only needed while choosing the best hit; ~0.3 GB at 138 genomes
    return HitTables(best_hit=best_hit, best_pident=best_pident, hit_from=hit_from)


# --------------------------------------------------------------------------------------
# orthology
# --------------------------------------------------------------------------------------

def reciprocal_best_hits(proteomes, hits, pident_min):
    """Find every reciprocal best hit.

    Returns (partner, left, right, identity):

      partner   int32 [protein, genome] -> the protein's reciprocal best hit there, or -1.
                This is the structure the core genome is computed from.
      left/right/identity
                the same pairs flattened, used for AAI. These are NOT identity-filtered,
                because filtering them would bias AAI upwards; `partner` is.

    Vectorised per genome pair, so the cost is one array operation per pair rather than one
    Python-level test per protein per genome.
    """
    best_hit = hits.best_hit
    best_pident = hits.best_pident
    genome_count = len(proteomes.genomes)
    partner = np.full((len(proteomes.names), genome_count), -1, dtype=np.int32)
    left, right, identity = [], [], []
    for a in range(genome_count):
        start, stop = proteomes.offsets[a]
        proteins_a = np.arange(start, stop, dtype=np.int32)
        for b in range(a + 1, genome_count):
            forward = best_hit[proteins_a, b]
            has_forward = forward >= 0
            if not has_forward.any():
                continue
            sources = proteins_a[has_forward]
            targets = forward[has_forward]
            mutual = best_hit[targets, a] == sources
            if not mutual.any():
                continue
            sources = sources[mutual]
            targets = targets[mutual]
            # the two directions can differ slightly; the mean is the conventional choice
            pair_identity = (best_pident[sources, b] + best_pident[targets, a]) / 2.0
            left.append(sources)
            right.append(targets)
            identity.append(pair_identity)
            # orthology, unlike AAI, does obey the identity threshold
            strong = pair_identity >= pident_min
            partner[sources[strong], b] = targets[strong]
            partner[targets[strong], a] = sources[strong]
    if not left:
        raise SystemExit('No reciprocal best hits were found between any pair of genomes.')
    left = np.concatenate(left)
    right = np.concatenate(right)
    identity = np.concatenate(identity)
    kept = int((partner >= 0).sum() // 2)
    LOG.info('reciprocal best hits: %s pairs, %s at or above %g%% identity',
             human_count(len(left)), human_count(kept), pident_min)
    if kept == 0:
        raise SystemExit('No reciprocal best hit reached --pident %g%%. These genomes may be more '
                         'divergent than that threshold allows. Try a lower --pident (50-70 is '
                         'usual for cross-species work) and consult reports/aai_matrix.tsv.'
                         % pident_min)
    return partner, left, right, identity


# Conventional boundaries from the prokaryotic taxonomy literature. They are guidelines that
# describe where most described taxa fall, not laws, and are reported as such.
AAI_SPECIES = 95.0    # at or above: same species
AAI_GENUS = 65.0      # at or above: same genus; below: usually a different genus
AAI_DIVERGENT = 80.0  # between GENUS and this: same genus but clearly distinct species


def pairwise_identity_matrices(proteomes, left, right, identity, threads=1):
    """Return (AAI matrix, shared-ortholog-count matrix) over all genome pairs.

    AAI - average amino acid identity - is the mean percentage identity across the reciprocal
    best hits between two genomes. It is one of the standard measures used to decide whether
    two organisms belong to the same species or genus, and it costs nothing extra here because
    the identities were already parsed.

    This is the O(genomes^2) step: every pair does a full boolean-mask pass over the hit
    array. Each numpy comparison releases the GIL for its C-level pass, so spreading pairs
    across `threads` worker threads gives real wall-clock speedup even though this is pure
    Python-level parallelism. Threads are split by row `a`; every (a, b) cell with a < b is
    written by exactly one thread, so no locking is needed.
    """
    genome_count = len(proteomes.genomes)
    aai = np.full((genome_count, genome_count), np.nan, dtype=np.float64)
    shared = np.zeros((genome_count, genome_count), dtype=np.int64)
    genome_of = proteomes.genome_of
    genome_left = genome_of[left]
    genome_right = genome_of[right]
    for a in range(genome_count):
        aai[a, a] = 100.0

    def compute_row(a):
        for b in range(a + 1, genome_count):
            selected = ((genome_left == a) & (genome_right == b)) | \
                       ((genome_left == b) & (genome_right == a))
            count = int(selected.sum())
            shared[a, b] = shared[b, a] = count
            if count:
                mean = float(identity[selected].mean())
                aai[a, b] = aai[b, a] = mean

    workers = max(1, min(threads, genome_count))
    if workers > 1:
        with ThreadPoolExecutor(max_workers=workers) as pool:
            list(pool.map(compute_row, range(genome_count)))
    else:
        for a in range(genome_count):
            compute_row(a)
    return aai, shared


def classify_pair(aai_value):
    """Turn one AAI value into a conventional taxonomic reading."""
    if np.isnan(aai_value):
        return 'no orthologs', 'unrelated at this threshold, or the search missed them'
    if aai_value >= AAI_SPECIES:
        return 'same species', 'AAI >= %g%%' % AAI_SPECIES
    if aai_value >= AAI_DIVERGENT:
        return 'same genus, different species', '%g%% <= AAI < %g%%' % (AAI_DIVERGENT, AAI_SPECIES)
    if aai_value >= AAI_GENUS:
        return 'same genus, divergent', '%g%% <= AAI < %g%%' % (AAI_GENUS, AAI_DIVERGENT)
    return 'likely different genus', 'AAI < %g%%' % AAI_GENUS


def predict_speciation(proteomes, aai, shared, matrix, genome_ids, leave_out_rows):
    """Per genome, a prediction about where it sits taxonomically.

    Two independent lines of evidence are combined:

      * AAI to its nearest neighbour, read against the conventional boundaries above
      * how much the shared core grows when that genome is removed, which identifies a
        genome that is single-handedly collapsing the core

    Agreement between the two is what makes a claim worth pursuing. Either alone is weak, and
    the report says so rather than pretending otherwise.
    """
    names = proteomes.genomes
    baseline = len(core_group_indices(matrix, genome_ids))
    recovery = {name: size - baseline for name, size, _delta in
                [(n, s, d) for n, s, d in leave_out_rows]}
    best_gain = max(recovery.values()) if recovery else 0

    rows = []
    for genome_id in genome_ids:
        others = [g for g in genome_ids if g != genome_id]
        if not others:
            continue
        values = np.array([aai[genome_id, g] for g in others], dtype=np.float64)
        if np.all(np.isnan(values)):
            nearest, nearest_aai = '-', float('nan')
        else:
            position = int(np.nanargmax(values))
            nearest = names[others[position]]
            nearest_aai = float(values[position])
        reading, basis = classify_pair(nearest_aai)
        gain = recovery.get(names[genome_id], 0)
        collapses_core = bool(best_gain > 0 and gain == best_gain and gain >= max(1, baseline))

        if np.isnan(nearest_aai) or nearest_aai < AAI_GENUS:
            verdict = 'CANDIDATE NOVEL TAXON'
        elif collapses_core and nearest_aai < AAI_SPECIES:
            verdict = 'REVIEW - collapses the shared core'
        elif nearest_aai < AAI_DIVERGENT:
            verdict = 'distinct species within the genus'
        elif nearest_aai < AAI_SPECIES:
            verdict = 'distinct species'
        else:
            verdict = 'consistent with the group'
        rows.append((names[genome_id], nearest,
                     'NA' if np.isnan(nearest_aai) else round(nearest_aai, 2),
                     int(shared[genome_id, others[int(np.nanargmax(values))]])
                     if not np.all(np.isnan(values)) else 0,
                     gain, reading, verdict))
    rows.sort(key=lambda r: (r[2] if r[2] != 'NA' else -1))
    return rows


def homolog_counts(proteomes, hits, genome_ids):
    """Per genome, how many of its proteins are hit from EVERY other genome in the set.

    This is the author's original Process D measure. It is not a diagnostic afterthought: the
    genome with the FEWEST homologs is the one that anchors the core, below.
    """
    counts = {}
    for genome_id in genome_ids:
        start, stop = proteomes.offsets[genome_id]
        others = [g for g in genome_ids if g != genome_id]
        if not others:
            counts[genome_id] = 0
            continue
        counts[genome_id] = int(hits.hit_from[start:stop, others].all(axis=1).sum())
    return counts


def select_reference_genome(proteomes, counts, genome_ids):
    """The genome with the fewest homologs anchors the core. This is the author's rule.

    It is not an arbitrary pick. The core genome is by definition contained in every genome's
    gene set, so the genome sharing the fewest genes gives the TIGHTEST superset of the core -
    the most efficient pivot for the intersection, and the most conservative one. Anchoring the
    core on a direct reciprocal hit from this genome to every other is also stricter than
    joining orthologs transitively through a graph, which can chain through errors.

    Ties are broken by name so that a run is reproducible.
    """
    if not genome_ids:
        raise SystemExit('No genomes to choose a reference from.')
    reference = min(sorted(genome_ids, key=lambda g: proteomes.genomes[g]),
                    key=lambda g: counts[g])
    return reference


def core_proteins(proteomes, partner, genome_ids, reference):
    """Reference proteins having a reciprocal best hit in every other genome of the set.

    Returns the reference protein indices, in genome order, which fixes the ortholog group
    order and so the index alignment of every table written later.
    """
    start, stop = proteomes.offsets[reference]
    others = [g for g in genome_ids if g != reference]
    indices = np.arange(start, stop, dtype=np.int64)
    if not others:
        return indices
    present = partner[start:stop, others] >= 0
    return indices[present.all(axis=1)]


def core_for_subset(proteomes, hits, partner, genome_ids, reference=None):
    """Apply the author's rule to any subset: pick its reference, then take the core.

    Re-selecting the reference within the subset is what lets leave-one-out and greedy removal
    work at all - a reference chosen on the full set cannot answer a question about a subset
    that excludes it.
    """
    if reference is None:
        counts = homolog_counts(proteomes, hits, genome_ids)
        reference = select_reference_genome(proteomes, counts, genome_ids)
    return reference, core_proteins(proteomes, partner, genome_ids, reference)


def build_presence_matrix(proteomes, partner, genome_ids, reference):
    """bool [genome, reference protein]: does this genome have an ortholog of that protein?

    Anchored on the reference, this is the pangenome matrix. Every column is one ortholog
    group, named by the reference gene that defines it.
    """
    start, stop = proteomes.offsets[reference]
    matrix = np.zeros((len(proteomes.genomes), stop - start), dtype=bool)
    for genome_id in genome_ids:
        if genome_id == reference:
            matrix[genome_id, :] = True
        else:
            matrix[genome_id, :] = partner[start:stop, genome_id] >= 0
    return matrix


def core_group_indices(matrix, genome_ids, fraction=1.0):
    """Columns present in at least `fraction` of the given genomes.

    fraction=1.0 is the strict core; a lower value gives a soft core, standard practice
    precisely because a strict intersection is brittle.
    """
    if not len(genome_ids):
        return np.array([], dtype=np.int64)
    submatrix = matrix[np.asarray(genome_ids, dtype=np.int64), :]
    needed = len(genome_ids) if fraction >= 1.0 else int(np.ceil(fraction * len(genome_ids)))
    return np.flatnonzero(submatrix.sum(axis=0) >= needed)


def soft_core_curve(matrix, genome_ids):
    """Core size as the presence requirement is relaxed from 100% downwards."""
    rows = []
    total = len(genome_ids)
    for present in range(total, 0, -1):
        fraction = present / total
        size = len(core_group_indices(matrix, genome_ids, fraction))
        rows.append((present, total, round(fraction * 100, 1), size))
    return rows


def leave_one_out_core(proteomes, hits, partner, genome_ids):
    """Core size with each genome excluded in turn.

    The reference is re-selected inside every subset by the author's rule, so excluding the
    reference genome itself is handled naturally rather than being a special case.

    If dropping one genome restores a large core while dropping any other does not, that
    genome is the outlier - the cheapest direct test available for the question this program
    exists to answer.
    """
    _, baseline_core = core_for_subset(proteomes, hits, partner, genome_ids)
    baseline = len(baseline_core)
    rows = []
    for excluded in genome_ids:
        remaining = [g for g in genome_ids if g != excluded]
        if len(remaining) < 2:
            continue
        _, core = core_for_subset(proteomes, hits, partner, remaining)
        rows.append((proteomes.genomes[excluded], len(core), len(core) - baseline))
    rows.sort(key=lambda row: -row[1])
    return baseline, rows


def greedy_removal_curve(proteomes, hits, partner, genome_ids):
    """Repeatedly drop the genome whose removal most increases the core.

    Core size is monotonic - adding a genome can only shrink the core - so this ordering is
    principled, and it sidesteps the impossibility of enumerating subsets (2^100 for a hundred
    genomes).

    Evaluating every candidate by rebuilding its core would be O(N^3) in protein comparisons.
    Instead each round makes ONE pass to count, per reference protein, how many of the
    remaining genomes carry it. A protein survives the removal of genome c when either every
    genome still has it, or exactly one lacks it and that one is c. Both are counted from the
    row sums, so each candidate then costs a single O(proteins) test.
    """
    remaining = list(genome_ids)
    reference, core = core_for_subset(proteomes, hits, partner, remaining)
    rows = [(len(remaining), len(core), '')]
    while len(remaining) > 2:
        counts = homolog_counts(proteomes, hits, remaining)
        reference = select_reference_genome(proteomes, counts, remaining)
        start, stop = proteomes.offsets[reference]
        others = [g for g in remaining if g != reference]
        present = partner[start:stop, others] >= 0
        row_sums = present.sum(axis=1)
        total = len(others)
        complete = int((row_sums == total).sum())          # kept whatever is removed
        one_short = row_sums == (total - 1)                 # kept only if the missing one goes
        best_genome, best_size = None, -1
        for position, candidate in enumerate(others):
            size = complete + int((one_short & ~present[:, position]).sum())
            if size > best_size:
                best_genome, best_size = candidate, size
        if best_genome is None:
            break
        remaining.remove(best_genome)
        _, core = core_for_subset(proteomes, hits, partner, remaining)
        rows.append((len(remaining), len(core), proteomes.genomes[best_genome]))
    return rows


# --------------------------------------------------------------------------------------
# output
# --------------------------------------------------------------------------------------

def prepare_output_directory(output_dir, force, keep_working=True):
    """Create the results directory, refusing to clobber existing data unless forced.

    With --force the results are replaced but the working directory is kept, so a finished
    search is not thrown away merely because the analysis is being rerun with new thresholds.
    """
    if os.path.exists(output_dir):
        if not os.path.isdir(output_dir):
            raise SystemExit('Output path exists and is not a directory: ' + output_dir)
        if os.listdir(output_dir):
            if not force:
                raise SystemExit('Output directory is not empty: %s\n'
                                 'Use --force to delete and recreate it, or choose another '
                                 '--output.' % output_dir)
            LOG.warning('--force given, replacing results in %s', output_dir)
            for entry in os.listdir(output_dir):
                if keep_working and entry == 'working':
                    continue
                target = os.path.join(output_dir, entry)
                shutil.rmtree(target) if os.path.isdir(target) else os.remove(target)
    os.makedirs(output_dir, exist_ok=True)
    return output_dir


def write_table(path, header, rows):
    """Write one tab-separated table."""
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, 'w') as handle:
        handle.write('\t'.join(header) + '\n')
        for row in rows:
            handle.write('\t'.join(str(value) for value in row) + '\n')
    return path


def write_matrix_table(path, labels, matrix, value_format):
    """Write a square genome-by-genome matrix as a labelled TSV."""
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, 'w') as handle:
        handle.write('genome\t' + '\t'.join(labels) + '\n')
        for row_index, name in enumerate(labels):
            values = []
            for column in range(len(labels)):
                value = matrix[row_index, column]
                values.append('NA' if value != value else value_format % value)
            handle.write(name + '\t' + '\t'.join(values) + '\n')
    return path


def write_presence_matrix(path, proteomes, matrix):
    """Write the pangenome matrix: one row per genome, one column per ortholog group."""
    with open(path, 'w') as handle:
        handle.write('genome\t' + '\t'.join('group%d' % i for i in range(matrix.shape[1])) + '\n')
        for genome_id, name in enumerate(proteomes.genomes):
            handle.write(name + '\t' + '\t'.join('1' if v else '0'
                                                 for v in matrix[genome_id]) + '\n')
    return path


def write_summary(path, proteomes, settings, reference, core_size, homologs,
                  leave_out_rows, greedy_rows, speciation_rows, bin_rows):
    """Write the human-readable verdict: what was found, and what to do next.

    The report is deliberately written as a STARTING POINT. Its job is to let a user say
    "the screen found this, which is why we then did X" - where X is ANI, dDDH, or wet-lab
    characterisation. It never claims to have settled a taxonomic question by itself.
    """
    genome_count = len(proteomes.genomes)
    lines = []
    add = lines.append
    add('RaGCAn %s - core genome screen' % __version__)
    add('=' * 72)
    add('')
    add('genomes    : %d' % genome_count)
    add('proteins   : %s' % human_count(len(proteomes.names)))
    add('reference  : %s (fewest homologs, %d - it anchors the core)'
        % (proteomes.genomes[reference], homologs[proteomes.genomes[reference]]))
    add('')
    add('thresholds : identity >= %g%%, gap opens <= %d, e-value < %g, coverage >= %g%%'
        % (settings.pident, settings.gapopen, settings.evalue, settings.min_coverage * 100))
    add('binning    : AAI >= %g%% (complete linkage)' % settings.bin_aai)
    add('')
    add('CORE GENOME: %d genes shared by all %d genomes' % (core_size, genome_count))
    add('')
    if core_size == 0:
        add('*** THE CORE COLLAPSED - no gene is shared by every genome. ***')
        add('')
        add('This is a result, not a failure. Before reading it as biology, rule out the')
        add('mundane causes: a fragmented or contaminated assembly, a proteome called with a')
        add('different gene caller, or an identity threshold too strict for cross-species work.')
        add('')
        add('If those are excluded, the set as given is not one coherent taxon. The proposed')
        add('bins below are the re-categorisation the data supports.')
        add('')
    add('-' * 72)
    add('PROPOSED BINS - each bin is a candidate sub-taxon')
    add('')
    add('  %3s %7s %10s  %8s %8s  %s'
        % ('bin', 'genomes', 'core genes', 'min AAI', 'mean AAI', 'members'))
    for index, size, core, _ref, lowest, average, members in bin_rows:
        add('  %3d %7d %10d  %8s %8s  %s' % (index, size, core, lowest, average, members[:80]))
    add('')
    if len(bin_rows) == 1:
        add('  One bin: the set is coherent at this threshold. No re-categorisation implied.')
    else:
        add('  %d bins: the set as given does NOT form one coherent group at AAI >= %g%%.'
            % (len(bin_rows), settings.bin_aai))
        singletons = [row for row in bin_rows if row[1] == 1]
        if singletons:
            add('  %d genome(s) bin alone, belonging with nothing else in the set:'
                % len(singletons))
            for row in singletons:
                add('      %s' % row[6])
    add('')
    add('-' * 72)
    add('SPECIATION PREDICTION')
    add('')
    add('Two independent lines of evidence per genome: average amino acid identity (AAI) to')
    add('its nearest neighbour, and how much the shared core grows when it is removed. A')
    add('claim is worth pursuing when the two AGREE; either on its own is weak evidence.')
    add('')
    add('  %-28s %-22s %7s %9s  %s'
        % ('genome', 'nearest', 'AAI%', 'core gain', 'verdict'))
    for name, nearest, aai_value, _shared, gain, _reading, verdict in speciation_rows:
        add('  %-28s %-22s %7s %+9d  %s'
            % (name[:28], nearest[:22], aai_value, gain, verdict))
    add('')
    add('  Conventional boundaries (guidelines from the literature, not laws):')
    add('    AAI >= %g%%   same species' % AAI_SPECIES)
    add('    %g - %g%%    same genus, different species' % (AAI_DIVERGENT, AAI_SPECIES))
    add('    %g - %g%%    same genus, divergent' % (AAI_GENUS, AAI_DIVERGENT))
    add('    AAI <  %g%%   likely a different genus' % AAI_GENUS)
    add('')
    add('-' * 72)
    add('HOMOLOGS PER GENOME (proteins hit from every other genome)')
    add('')
    for name, count in sorted(homologs.items(), key=lambda kv: kv[1]):
        add('  %-48s %6d' % (name[:48], count))
    add('')
    add('-' * 72)
    add('LEAVE-ONE-OUT (core size with each genome excluded)')
    add('')
    add('  %-48s %8s %8s' % ('excluded genome', 'core', 'change'))
    for name, size, delta in leave_out_rows:
        add('  %-48s %8d %+8d' % (name[:48], size, delta))
    add('')
    add('-' * 72)
    add('GREEDY REMOVAL (core size as the worst genome is dropped repeatedly)')
    add('')
    add('  %8s %8s   %s' % ('genomes', 'core', 'genome removed'))
    for count, size, removed in greedy_rows:
        add('  %8d %8d   %s' % (count, size, removed))
    add('')
    add('=' * 72)
    add('WHAT THIS IS, AND WHAT TO DO NEXT')
    add('')
    add('This is a SCREEN, not a taxonomic conclusion. It is designed to be the first step:')
    add('to let you say "the screen indicated this, which is why we then did the following".')
    add('')
    add('If a genome is flagged as a candidate novel taxon, or the bins disagree with the')
    add('current classification, the conventional next steps are:')
    add('')
    add('  dry lab   ANI against type strains (< 95-96% supports a distinct species)')
    add('            dDDH (< 70% supports a distinct species)')
    add('            16S rRNA identity (< 98.7% supports a candidate new species)')
    add('            a phylogeny: align the ortholog groups written by this run with MAFFT,')
    add('            concatenate, and infer a tree with IQ-TREE')
    add('')
    add('  wet lab   phenotypic and chemotaxonomic characterisation of the candidate:')
    add('            morphology, growth range, biochemistry, fatty acid and polar lipid')
    add('            profiles, respiratory quinones - as required to describe a novel taxon')
    add('')
    add('Core gene sequences are written unaligned, per genome and per ortholog group, so')
    add('the phylogenetic step can start directly from this output.')
    text = '\n'.join(lines) + '\n'
    with open(path, 'w') as handle:
        handle.write(text)
    return text


def group_members(proteomes, partner, ref_protein, genome_ids, reference):
    """The ortholog group defined by one reference protein: {genome id: protein index}."""
    members = {reference: int(ref_protein)}
    for genome_id in genome_ids:
        if genome_id == reference:
            continue
        mate = int(partner[ref_protein, genome_id])
        if mate >= 0:
            members[genome_id] = mate
    return members


def build_core_table(partner, core_indices, genome_ids, reference):
    """{genome id: [protein index per core gene]}, in core gene order.

    Computed once and shared by every writer. Doing the lookup inside each writer's nested
    loop instead would repeat it genomes x genes x genomes times.
    """
    columns = np.asarray(core_indices, dtype=np.int64)
    table = {}
    for genome_id in genome_ids:
        if genome_id == reference:
            table[genome_id] = columns.tolist()
        else:
            table[genome_id] = partner[columns, genome_id].astype(np.int64).tolist()
    return table


def write_core_gene_tables(out_dir, proteomes, core_table, genome_ids):
    """Write the index-aligned core gene name table for each genome.

    Row i of every table is the same ortholog group - the property the original program was
    built around, preserved here.
    """
    os.makedirs(out_dir, exist_ok=True)
    names = proteomes.names
    for genome_id in genome_ids:
        path = os.path.join(out_dir, proteomes.genomes[genome_id] + '.txt')
        with open(path, 'w') as handle:
            handle.write(''.join(names[i] + '\n' for i in core_table[genome_id]))
    return out_dir


def write_core_sequences(out_dir, proteomes, core_table, genome_ids):
    """Write each genome's core gene sequences, in ortholog group order."""
    os.makedirs(out_dir, exist_ok=True)
    names, sequences = proteomes.names, proteomes.sequences
    for genome_id in genome_ids:
        path = os.path.join(out_dir, proteomes.genomes[genome_id] + '.fasta')
        with open(path, 'w') as handle:
            handle.write(''.join('>%s\n%s\n' % (names[i], sequences[i])
                                 for i in core_table[genome_id]))
    return out_dir


def write_ortholog_group_fasta(out_dir, proteomes, core_table, genome_ids, gene_count):
    """Write one unaligned FASTA per core ortholog group, ready for MAFFT.

    Headers are the genome name, so an alignment of these files concatenates into a supermatrix
    directly and a tree built from it already carries sensible tip labels.
    """
    os.makedirs(out_dir, exist_ok=True)
    names, sequences, genomes = proteomes.names, proteomes.sequences, proteomes.genomes
    width = len(str(max(1, gene_count - 1)))
    for position in range(gene_count):
        path = os.path.join(out_dir, 'group%s.fasta' % str(position).zfill(width))
        with open(path, 'w') as handle:
            handle.write(''.join(
                '>%s %s\n%s\n' % (genomes[g], names[core_table[g][position]],
                                   sequences[core_table[g][position]])
                for g in genome_ids))
    return out_dir


# --------------------------------------------------------------------------------------
# binning - the proposed re-categorisation
# --------------------------------------------------------------------------------------

def bin_genomes(aai, genome_ids, threshold):
    """Partition genomes into bins, each a candidate sub-taxon.

    Complete-linkage agglomerative clustering on AAI: two bins merge only when EVERY genome in
    one is at or above `threshold` identity to EVERY genome in the other. Complete linkage,
    rather than single linkage, is what makes a bin defensible as a taxon - single linkage
    would let a chain of intermediates drag two unrelated genomes into one group.

    Merging stops when no admissible merge remains, so the number of bins is discovered from
    the data rather than fixed in advance. One bin means the set is coherent; more than one is
    a proposed re-categorisation; a singleton bin is a genome that belongs with nothing else.

    Also returns the merge trace: one row per merge, in the order it happened, recording the
    linkage value (the worst-case AAI that justified it) and the two sides merged. A single
    threshold collapses this whole history into one yes/no call per bin; the trace is what
    lets that call be audited - e.g. seeing that a large bin exists only because of one pair
    sitting barely above threshold, rather than a clean gap in the AAI distribution.
    """
    bins = [[g] for g in genome_ids]
    # linkage[i, j] is the WORST AAI between any member of bin i and any member of bin j.
    # Complete linkage lets it be updated on merge as the minimum of the two rows, so the
    # whole clustering costs O(N^2) instead of rescanning every cross pair every round.
    order = np.asarray(genome_ids, dtype=np.int64)
    linkage = aai[np.ix_(order, order)].astype(np.float64, copy=True)
    linkage[np.isnan(linkage)] = -np.inf
    np.fill_diagonal(linkage, -np.inf)
    alive = np.ones(len(bins), dtype=bool)
    trace = []
    while alive.sum() > 1:
        masked = np.where(alive[:, None] & alive[None, :], linkage, -np.inf)
        flat = int(np.argmax(masked))
        i, j = divmod(flat, masked.shape[1])
        if masked[i, j] < threshold:
            break
        if i > j:
            i, j = j, i
        trace.append((float(masked[i, j]), list(bins[i]), list(bins[j])))
        bins[i] = sorted(bins[i] + bins[j])
        bins[j] = []
        alive[j] = False
        linkage[i, :] = np.minimum(linkage[i, :], linkage[j, :])
        linkage[:, i] = linkage[i, :]
        linkage[i, i] = -np.inf
    bins = [members for members in bins if members]
    bins.sort(key=lambda members: (-len(members), members[0]))
    return bins, trace


def describe_bins(proteomes, hits, partner, bins, aai):
    """Summarise each proposed bin: size, members, internal AAI, and its own core genome."""
    rows = []
    for index, members in enumerate(bins, start=1):
        reference, core = core_for_subset(proteomes, hits, partner, members)
        values = [aai[a, b] for a in members for b in members
                  if a < b and aai[a, b] == aai[a, b]]
        lowest = round(min(values), 2) if values else 'NA'
        average = round(sum(values) / len(values), 2) if values else 'NA'
        rows.append((index, len(members), len(core), proteomes.genomes[reference],
                     lowest, average,
                     ' '.join(proteomes.genomes[g] for g in members)))
    return rows


def merge_trace_rows(proteomes, trace):
    """Flatten bin_genomes' merge trace into table rows, one per merge step in order.

    side_a/side_b list every genome on each side of the merge (semicolon-separated) so the
    exact bridging pair behind any surprising merge can be found without recomputing anything.
    """
    rows = []
    for step, (linkage_aai, side_a, side_b) in enumerate(trace, start=1):
        rows.append((step, round(linkage_aai, 2), len(side_a), len(side_b),
                     len(side_a) + len(side_b),
                     ';'.join(proteomes.genomes[g] for g in side_a),
                     ';'.join(proteomes.genomes[g] for g in side_b)))
    return rows


def run_screen(settings, tools):
    """Drive every stage in order and return the summary text."""
    out = settings.output_dir
    prepare_output_directory(out, settings.force, keep_working=not settings.rerun_search)
    work_dir = os.path.join(out, 'working')
    os.makedirs(work_dir, exist_ok=True)

    LOG.info('[1] loading proteomes')
    proteomes = load_proteomes(settings.input_dir)
    pooled = write_pooled_fasta(proteomes, os.path.join(work_dir, 'pooled.faa'))

    estimate = len(proteomes.names) * len(proteomes.genomes) * 13 / 1e9
    LOG.info('hit tables will need about %.2f GB of RAM (%s proteins x %d genomes)',
             estimate, human_count(len(proteomes.names)), len(proteomes.genomes))
    if estimate > 4.0:
        LOG.warning('that is a large allocation; if the run is killed, reduce the number of '
                    'genomes or run on a machine with more memory')

    LOG.info('[2] all-against-all search')
    hits_path = run_all_vs_all_search(tools, pooled, work_dir, settings)

    LOG.info('[3] parsing hits')
    hits = parse_best_hits(hits_path, proteomes, settings)

    LOG.info('[4] reciprocal best hits')
    partner, left, right, identity = reciprocal_best_hits(proteomes, hits, settings.pident)

    LOG.info('[5] average amino acid identity')
    aai, shared = pairwise_identity_matrices(proteomes, left, right, identity, settings.threads)

    LOG.info('[6] core genome')
    genome_ids = list(range(len(proteomes.genomes)))
    counts = homolog_counts(proteomes, hits, genome_ids)
    reference = select_reference_genome(proteomes, counts, genome_ids)
    LOG.info('reference genome (fewest homologs): %s (%d homologs)',
             proteomes.genomes[reference], counts[reference])
    matrix = build_presence_matrix(proteomes, partner, genome_ids, reference)
    core_indices = core_proteins(proteomes, partner, genome_ids, reference)
    LOG.info('core genome: %d genes shared by all %d genomes',
             len(core_indices), len(genome_ids))

    LOG.info('[7] diagnostics')
    homologs = {proteomes.genomes[g]: counts[g] for g in genome_ids}
    baseline, leave_out_rows = leave_one_out_core(proteomes, hits, partner, genome_ids)
    greedy_rows = greedy_removal_curve(proteomes, hits, partner, genome_ids)
    speciation_rows = predict_speciation(proteomes, aai, shared, matrix, genome_ids,
                                         leave_out_rows)

    LOG.info('[8] binning into candidate sub-taxa')
    bins, merge_trace = bin_genomes(aai, genome_ids, settings.bin_aai)
    bin_rows = describe_bins(proteomes, hits, partner, bins, aai)
    LOG.info('%d bin(s) proposed at AAI >= %g%%', len(bins), settings.bin_aai)

    LOG.info('[9] writing results')
    reports = os.path.join(out, 'reports')
    write_presence_matrix(os.path.join(out, 'pangenome_matrix.tsv'), proteomes, matrix)
    write_table(os.path.join(reports, 'homolog_counts.tsv'), ['genome', 'homologs'],
                sorted(homologs.items(), key=lambda kv: kv[1]))
    write_table(os.path.join(reports, 'leave_one_out.tsv'),
                ['excluded_genome', 'core_size', 'change_vs_baseline'], leave_out_rows)
    write_table(os.path.join(reports, 'greedy_removal.tsv'),
                ['genomes_retained', 'core_size', 'genome_removed'], greedy_rows)
    write_table(os.path.join(reports, 'soft_core_curve.tsv'),
                ['genomes_present', 'genomes_total', 'percent', 'group_count'],
                soft_core_curve(matrix, genome_ids))
    write_table(os.path.join(reports, 'speciation_prediction.tsv'),
                ['genome', 'nearest_genome', 'aai_percent', 'shared_orthologs',
                 'core_gain_if_removed', 'aai_reading', 'verdict'], speciation_rows)
    write_table(os.path.join(reports, 'proposed_bins.tsv'),
                ['bin', 'genomes', 'core_genes', 'bin_reference', 'min_aai', 'mean_aai',
                 'members'], bin_rows)
    write_table(os.path.join(reports, 'merge_trace.tsv'),
                ['step', 'linkage_aai', 'size_a', 'size_b', 'merged_size', 'side_a', 'side_b'],
                merge_trace_rows(proteomes, merge_trace))
    write_matrix_table(os.path.join(reports, 'aai_matrix.tsv'), proteomes.genomes, aai, '%.2f')
    write_matrix_table(os.path.join(reports, 'shared_orthologs.tsv'), proteomes.genomes,
                       shared, '%d')

    if len(core_indices):
        core_table = build_core_table(partner, core_indices, genome_ids, reference)
        write_core_gene_tables(os.path.join(out, 'core_genes', 'gene_names'),
                               proteomes, core_table, genome_ids)
        write_core_sequences(os.path.join(out, 'core_genes', 'per_genome'),
                             proteomes, core_table, genome_ids)
        write_ortholog_group_fasta(os.path.join(out, 'core_genes', 'ortholog_groups'),
                                   proteomes, core_table, genome_ids, len(core_indices))
    else:
        LOG.warning('core is empty - no core gene sequences written; '
                    'see reports/proposed_bins.tsv for the suggested re-categorisation')

    summary = write_summary(os.path.join(reports, 'summary.txt'), proteomes, settings,
                            reference, len(core_indices), homologs, leave_out_rows,
                            greedy_rows, speciation_rows, bin_rows)
    if settings.clean:
        shutil.rmtree(work_dir, ignore_errors=True)
    return summary


# --------------------------------------------------------------------------------------
# command line
# --------------------------------------------------------------------------------------

class Settings(NamedTuple):
    """Everything the screen needs, fixed once at start-up and never mutated."""
    input_dir: str
    output_dir: str
    threads: int
    pident: float
    gapopen: int
    evalue: float
    min_coverage: float
    bin_aai: float
    block_size: float
    index_chunks: int
    rerun_search: bool
    sensitivity: str
    max_target_seqs: int
    force: bool
    clean: bool


def build_parser():
    """Define the command-line interface."""
    parser = argparse.ArgumentParser(
        prog='RaGCAn',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='Rapid core-genome screen for taxonomic coherence.',
        epilog='example:\n'
               '  RaGCAn.py -i proteomes/ -o results/ -t 16\n\n'
               'Needs only DIAMOND on PATH:\n'
               '  micromamba env create -f environment.yml && micromamba activate ragcan\n')
    parser.add_argument('-i', '--input', required=True, metavar='DIR',
                        help='directory of protein FASTA files, one per genome (never modified)')
    parser.add_argument('-o', '--output', metavar='DIR',
                        help='directory for results (default: <input>_Result)')
    parser.add_argument('-t', '--threads', type=int, default=max(1, (os.cpu_count() or 2) - 1),
                        metavar='N', help='threads for DIAMOND (default: all but one, here %d)'
                                          % max(1, (os.cpu_count() or 2) - 1))
    group = parser.add_argument_group('orthology thresholds')
    group.add_argument('--pident', type=float, default=85.0, metavar='PCT',
                       help='minimum percentage identity (default: 85, the original value; '
                            'consider 50-70 for cross-species work)')
    group.add_argument('--gapopen', type=int, default=3, metavar='N',
                       help='maximum gap opens in a hit (default: 3)')
    group.add_argument('--evalue', type=float, default=1e-4, metavar='E',
                       help='maximum e-value (default: 1e-4, what the original\'s "10e-5" '
                            'actually evaluates to)')
    group.add_argument('--min-coverage', type=float, default=0.0, metavar='FRAC',
                       help='minimum aligned fraction of BOTH sequences, 0-1 '
                            '(default: 0, off; 0.5 is a common choice)')
    group.add_argument('--bin-aai', type=float, default=AAI_GENUS, metavar='PCT',
                       help='AAI at or above which two genomes may share a bin '
                            '(default: %g, the conventional genus boundary; use %g for '
                            'species-level bins)' % (AAI_GENUS, AAI_SPECIES))
    search = parser.add_argument_group('search')
    search.add_argument('--sensitivity', default='very-sensitive',
                        choices=['default', 'mid-sensitive', 'sensitive',
                                 'more-sensitive', 'very-sensitive', 'ultra-sensitive'],
                        help='DIAMOND sensitivity (default: very-sensitive; the DIAMOND '
                             'default is tuned for short reads and misses diverged orthologs)')
    search.add_argument('--block-size', type=float, default=0.4, metavar='GB',
                        help='DIAMOND block size (default: 0.4). Peak memory is roughly six '
                             'times this in GB, so 0.4 needs about 2.4 GB. Raise it for speed '
                             'if you have RAM to spare')
    search.add_argument('--index-chunks', type=int, default=4, metavar='N',
                        help='DIAMOND index chunks (default: 4). More chunks means less '
                             'memory and a slower search')
    search.add_argument('--max-target-seqs', type=int, default=500, metavar='N',
                        help='hits kept per query (default: 500; DIAMOND\'s own default of 25 '
                             'can truncate all-against-all runs)')
    behaviour = parser.add_argument_group('behaviour')
    behaviour.add_argument('--rerun-search', action='store_true',
                           help='discard a previously finished search and run it again. '
                                'Without this a finished search in the output directory is '
                                'reused, so thresholds can be changed without paying for it '
                                'twice')
    behaviour.add_argument('--force', action='store_true',
                           help='delete and recreate the output directory if it is not empty')
    behaviour.add_argument('--clean', action='store_true',
                           help='delete the working files when finished')
    behaviour.add_argument('--log-file', metavar='FILE', help='also write the log to FILE')
    verbosity = behaviour.add_mutually_exclusive_group()
    verbosity.add_argument('-v', '--verbose', action='store_true', help='show every command run')
    verbosity.add_argument('-q', '--quiet', action='store_true', help='warnings and errors only')
    parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)
    return parser


def configure_logging(verbose, quiet):
    """Send progress to stderr. Handlers are reset so a second call cannot double every line."""
    LOG.setLevel(logging.DEBUG)
    LOG.handlers.clear()
    stream = logging.StreamHandler()
    stream.setLevel(logging.DEBUG if verbose else logging.WARNING if quiet else logging.INFO)
    stream.setFormatter(logging.Formatter('%(levelname)-8s %(message)s'))
    LOG.addHandler(stream)


def add_log_file(log_file):
    """Tee the log to a file, always at full detail whatever -q/-v say."""
    os.makedirs(os.path.dirname(os.path.abspath(log_file)), exist_ok=True)
    handler = logging.FileHandler(log_file, mode='w')
    handler.setLevel(logging.DEBUG)
    handler.setFormatter(logging.Formatter('%(asctime)s %(levelname)-8s %(message)s'))
    LOG.addHandler(handler)


def settings_from_args(args):
    """Turn the parsed command line into validated Settings.

    Rejects any combination that would destroy data, resolving both directories with realpath
    so a symlink cannot slip past the string containment tests.
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
    # --force deletes the whole output tree, so an output directory that CONTAINS the input
    # would destroy the proteomes it is about to read
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
    if not 0.0 <= args.min_coverage <= 1.0:
        raise SystemExit('--min-coverage must be between 0 and 1.')
    if args.block_size <= 0:
        raise SystemExit('--block-size must be greater than 0.')
    if args.index_chunks < 1:
        raise SystemExit('--index-chunks must be 1 or more.')
    if not 0 < args.bin_aai <= 100:
        raise SystemExit('--bin-aai must be between 0 and 100.')
    if args.max_target_seqs < 1:
        raise SystemExit('--max-target-seqs must be 1 or more.')
    return Settings(input_dir=input_dir, output_dir=output_dir, threads=args.threads,
                    pident=args.pident, gapopen=args.gapopen, evalue=args.evalue,
                    min_coverage=args.min_coverage, bin_aai=args.bin_aai,
                    block_size=args.block_size, index_chunks=args.index_chunks,
                    rerun_search=args.rerun_search,
                    sensitivity='' if args.sensitivity == 'default' else args.sensitivity,
                    max_target_seqs=args.max_target_seqs, force=args.force, clean=args.clean)


def validated_log_path(log_file, output_dir):
    """Resolve --log-file and keep it out of the output tree."""
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
        summary = run_screen(settings, tools)
    except SystemExit as exc:
        if isinstance(exc.code, str):
            LOG.error('%s', exc.code)
            return 1
        raise
    except KeyboardInterrupt:
        LOG.error('interrupted')
        return 130
    print()
    print(summary)
    return 0


if __name__ == '__main__':
    sys.exit(main())
