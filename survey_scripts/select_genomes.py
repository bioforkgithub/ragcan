#!/usr/bin/env python3
"""
Phase 1 genome selection for the LPSN genus-coherence survey.

Selection rule, applied uniformly to every genus (no exceptions, no per-genus tuning):

  1. version_status == "latest"
  2. refseq_category in {"reference genome", "representative genome"}
  3. organism_name is a proper binomial: token[0] is the genus (capitalised),
     token[1] is a lowercase alphabetic specific epithet.
     Rejects "sp.", "cf.", strain-only designations, epithets with digits or
     punctuation, and anything not resolving to a named species.
  4. ONE genome per named species (genus + epithet). Ties broken deterministically:
     "reference genome" beats "representative genome"; then assembly level
     (Complete genome > Chromosome > Scaffold > Contig); then highest accession.

*** NO QUALITY FILTER IS APPLIED. THIS IS DELIBERATE. ***

PLAN.md section 3 says to drop proteomes more than 1 SD below the genus median
protein count. That instruction is SUPERSEDED - see PROJECT_STATE.md sections 4.1,
12.2 and 12.3. Measured on 2026-08-26: that filter removes all seven
*Halopseudomonas* from the 389-genome *Pseudomonas* set, i.e. it deletes exactly the
reduced-genome lineage the screen exists to find. Filtering on assembly level is
also wrong ("complete" does not mean "gene-rich": *C. taklimakanense* is a complete
assembly with 2,553 proteins). The agreed design is one run per genus with drafts
INCLUDED, and assembly level plus protein count recorded as COVARIATES, analysed
afterwards. If you are reading this because you are about to add a filter here: do
not. Read PROJECT_STATE.md section 12.3 first.

Writes reports/genome_selection.tsv - the manifest every later stage reads.
"""

import argparse
import csv
import os
import re
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SUMMARY = os.path.join(ROOT, 'reference', 'assembly_summary_bacteria.txt')
OUT = os.path.join(ROOT, 'reports', 'genome_selection.tsv')

PILOT_GENERA = [
    # expected coherent
    'Chryseobacterium', 'Stutzerimonas', 'Halopseudomonas',
    'Bordetella', 'Brucella', 'Yersinia',
    # expected problematic / split by GTDB
    'Pseudomonas', 'Clostridium', 'Bacillus',
    'Lactobacillus', 'Streptococcus', 'Mycobacterium',
    # recently split, documented answer
    'Aeromonas', 'Ralstonia', 'Paraburkholderia', 'Burkholderia',
]

# column indices in assembly_summary_bacteria.txt (0-based)
C_ACC, C_CAT, C_NAME, C_VER, C_LEVEL, C_FTP = 0, 4, 7, 10, 11, 19
C_TOTAL_GENE, C_PROT_GENE = 34, 35

CATEGORIES = {'reference genome', 'representative genome'}
CAT_RANK = {'reference genome': 0, 'representative genome': 1}
LEVEL_RANK = {'Complete Genome': 0, 'Chromosome': 1, 'Scaffold': 2, 'Contig': 3}

EPITHET = re.compile(r'^[a-z]+$')


def safe_name(organism_name):
    """Reproduce the file-naming convention already used by this project.

    Existing proteome files are named <organism_name with unsafe characters
    replaced by _>__<accession>.faa, e.g.
      'Pseudomonas amygdali pv. tabaci str. ATCC 11528'
        -> Pseudomonas_amygdali_pv._tabaci_str._ATCC_11528__GCF_000145945.2.faa
      'Pseudomonas cremoricolorata DSM 17059 = NBRC 16634'
        -> Pseudomonas_cremoricolorata_DSM_17059___NBRC_16634__GCF_000425745.1.faa
    Verified against Stutzerimonas_Validation/proteomes_pseudomonas_full/.
    """
    return re.sub(r'[^A-Za-z0-9._-]', '_', organism_name)


def parse(genera):
    wanted = set(genera)
    picked = {}          # (genus, epithet) -> row dict
    rejected_no_ftp = []
    for_genus_rows = defaultdict(int)

    with open(SUMMARY, newline='') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) <= C_FTP:
                continue
            name = f[C_NAME]
            sp = name.split(' ', 2)
            if len(sp) < 2:
                continue
            genus, epithet = sp[0], sp[1]
            if genus not in wanted:
                continue
            for_genus_rows[genus] += 1
            if f[C_VER] != 'latest':
                continue
            if f[C_CAT] not in CATEGORIES:
                continue
            if not EPITHET.match(epithet):
                continue        # rejects sp., cf., 'Sp1', epithets with digits/dots
            key = (genus, epithet)
            row = {
                'genus': genus,
                'species': '%s %s' % (genus, epithet),
                'accession': f[C_ACC],
                'organism_name': name,
                'refseq_category': f[C_CAT],
                'assembly_level': f[C_LEVEL],
                'protein_coding_gene_count': f[C_PROT_GENE] if len(f) > C_PROT_GENE else '',
                'total_gene_count': f[C_TOTAL_GENE] if len(f) > C_TOTAL_GENE else '',
                'ftp_path': f[C_FTP],
            }
            if row['ftp_path'] in ('', 'na'):
                rejected_no_ftp.append(row)
                continue
            prev = picked.get(key)
            if prev is None or _rank(row) < _rank(prev):
                picked[key] = row
    return picked, rejected_no_ftp


def _rank(row):
    return (CAT_RANK.get(row['refseq_category'], 9),
            LEVEL_RANK.get(row['assembly_level'], 9),
            # highest accession wins -> negate by using reverse string compare
            tuple(-ord(c) for c in row['accession']))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out', default=OUT)
    ap.add_argument('--genera', nargs='*', default=PILOT_GENERA)
    args = ap.parse_args()

    picked, no_ftp = parse(args.genera)

    rows = sorted(picked.values(), key=lambda r: (r['genus'], r['species']))
    for r in rows:
        r['faa_name'] = '%s__%s.faa' % (safe_name(r['organism_name']), r['accession'])

    tmp = args.out + '.tmp'
    cols = ['genus', 'species', 'accession', 'organism_name', 'refseq_category',
            'assembly_level', 'protein_coding_gene_count', 'total_gene_count',
            'ftp_path', 'faa_name']
    with open(tmp, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter='\t', lineterminator='\n')
        w.writeheader()
        for r in rows:
            w.writerow({c: r[c] for c in cols})
    os.replace(tmp, args.out)

    counts = defaultdict(int)
    for r in rows:
        counts[r['genus']] += 1
    print('genus\tselected')
    for g in args.genera:
        print('%s\t%d' % (g, counts[g]))
    print('TOTAL\t%d' % len(rows))
    if no_ftp:
        print('\nrows dropped for missing ftp_path: %d' % len(no_ftp), file=sys.stderr)
        for r in no_ftp:
            print('  %s\t%s' % (r['species'], r['accession']), file=sys.stderr)
    print('\nwrote %s' % args.out)


if __name__ == '__main__':
    main()
