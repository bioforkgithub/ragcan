#!/usr/bin/env python3
"""
Turn one finished RaGCAn run into (a) a row of reports/survey_results.tsv, (b) the
genus's threshold curve, and (c) the boundary-crossing assertion that PROJECT_STATE.md
section 6.2a makes a CONDITION OF USE of the fast AAI path.

Everything here is read back off disk from results/<Genus>/reports/. Nothing is
recomputed from memory and nothing is estimated. If a file is missing, this raises
rather than filling in a plausible number - see the project rule, "verify, don't recall".

THE THRESHOLD CURVE
  Re-binning at other --bin-aai values does not need a re-run: it calls pygemini's own
  bin_genomes() (pygemini.py:944), the same function the pipeline calls at line 1071,
  on the AAI matrix the run wrote. Thresholds 65..85 inclusive, step 1.
  CAVEAT recorded in the output: reports/aai_matrix.tsv stores AAI at 2 decimal places,
  so the curve is computed at 2 dp. The bin count at 65 is cross-checked against
  pygemini's own proposed_bins.tsv (full precision) and any disagreement is recorded in
  the results row rather than hidden.

THE BOUNDARY ASSERTION
  PROJECT_STATE.md section 6.2a: the fast AAI path may be used in production only with a
  boundary-crossing check kept live, because the tightest observed approach to a
  threshold anywhere in the project was 0.0011 AAI units against a max original-vs-fast
  disagreement of 1.36e-05 (81x margin), and across many genera a tighter pair is
  possible.
  Running both implementations per genus would cost the speedup it buys, so the check
  used here is a sound SUPERSET test on the stored matrix. A pair can only flip between
  implementations if its true AAI lies within ~1.4e-05 of 65.000000 (or 95.000000). Any
  such value necessarily lies in [64.995, 65.005) and therefore stores as exactly
  "65.00" at 2 dp. So: flag every stored 65.00 and every stored 95.00. Zero flagged
  means no pair can have flipped. Anything flagged is written to
  reports/boundary_flags.tsv for a full-precision follow-up; it is never suppressed.
"""

import argparse
import csv
import math
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
PROJ = os.path.dirname(ROOT)
PYG = os.path.join(PROJ, 'pygemini_work', 'PY-GEMINI')
sys.path.insert(0, PYG)

import pygemini as pg   # noqa: E402  imported read-only, never modified

REPORTS = os.path.join(ROOT, 'reports')

THRESHOLDS = list(range(65, 86))

RESULT_COLS = [
    'genus', 'expectation', 'expectation_basis', 'verdict', 'agrees_with_expectation',
    'species_selected', 'genomes_used', 'total_proteins', 'core_genes',
    'bins_at_65', 'largest_bin', 'min_aai', 'mean_aai', 'pairs_ge_95',
    'loo_max_core_gain', 'loo_max_core_gain_genome',
    'first_split_threshold', 'bins_at_70', 'bins_at_75', 'bins_at_80', 'bins_at_85',
    'asm_complete', 'asm_chromosome', 'asm_scaffold', 'asm_contig', 'asm_other',
    'proteins_min', 'proteins_median', 'proteins_max',
    'proteomes_reused', 'proteomes_downloaded', 'species_no_annotation',
    'boundary_flags_65', 'boundary_flags_95',
    'bin65_matches_pygemini', 'wall_seconds', 'finished_utc',
]

# The pilot's registered expectations (PLAN.md Phase 1). Recorded so that a
# disagreement is visible the moment it happens - that disagreement IS the go/no-go
# signal, not a bug to be tidied away.
EXPECTED = {
    'Chryseobacterium': 'coherent', 'Stutzerimonas': 'coherent',
    'Halopseudomonas': 'coherent', 'Bordetella': 'coherent',
    'Brucella': 'coherent', 'Yersinia': 'coherent',
    'Pseudomonas': 'split', 'Clostridium': 'split', 'Bacillus': 'split',
    'Lactobacillus': 'split', 'Streptococcus': 'split', 'Mycobacterium': 'split',
    'Aeromonas': 'split', 'Ralstonia': 'split',
    'Paraburkholderia': 'split', 'Burkholderia': 'split',
}

# Where each expectation comes from, recorded because the three groups are NOT equally
# well specified and a MISMATCH must not be read as if they were.
#   PLAN.md Phase 1 states two firm groups ("should be coherent", "known to be
#   problematic / split by GTDB") and a third, "recently split, so the answer is
#   documented", which names four genera WITHOUT saying which way the bin count should
#   go. Mapping that third group to "expect >1 bin" is an interpretation made here, not
#   something PLAN.md asserts. It is labelled so, so that a mismatch in that group is
#   read as "the expectation was under-specified", not as "the screen failed".
EXPECT_BASIS = {
    'coherent': 'PLAN.md Phase 1: "should be coherent" - firm',
    'split': 'PLAN.md Phase 1: "known to be problematic / split by GTDB" - firm',
}
INTERPRETED = {'Aeromonas', 'Ralstonia', 'Paraburkholderia', 'Burkholderia'}
INTERPRETED_BASIS = ('PLAN.md Phase 1: "recently split, so the answer is documented" - '
                     'PLAN.md does not state a bin count; ">1 bin" is an interpretation '
                     'made by this pipeline, NOT an assertion of PLAN.md')


def read_matrix(path):
    names, rows = [], []
    with open(path) as fh:
        header = fh.readline().rstrip('\n').split('\t')[1:]
        for line in fh:
            f = line.rstrip('\n').split('\t')
            names.append(f[0])
            rows.append([float(x) if x not in ('', 'NA') else float('nan')
                         for x in f[1:]])
    aai = np.array(rows, dtype=np.float64)
    assert list(header) == names, 'aai_matrix.tsv is not square/symmetric in labels'
    return names, aai


def parse_summary(path):
    txt = open(path).read()
    out = {}
    m = re.search(r'^genomes\s*:\s*([\d,]+)', txt, re.M)
    out['genomes'] = int(m.group(1).replace(',', '')) if m else None
    m = re.search(r'^proteins\s*:\s*([\d,]+)', txt, re.M)
    out['proteins'] = int(m.group(1).replace(',', '')) if m else None
    m = re.search(r'CORE GENOME:\s*([\d,]+)\s*genes', txt)
    out['core'] = int(m.group(1).replace(',', '')) if m else None
    return out


def parse_bins(path):
    with open(path, newline='') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def parse_loo(path):
    best = (None, None)
    with open(path, newline='') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            try:
                d = int(r['change_vs_baseline'])
            except (KeyError, ValueError):
                continue
            if best[0] is None or d > best[0]:
                best = (d, r['excluded_genome'])
    return best


def covariates(genus):
    """Assembly level and protein-count covariates, from the fetch manifest.

    These are recorded, NOT filtered on. PROJECT_STATE.md sections 4.1 / 12.2 / 12.3:
    a protein-count filter deletes reduced-genome lineages (all 7 Halopseudomonas), and
    an assembly-level filter does not control gene count. Neither is safe, so the
    design is no pre-filter with quality as a covariate analysed afterwards.
    """
    path = os.path.join(REPORTS, 'fetch_%s.tsv' % genus)
    pdir = os.path.join(ROOT, 'proteomes', genus)
    lv = {'Complete Genome': 0, 'Chromosome': 0, 'Scaffold': 0, 'Contig': 0, 'other': 0}
    prot, reused, downloaded, noann = [], 0, 0, 0
    with open(path, newline='') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            if r['status'] == 'no_protein_annotation':
                noann += 1
                continue
            # Reused vs downloaded is decided by OBSERVING the file, not by the status
            # word: a symlink points at a proteome that already existed elsewhere in the
            # project, a regular file was fetched from NCBI by this pipeline. The status
            # word 'already_present' is ambiguous across re-runs and would over-report
            # reuse.
            if os.path.islink(os.path.join(pdir, r['faa_name'])):
                reused += 1
            else:
                downloaded += 1
            lv[r['assembly_level'] if r['assembly_level'] in lv else 'other'] += 1
            try:
                prot.append(int(r['protein_coding_gene_count']))
            except (TypeError, ValueError):
                pass
    prot.sort()
    med = prot[len(prot) // 2] if prot else ''
    return {
        'asm_complete': lv['Complete Genome'], 'asm_chromosome': lv['Chromosome'],
        'asm_scaffold': lv['Scaffold'], 'asm_contig': lv['Contig'],
        'asm_other': lv['other'],
        'proteins_min': prot[0] if prot else '', 'proteins_median': med,
        'proteins_max': prot[-1] if prot else '',
        'proteomes_reused': reused, 'proteomes_downloaded': downloaded,
        'species_no_annotation': noann,
    }


def boundary_check(genus, names, aai):
    """See module docstring. Returns (n65, n95) and writes any flags to disk."""
    n = len(names)
    iu = np.triu_indices(n, 1)
    vals = aai[iu]
    flags = []
    for boundary in (65.00, 95.00):
        hit = np.where(np.abs(vals - boundary) < 1e-9)[0]
        for k in hit:
            flags.append((genus, names[iu[0][k]], names[iu[1][k]],
                          '%.2f' % vals[k], boundary))
    # Replace this genus's rows rather than appending, so re-summarising a genus
    # cannot duplicate its flags.
    path = os.path.join(REPORTS, 'boundary_flags.tsv')
    keep = []
    if os.path.exists(path):
        with open(path, newline='') as fh:
            keep = [r for r in csv.DictReader(fh, delimiter='\t') if r['genus'] != genus]
    if flags or keep:
        tmp = path + '.tmp'
        with open(tmp, 'w', newline='') as fh:
            w = csv.writer(fh, delimiter='\t', lineterminator='\n')
            w.writerow(['genus', 'genome_a', 'genome_b', 'stored_aai', 'boundary'])
            for r in keep:
                w.writerow([r['genus'], r['genome_a'], r['genome_b'],
                            r['stored_aai'], r['boundary']])
            w.writerows(flags)
        os.replace(tmp, path)
    return (sum(1 for f in flags if f[4] == 65.00),
            sum(1 for f in flags if f[4] == 95.00))


def threshold_curve(genus, names, aai):
    """bin_genomes() at every threshold 65..85. Uses pygemini's own function."""
    ids = list(range(len(names)))
    rows = []
    for t in THRESHOLDS:
        bins, _trace = pg.bin_genomes(aai, ids, float(t))
        rows.append({'genus': genus, 'bin_aai': t, 'bins': len(bins),
                     'largest_bin': max(len(b) for b in bins),
                     'singletons': sum(1 for b in bins if len(b) == 1)})
    path = os.path.join(REPORTS, 'threshold_curves.tsv')
    existing = []
    if os.path.exists(path):
        with open(path, newline='') as fh:
            existing = [r for r in csv.DictReader(fh, delimiter='\t')
                        if r['genus'] != genus]
    tmp = path + '.tmp'
    with open(tmp, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=['genus', 'bin_aai', 'bins',
                                           'largest_bin', 'singletons'],
                           delimiter='\t', lineterminator='\n')
        w.writeheader()
        for r in existing + rows:
            w.writerow(r)
    os.replace(tmp, path)
    return rows


def summarise(genus, wall_seconds=''):
    rdir = os.path.join(ROOT, 'results', genus, 'reports')
    s = parse_summary(os.path.join(rdir, 'summary.txt'))
    bins = parse_bins(os.path.join(rdir, 'proposed_bins.tsv'))
    loo_gain, loo_genome = parse_loo(os.path.join(rdir, 'leave_one_out.tsv'))
    names, aai = read_matrix(os.path.join(rdir, 'aai_matrix.tsv'))

    n = len(names)
    iu = np.triu_indices(n, 1)
    vals = aai[iu]
    finite = vals[~np.isnan(vals)]
    min_aai = float(np.min(finite)) if finite.size else float('nan')
    mean_aai = float(np.mean(finite)) if finite.size else float('nan')
    pairs95 = int(np.sum(finite >= 95.0))

    curve = threshold_curve(genus, names, aai)
    by_t = {r['bin_aai']: r['bins'] for r in curve}
    first_split = next((t for t in THRESHOLDS if by_t[t] > 1), '')

    n65, n95 = boundary_check(genus, names, aai)

    bins_at_65 = len(bins)
    verdict = 'coherent' if bins_at_65 == 1 else 'split'
    exp = EXPECTED.get(genus, '')
    agrees = 'MATCH' if (exp and verdict == exp) else ('MISMATCH' if exp else 'NA')

    cov = covariates(genus)
    sel = sum(1 for _ in open(os.path.join(REPORTS, 'fetch_%s.tsv' % genus))) - 1

    row = {
        'genus': genus, 'expectation': exp, 'verdict': verdict,
        'expectation_basis': (INTERPRETED_BASIS if genus in INTERPRETED
                              else EXPECT_BASIS.get(exp, '')),
        'agrees_with_expectation': agrees,
        'species_selected': sel, 'genomes_used': s['genomes'],
        'total_proteins': s['proteins'], 'core_genes': s['core'],
        'bins_at_65': bins_at_65,
        'largest_bin': max(int(b['genomes']) for b in bins),
        'min_aai': '%.2f' % min_aai, 'mean_aai': '%.2f' % mean_aai,
        'pairs_ge_95': pairs95,
        'loo_max_core_gain': loo_gain, 'loo_max_core_gain_genome': loo_genome,
        'first_split_threshold': first_split,
        'bins_at_70': by_t[70], 'bins_at_75': by_t[75],
        'bins_at_80': by_t[80], 'bins_at_85': by_t[85],
        'boundary_flags_65': n65, 'boundary_flags_95': n95,
        'bin65_matches_pygemini': 'yes' if by_t[65] == bins_at_65 else
                                  'NO (%d recomputed vs %d reported)' % (by_t[65], bins_at_65),
        'wall_seconds': wall_seconds,
        'finished_utc': __import__('datetime').datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%SZ'),
    }
    row.update(cov)
    write_result_row(row)
    return row


def write_result_row(row):
    """Rewrite survey_results.tsv atomically, replacing this genus's row."""
    path = os.path.join(REPORTS, 'survey_results.tsv')
    existing = []
    if os.path.exists(path):
        with open(path, newline='') as fh:
            existing = [r for r in csv.DictReader(fh, delimiter='\t')
                        if r.get('genus') != row['genus']]
    tmp = path + '.tmp'
    with open(tmp, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=RESULT_COLS, delimiter='\t',
                           lineterminator='\n', extrasaction='ignore')
        w.writeheader()
        for r in existing:
            w.writerow({c: r.get(c, '') for c in RESULT_COLS})
        w.writerow({c: row.get(c, '') for c in RESULT_COLS})
    os.replace(tmp, path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('genus')
    ap.add_argument('--wall-seconds', default='')
    args = ap.parse_args()
    row = summarise(args.genus, args.wall_seconds)
    for c in RESULT_COLS:
        print('%-28s %s' % (c, row.get(c, '')))


if __name__ == '__main__':
    main()
