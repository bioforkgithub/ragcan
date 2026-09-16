#!/usr/bin/env python3
"""
PHASE 3 — RaGCAn bins vs GTDB genus assignments. PLAN.md: "this is what makes it a paper."

Not a score against LPSN: PLAN.md Phase 7 registers that agreement with existing taxonomy is
the WRONG target, since a tool that reproduces the taxonomy cannot claim the taxonomy is wrong.
This asks the right question instead - when RaGCAn disagrees with the current genus, does GTDB
disagree in the same place?

A GTDB genus counts as REAL only if it holds >= MIN_MEMBERS of that genus's genomes. Singletons
are usually one misfiled genome, not a genuine split, and counting them inflates "GTDB splits it".

Reads only files already on disk. Nothing is downloaded and nothing is re-run.
  python3 phase3/gtdb_compare.py [--min-members 2] [--csv out.tsv]
"""
import csv, pathlib, re, sys, argparse
from collections import Counter

def _taxonomy_lines():
    """Yield lines from BOTH GTDB taxonomy files.

    bac120_taxonomy.tsv is bacteria only. The 65 archaeal genera added 2026-08-27 live in
    ar53_taxonomy.tsv, and omitting it silently dropped every archaeal genome from the
    comparison - 0 of 65 genera scored, with no error. Both files have the same
    'accession<TAB>d__...;g__...;s__...' format.
    """
    import pathlib as _pl
    for name in ('bac120_taxonomy.tsv', 'ar53_taxonomy.tsv'):
        f = _pl.Path(ROOT)/'reference'/name
        if f.is_file():
            with open(f) as fh:
                for line in fh:
                    yield line



ROOT = pathlib.Path(__file__).resolve().parent.parent
MIN_MEMBERS = 2

def load_gtdb():
    g = {}
    with _taxonomy_lines() as fh:
        for line in fh:
            acc, tax = line.rstrip('\n').split('\t', 1)
            if acc[:3] in ('RS_', 'GB_'):
                acc = acc.split('_', 1)[1]
            m = re.search(r'g__([^;]*)', tax)
            g[acc] = m.group(1) if m else ''
    return g

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--min-members', type=int, default=MIN_MEMBERS)
    ap.add_argument('--csv')
    a = ap.parse_args()
    gtdb = load_gtdb()
    rows = []
    for d in sorted((ROOT/'results').iterdir()):
        pb = d/'reports'/'proposed_bins.tsv'
        if not pb.is_file():
            continue
        bins = [[m for m in re.split(r'[;\s]+', r['members'].strip()) if m]
                for r in csv.DictReader(open(pb), delimiter='\t')]
        if not bins:
            continue
        seen = Counter()
        for acc in (x for b in bins for x in b):
            k = re.search(r'(GC[AF]_\d+\.\d+)', acc)
            gg = gtdb.get(k.group(1)) if k else None
            if gg:
                seen[gg] += 1
        if not seen:
            continue
        real = sum(1 for n in seen.values() if n >= a.min_members)
        rs, gs = len(bins) > 1, real > 1
        out = ('split/split corroborated'      if rs and gs else
               'coherent/coherent agree'       if not rs and not gs else
               'split/coherent RaGCAn only'    if rs else
               'coherent/split too permissive')
        rows.append(dict(genus=d.name, genomes=sum(seen.values()), bins=len(bins),
                         gtdb_genera=len(seen), gtdb_real=real,
                         dominant_pct=round(max(seen.values())/sum(seen.values())*100, 1),
                         outcome=out))
    c = Counter(r['outcome'] for r in rows)
    tot = len(rows)
    print(f"PHASE 3 — {tot} genera scored against GTDB r2026-04 (min {a.min_members} members per GTDB genus)\n")
    for k, v in c.most_common():
        print(f"  {v:>5}  {v/tot*100:>5.1f}%   {k}")
    agree = c['coherent/coherent agree'] + c['split/split corroborated']
    print(f"\n  {agree:>5}  {agree/tot*100:>5.1f}%   AGREE WITH GTDB (either both coherent or both split)")
    print(f"  {c['coherent/split too permissive']:>5}  {c['coherent/split too permissive']/tot*100:>5.1f}%   "
          f"RaGCAn LUMPED where GTDB splits  <- the known 65% limitation")
    print(f"  {c['split/coherent RaGCAn only']:>5}  {c['split/coherent RaGCAn only']/tot*100:>5.1f}%   "
          f"RaGCAn SPLIT where GTDB does not <- candidate false positives")
    if a.csv:
        with open(a.csv, 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
            w.writeheader(); w.writerows(rows)
        print(f"\n  per-genus table -> {a.csv}")

if __name__ == '__main__':
    main()
