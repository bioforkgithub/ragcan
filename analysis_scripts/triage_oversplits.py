#!/usr/bin/env python3
"""
Triage the 'RaGCAn split, GTDB did not' genera: threshold noise, or a real gap?

PLAN.md Phase 3 names this ambiguity - "candidate false positive, OR a real finding GTDB has
not made" - and does not resolve it. This resolves the cheap half of it, using only files
already on disk. Nothing is re-run.

THE TEST. RaGCAn merges two bins only if EVERY cross pair is >= 65% AAI (complete linkage).
So for a genus it split, there is a best rejected merge: the bin pair whose weakest cross pair
sits highest, but still below 65. That value is the MARGIN.

  margin 64.9  -> the split is a coin flip. One pair, 0.1 points below an arbitrary line.
  margin 55.0  -> a genuine 10-point gap. RaGCAn is reporting real structure.

A margin close to 65 does NOT prove RaGCAn wrong - it proves the CALL WAS ARBITRARY, which is
the honest thing to report and is exactly the threshold argument. A wide margin means GTDB and
RaGCAn genuinely disagree about the organisms, and that case needs a human.
"""
import csv, pathlib, re, sys, argparse
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parent.parent
THRESH = 65.0

def load_matrix(p):
    rows = list(csv.reader(open(p), delimiter='\t'))
    ids = rows[0][1:]
    m = np.array([[float(x) if x not in ('', 'NA', 'nan') else np.nan for x in r[1:]]
                  for r in rows[1:]], dtype=float)
    return ids, m

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--per-genus', default='phase3/per_genus.tsv')
    ap.add_argument('--out', default='phase3/oversplit_triage.tsv')
    a = ap.parse_args()

    targets = [r for r in csv.DictReader(open(a.per_genus), delimiter='\t')
               if r['outcome'] == 'split/coherent RaGCAn only']
    out = []
    for r in targets:
        d = ROOT / 'results' / r['genus']
        mp, bp = d/'reports'/'aai_matrix.tsv', d/'reports'/'proposed_bins.tsv'
        if not (mp.is_file() and bp.is_file()):
            continue
        ids, m = load_matrix(mp)
        idx = {g: i for i, g in enumerate(ids)}
        bins = []
        for row in csv.DictReader(open(bp), delimiter='\t'):
            members = [idx[x] for x in re.split(r'[;\s]+', row['members'].strip())
                       if x and x in idx]
            if members:
                bins.append(members)
        if len(bins) < 2:
            continue
        best = -np.inf                      # highest rejected complete-linkage value
        for i in range(len(bins)):
            for j in range(i+1, len(bins)):
                sub = m[np.ix_(bins[i], bins[j])]
                if np.all(np.isnan(sub)):
                    continue
                link = np.nanmin(sub)       # complete linkage = worst cross pair
                best = max(best, link)
        if not np.isfinite(best):
            continue
        out.append(dict(genus=r['genus'], genomes=r['genomes'], bins=r['bins'],
                        best_rejected_linkage=round(float(best), 2),
                        margin_below_65=round(THRESH - float(best), 2),
                        gtdb_dominant_pct=r['dominant_pct']))

    out.sort(key=lambda x: x['margin_below_65'])
    with open(a.out, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]), delimiter='\t', lineterminator='\n')
        w.writeheader(); w.writerows(out)

    coin = [r for r in out if r['margin_below_65'] <= 1.0]
    near = [r for r in out if 1.0 < r['margin_below_65'] <= 5.0]
    real = [r for r in out if r['margin_below_65'] > 5.0]
    n = len(out)
    print(f"Triage of {n} 'RaGCAn split, GTDB did not' genera\n")
    print(f"  {len(coin):>4}  {len(coin)/n*100:>5.1f}%   margin <= 1.0 pt   ARBITRARY - the split is a coin flip at the threshold")
    print(f"  {len(near):>4}  {len(near)/n*100:>5.1f}%   margin 1-5 pt      marginal")
    print(f"  {len(real):>4}  {len(real)/n*100:>5.1f}%   margin > 5.0 pt    REAL GAP - needs a human; may be a finding GTDB has not made")
    print(f"\n  the 10 with the WIDEST gaps (look at these first):")
    for r in sorted(out, key=lambda x: -x['margin_below_65'])[:10]:
        print(f"    {r['genus']:<24} {r['genomes']:>4} genomes, {r['bins']:>2} bins, "
              f"best rejected linkage {r['best_rejected_linkage']:>6}  (gap {r['margin_below_65']:>5} pt)")
    print(f"\n  full table -> {a.out}")

if __name__ == '__main__':
    main()
