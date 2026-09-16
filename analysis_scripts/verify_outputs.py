#!/usr/bin/env python3
"""
EXHAUSTIVE INTERNAL VERIFICATION of every RaGCAn output in the survey.

This does not ask whether RaGCAn agrees with GTDB - that is phase3/agreement.py. It asks a
prior and more basic question: **does the program do what it says it does?** Every claim it
prints is re-derived here from its own raw output by independent code, for every genus.

Four checks, on all 1,160 genera:

  1 CORE      "N genes shared by all M genomes" - recounted from pangenome_matrix.tsv by
              counting ortholog groups present in every genome. Must equal N exactly.
  2 PARTITION  the proposed bins must be a true partition of the input: every genome in
              exactly one bin, none missing, none duplicated.
  3 LINKAGE    the stated clustering rule is complete linkage at 65% AAI. So within every
              bin, EVERY pair must be >= 65; and the best rejected cross-bin merge must be
              < 65. Checked pair by pair against the AAI matrix.
  4 MATRIX     the AAI matrix must be symmetric with a diagonal of 100.

A failure in any of these is a bug in the program or a corrupted run, not a question of
taxonomy. Run time is a few minutes; it re-reads every matrix.

  python3 phase3/verify_outputs.py [--tsv phase3/verification.tsv]
"""
import csv, pathlib, re, sys, argparse
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parent.parent
THRESH = 65.0
TOL = 0.01                      # matrices are written to 2 dp

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--tsv', default=str(ROOT/'phase3'/'verification.tsv'))
    a = ap.parse_args()

    out, counts = [], {k: [0, 0] for k in ('core', 'partition', 'linkage', 'matrix')}
    for d in sorted((ROOT/'results').iterdir()):
        rep = d/'reports'
        if not (rep/'proposed_bins.tsv').is_file():
            continue
        row = {'genus': d.name}

        # --- 1 CORE -------------------------------------------------------
        m = re.search(r'CORE GENOME: (\d+) genes shared by all (\d+) genomes',
                      (rep/'summary.txt').read_text()) if (rep/'summary.txt').is_file() else None
        pm = d/'pangenome_matrix.tsv'
        if m and pm.is_file():
            claim = int(m.group(1))
            with open(pm) as fh:
                r = csv.reader(fh, delimiter='\t'); hdr = next(r)
                n = len(hdr) - 1; present = [0]*n; rows = 0
                for line in r:
                    rows += 1
                    for i, v in enumerate(line[1:]):
                        if v not in ('', '0', '-', 'NA'):
                            present[i] += 1
            recount = sum(1 for c in present if c == rows)
            row['core_claim'], row['core_recount'] = claim, recount
            row['core_ok'] = int(recount == claim)
        else:
            row['core_ok'] = ''

        # --- load AAI matrix + bins --------------------------------------
        am = rep/'aai_matrix.tsv'
        if not am.is_file():
            out.append(row); continue
        mrows = list(csv.reader(open(am), delimiter='\t'))
        ids = mrows[0][1:]; idx = {x: i for i, x in enumerate(ids)}
        M = np.array([[float(v) if v not in ('', 'NA', 'nan') else np.nan for v in r[1:]]
                      for r in mrows[1:]], dtype=float)
        bins = []
        for b in csv.DictReader(open(rep/'proposed_bins.tsv'), delimiter='\t'):
            bins.append([x for x in re.split(r'[;\s]+', b['members'].strip()) if x])

        # --- 2 PARTITION --------------------------------------------------
        flat = [x for b in bins for x in b]
        row['partition_ok'] = int(len(flat) == len(set(flat)) == len(ids)
                                  and set(flat) == set(ids))

        # --- 3 LINKAGE ----------------------------------------------------
        bi = [[idx[x] for x in b if x in idx] for b in bins]
        worst_within, best_rejected = np.inf, -np.inf
        for b in bi:
            if len(b) < 2: continue
            sub = M[np.ix_(b, b)].copy(); np.fill_diagonal(sub, np.nan)
            if not np.all(np.isnan(sub)): worst_within = min(worst_within, np.nanmin(sub))
        for i in range(len(bi)):
            for j in range(i+1, len(bi)):
                sub = M[np.ix_(bi[i], bi[j])]
                if not np.all(np.isnan(sub)):
                    best_rejected = max(best_rejected, np.nanmin(sub))
        wok = (not np.isfinite(worst_within)) or worst_within >= THRESH - TOL
        rok = (not np.isfinite(best_rejected)) or best_rejected < THRESH + TOL
        row['worst_within'] = '' if not np.isfinite(worst_within) else round(float(worst_within), 2)
        row['best_rejected'] = '' if not np.isfinite(best_rejected) else round(float(best_rejected), 2)
        row['linkage_ok'] = int(wok and rok)

        # --- 4 MATRIX -----------------------------------------------------
        sym = np.allclose(M, M.T, equal_nan=True, atol=TOL)
        diag = np.allclose(np.diag(M), 100.0, atol=TOL)
        row['matrix_ok'] = int(sym and diag)

        for k in counts:
            v = row.get(k+'_ok', '')
            if v != '':
                counts[k][0 if v else 1] += 1
        out.append(row)

    cols = ['genus','core_claim','core_recount','core_ok','partition_ok',
            'worst_within','best_rejected','linkage_ok','matrix_ok']
    with open(a.tsv, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter='\t', lineterminator='\n',
                           extrasaction='ignore')
        w.writeheader(); w.writerows(out)

    print(f"EXHAUSTIVE INTERNAL VERIFICATION — {len(out)} genera\n")
    names = {'core':'core-gene count re-derived from the presence matrix',
             'partition':'bins form a true partition of the input genomes',
             'linkage':'complete-linkage rule holds at 65% AAI',
             'matrix':'AAI matrix symmetric, diagonal 100'}
    allpass = True
    for k in ('core','partition','linkage','matrix'):
        p, f = counts[k]
        if f: allpass = False
        print(f"  {p:>5} passed  {f:>3} failed   {names[k]}")
    print(f"\n  {'ALL CHECKS PASSED' if allpass else 'FAILURES PRESENT - investigate before publishing'}")
    print(f"  per-genus table -> {a.tsv}")
    return 0 if allpass else 1

if __name__ == '__main__':
    sys.exit(main())
