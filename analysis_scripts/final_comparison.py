#!/usr/bin/env python3
"""Final RaGCAn vs POCP comparison, all 1,160 genera. Run by finalise.sh when the POCP run ends."""
import argparse, csv, math, os, pathlib


def _data_root():
    """The directory holding POCP_run/ and LPSN_Genus_Survey/.

    Order: --root, then $RAGCAN_DATA_ROOT, then walk up from this file. No path is hardcoded.
    """
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root', metavar='DIR',
                    help='directory containing POCP_run/ and LPSN_Genus_Survey/')
    a = ap.parse_args()
    need = ('POCP_run', 'LPSN_Genus_Survey')

    def ok(d):
        return all((d / n).is_dir() for n in need)

    for how, d in (('--root', a.root), ('$RAGCAN_DATA_ROOT', os.environ.get('RAGCAN_DATA_ROOT'))):
        if d:
            d = pathlib.Path(d).resolve()
            if not ok(d):
                raise SystemExit(f'{how} is {d}, but it does not contain both {need[0]}/ and {need[1]}/.')
            return d
    here = pathlib.Path(__file__).resolve()
    for d in (here.parent, *here.parents):
        if ok(d):
            return d
    raise SystemExit('Could not find the result tables. Pass --root, or set RAGCAN_DATA_ROOT, to the\n'
                     'directory that contains POCP_run/ and LPSN_Genus_Survey/.')


R = _data_root()
pocp = {r['genus']: int(r['bins']) for r in
        csv.DictReader(open(R/'POCP_run'/'reports'/'pocp_status.tsv'), delimiter='\t')
        if r['state'] == 'done' and r['bins']}
ag = {r['genus']: r for r in
      csv.DictReader(open(R/'LPSN_Genus_Survey'/'phase3'/'agreement.tsv'), delimiter='\t')}
common = sorted(set(pocp) & set(ag))

def wilson(k, n, z=1.959963985):
    if n == 0: return (float('nan'),)*3
    p = k/n; d = 1+z*z/n; c = (p+z*z/(2*n))/d
    h = z*math.sqrt(p*(1-p)/n + z*z/(4*n*n))/d
    return p*100, max(0,(c-h))*100, min(1,(c+h))*100

def kappa(tp, fp, fn, tn):
    n = tp+fp+fn+tn; po = (tp+tn)/n
    pe = ((tp+fp)*(tp+fn)+(fn+tn)*(fp+tn))/(n*n)
    return (po-pe)/(1-pe)

print(f"FINAL COMPARISON — {len(common)} genera with both results")
print(f"POCP run complete: {len(pocp)} genera\n")
res = {}
print(f"{'method':<22}{'flag':>6}{'precision (95% CI)':>24}{'specificity (95% CI)':>24}{'recall (95% CI)':>22}{'kappa':>8}")
print('-'*108)
for lab, get in (('RaGCAn (AAI 65%)', lambda g: int(ag[g]['bins'])),
                 ('POCP (published)',  lambda g: pocp[g])):
    TP = sum(1 for g in common if get(g) > 1 and int(ag[g]['gtdb_genera']) > 1)
    FP = sum(1 for g in common if get(g) > 1 and int(ag[g]['gtdb_genera']) == 1)
    FN = sum(1 for g in common if get(g) == 1 and int(ag[g]['gtdb_genera']) > 1)
    TN = sum(1 for g in common if get(g) == 1 and int(ag[g]['gtdb_genera']) == 1)
    pr, sp, rc = wilson(TP, TP+FP), wilson(TN, TN+FP), wilson(TP, TP+FN)
    k = kappa(TP, FP, FN, TN); res[lab] = (pr, sp, rc, k)
    print(f"{lab:<22}{TP+FP:>6}{pr[0]:>9.1f} ({pr[1]:4.1f}-{pr[2]:4.1f}){sp[0]:>10.1f} ({sp[1]:4.1f}-{sp[2]:4.1f})"
          f"{rc[0]:>9.1f} ({rc[1]:4.1f}-{rc[2]:4.1f}){k:>8.3f}")
A, B = res['RaGCAn (AAI 65%)'], res['POCP (published)']
print()
for name, i in (('precision', 0), ('specificity', 1), ('recall', 2)):
    a, b = A[i], B[i]
    sep = 'SEPARATE - a real difference' if (a[1] > b[2] or b[1] > a[2]) else 'OVERLAP - not distinguishable'
    print(f"  {name:<12} RaGCAn {a[0]:5.1f}  POCP {b[0]:5.1f}   diff {a[0]-b[0]:+5.1f}   {sep}")
agree = sum(1 for g in common if (int(ag[g]['bins']) > 1) == (pocp[g] > 1))
print(f"\n  the two agree on {agree}/{len(common)} genera ({agree/len(common)*100:.1f}%)")
print("\n  ⚠ If the intervals OVERLAP, the manuscripts must say 'comparable', never 'better'.")
print("    Both currently say comparable. If any interval separates, that sentence can change.")
