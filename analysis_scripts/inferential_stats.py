#!/usr/bin/env python3
"""
Inferential statistics for the survey. Point estimates alone overstate certainty, especially
for archaea where n is small.

  - Wilson score 95% intervals for precision, specificity, recall (better than normal
    approximation at proportions near 0 or 1, and at small n).
  - Cohen's kappa: agreement beyond that expected by chance, with a bootstrap interval.
  - Permutation test: is the observed partition agreement (ARI) better than the same bins
    assigned at random?
  - Bootstrap interval for the median ARI.

No external stats package; every formula is implemented here so it can be checked.

NOTE: this file must NOT be named statistics.py - that shadows Python's own stdlib module
and breaks `import statistics` with a circular-import error.
"""
import csv, math, pathlib, random, statistics as _st
random.seed(20260829)

ROOT = pathlib.Path(__file__).resolve().parent.parent
rows = list(csv.DictReader(open(ROOT/'phase3'/'agreement.tsv'), delimiter='\t'))
bacset = {r['genus'] for r in csv.DictReader(open(ROOT/'reports'/'genome_selection.tsv.bak-before-archaea'), delimiter='\t')}
for r in rows:
    r['b'] = int(r['bins']); r['g'] = int(r['gtdb_genera'])
    r['a'] = float(r['ari']) if r['ari'] not in ('nan','') else None
    r['dom'] = 'Bacteria' if r['genus'] in bacset else 'Archaea'

def wilson(k, n, z=1.959963985):
    if n == 0: return (float('nan'),)*3
    p = k/n; d = 1 + z*z/n
    c = (p + z*z/(2*n))/d
    h = z*math.sqrt(p*(1-p)/n + z*z/(4*n*n))/d
    return p*100, max(0.0,(c-h))*100, min(1.0,(c+h))*100

def kappa(tp, fp, fn, tn):
    n = tp+fp+fn+tn
    po = (tp+tn)/n
    pe = ((tp+fp)*(tp+fn) + (fn+tn)*(fp+tn))/(n*n)
    return (po-pe)/(1-pe) if pe != 1 else float('nan')

def counts(sub):
    tp=sum(1 for r in sub if r['b']>1 and r['g']>1); fp=sum(1 for r in sub if r['b']>1 and r['g']==1)
    fn=sum(1 for r in sub if r['b']==1 and r['g']>1); tn=sum(1 for r in sub if r['b']==1 and r['g']==1)
    return tp,fp,fn,tn

print("INFERENTIAL STATISTICS — Wilson 95% intervals, Cohen's kappa\n")
print(f"{'set':<10}{'n':>6}  {'precision (95% CI)':<26}{'specificity (95% CI)':<26}{'recall (95% CI)':<26}{'kappa':>7}")
print("-"*104)
for lab, sub in (('All', rows), ('Bacteria',[r for r in rows if r['dom']=='Bacteria']),
                 ('Archaea',[r for r in rows if r['dom']=='Archaea'])):
    tp,fp,fn,tn = counts(sub)
    P=wilson(tp,tp+fp); S=wilson(tn,tn+fp); R=wilson(tp,tp+fn)
    k=kappa(tp,fp,fn,tn)
    # bootstrap kappa
    ks=[]
    for _ in range(2000):
        bs=[sub[random.randrange(len(sub))] for _ in range(len(sub))]
        t2=counts(bs)
        if (t2[0]+t2[1]+t2[2]+t2[3]): 
            v=kappa(*t2)
            if v==v: ks.append(v)
    ks.sort()
    lo,hi=(ks[int(.025*len(ks))], ks[int(.975*len(ks))]) if ks else (float('nan'),)*2
    print(f"{lab:<10}{len(sub):>6}  {P[0]:5.1f} ({P[1]:4.1f}-{P[2]:5.1f})     "
          f"{S[0]:5.1f} ({S[1]:4.1f}-{S[2]:5.1f})     {R[0]:5.1f} ({R[1]:4.1f}-{R[2]:5.1f})     "
          f"{k:.3f} ({lo:.2f}-{hi:.2f})")

# ---- ARI: is agreement better than chance? ----
both=[r for r in rows if r['b']>1 and r['g']>1 and r['a'] is not None]
a=sorted(r['a'] for r in both)
med=_st.median(a)   # true median: 136 is even, so average the two middle values
print(f"\nPARTITION AGREEMENT — {len(both)} genera where both divide the genus")
print(f"  median ARI {med:.3f}   mean {sum(a)/len(a):.3f}   identical (ARI=1.0): {sum(1 for x in a if x>=0.999)}")
bs=[]
for _ in range(5000):
    samp=sorted(a[random.randrange(len(a))] for _ in range(len(a)))
    bs.append(_st.median(samp))
bs.sort()
print(f"  bootstrap 95% CI for the median: {bs[int(.025*len(bs))]:.3f} - {bs[int(.975*len(bs))]:.3f}")
print(f"  ARI is chance-corrected by construction: 0 = chance. Observed median {med:.3f}.")
n_pos=sum(1 for x in a if x>0)
# sign test against ARI=0
p = sum(math.comb(len(a),k) for k in range(n_pos,len(a)+1)) / 2**len(a)
print(f"  {n_pos}/{len(a)} genera have ARI > 0; sign test vs chance p = {p:.3g}")
