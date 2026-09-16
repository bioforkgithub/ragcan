#!/usr/bin/env python3
"""
PHASE 3, done the way PLAN.md line 115 actually specifies:
"do RaGCAn's bins align with GTDB's groups? (adjusted Rand index / cluster agreement)"

WHY THIS REPLACES phase3/gtdb_compare.py's 2x2 TABLE
The earlier script asked "does RaGCAn split? does GTDB split?" and required a GTDB genus to hold
>= 2 genomes before it counted as a real split. That rule DISCARDS single-genome reassignments -
and those turn out to be most of what RaGCAn finds. Measured 2026-08-27: of 15 singleton bins in
the wide-gap genera, 11 are genomes GTDB ALSO assigns to a different genus. The 2x2 scored those
as false positives when they are agreements.

Adjusted Rand index compares the two PARTITIONS directly and has no such blind spot.
  ARI = 1.0  identical partitions
  ARI = 0.0  no better than chance
  ARI < 0    worse than chance
Genomes absent from GTDB are excluded, per genus, and counted.
"""
import csv, pathlib, re, sys
from collections import Counter
from itertools import combinations

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

def ari(a, b):
    """Adjusted Rand index between two labelings of the same items."""
    n = len(a)
    if n < 2: return float('nan')
    tab = Counter(zip(a, b))
    ra, rb = Counter(a), Counter(b)
    c2 = lambda x: x*(x-1)//2
    sij = sum(c2(v) for v in tab.values())
    sa, sb = sum(c2(v) for v in ra.values()), sum(c2(v) for v in rb.values())
    exp = sa*sb/c2(n) if c2(n) else 0
    mx = (sa+sb)/2
    return float('nan') if mx == exp else (sij-exp)/(mx-exp)

gt = {}
for line in _taxonomy_lines():
    a, t = line.rstrip('\n').split('\t', 1)
    if a[:3] in ('RS_','GB_'): a = a.split('_',1)[1]
    m = re.search(r'g__([^;]*)', t); gt[a] = m.group(1) if m else ''

rows = []
for d in sorted((ROOT/'results').iterdir()):
    pb = d/'reports'/'proposed_bins.tsv'
    if not pb.is_file(): continue
    binlab, gtlab, drop = [], [], 0
    for k, row in enumerate(csv.DictReader(open(pb), delimiter='\t')):
        for x in re.split(r'[;\s]+', row['members'].strip()):
            if not x: continue
            m = re.search(r'(GC[AF]_\d+\.\d+)', x)
            g = gt.get(m.group(1)) if m else None
            if not g: drop += 1; continue
            binlab.append(k); gtlab.append(g)
    if len(binlab) < 3: continue
    rows.append(dict(genus=d.name, n=len(binlab), not_in_gtdb=drop,
                     bins=len(set(binlab)), gtdb_genera=len(set(gtlab)),
                     ari=round(ari(binlab, gtlab), 4),
                     identical=int(len(set(binlab))==1 and len(set(gtlab))==1)))

with open(ROOT/'phase3'/'agreement.tsv','w',newline='') as fh:
    w=csv.DictWriter(fh,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n')
    w.writeheader(); w.writerows(rows)

both1=[r for r in rows if r['bins']==1 and r['gtdb_genera']==1]
rest=[r for r in rows if not (r['bins']==1 and r['gtdb_genera']==1)]
print(f"PHASE 3 by cluster agreement — {len(rows)} genera\n")
print(f"  {len(both1):>5}  {len(both1)/len(rows)*100:>5.1f}%   BOTH say one genus (ARI undefined - perfect trivial agreement)")
print(f"  {len(rest):>5}  {len(rest)/len(rows)*100:>5.1f}%   at least one of them splits -> ARI computed below\n")
vals=[r['ari'] for r in rest if r['ari']==r['ari']]
if vals:
    vals.sort()
    band=lambda lo,hi:[v for v in vals if lo<=v<hi]
    print(f"  of those {len(vals)} informative genera:")
    for lo,hi,lab in [(0.75,1.01,'ARI >= 0.75   strong agreement'),
                      (0.40,0.75,'ARI 0.40-0.75 moderate'),
                      (0.10,0.40,'ARI 0.10-0.40 weak'),
                      (-1.0,0.10,'ARI < 0.10    little or none')]:
        b=band(lo,hi); print(f"    {len(b):>5}  {len(b)/len(vals)*100:>5.1f}%   {lab}")
    print(f"\n    median ARI {vals[len(vals)//2]:.3f}   mean {sum(vals)/len(vals):.3f}")
print(f"\n  table -> phase3/agreement.tsv")
