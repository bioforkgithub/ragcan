#!/usr/bin/env python3
"""
THE DECISIVE CHECK. For each 'real gap' over-split genus, does the genome RaGCAn put in its own
bin match the genome GTDB assigns to a DIFFERENT genus?

If yes, RaGCAn is not over-splitting. It is agreeing with GTDB, and phase3/gtdb_compare.py hid
that by requiring a GTDB genus to hold >=2 genomes before it counted as a real split - which
discards exactly the one-genome reassignments this test is about.
"""
import csv, pathlib, re

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
gt = {}
for line in _taxonomy_lines():
    a,t=line.rstrip('\n').split('\t',1)
    if a[:3] in ('RS_','GB_'): a=a.split('_',1)[1]
    m=re.search(r'g__([^;]*)',t); gt[a]=m.group(1) if m else ''
def acc(l):
    m=re.search(r'(GC[AF]_\d+\.\d+)',l); return m.group(1) if m else None

hit=miss=nodata=0
print(f"{'genus':<20}{'RaGCAn singleton bin':<44}{'GTDB genus':<22}match?")
print("-"*104)
for r in csv.DictReader(open(ROOT/'phase3'/'oversplit_triage.tsv'), delimiter='\t'):
    if float(r['margin_below_65'])<=5.0: continue
    g=r['genus']; d=ROOT/'results'/g
    bins=[[x for x in re.split(r'[;\s]+',row['members'].strip()) if x]
          for row in csv.DictReader(open(d/'reports'/'proposed_bins.tsv'),delimiter='\t')]
    big=max(bins,key=len)
    major={gt.get(acc(x)) for x in big if gt.get(acc(x))}
    # the dominant GTDB genus among the LARGEST bin
    from collections import Counter
    cnt=Counter(gt.get(acc(x)) for x in big if gt.get(acc(x)))
    dom=cnt.most_common(1)[0][0] if cnt else None
    for b in bins:
        if len(b)!=1: continue
        a=acc(b[0]); gg=gt.get(a)
        if gg is None or gg=='':
            print(f"{g:<20}{b[0][:42]:<44}{'(not in GTDB)':<22}?"); nodata+=1
        elif gg!=dom:
            print(f"{g:<20}{b[0][:42]:<44}{gg:<22}YES - GTDB agrees it is elsewhere"); hit+=1
        else:
            print(f"{g:<20}{b[0][:42]:<44}{gg:<22}no - GTDB keeps it here"); miss+=1
print("-"*104)
tot=hit+miss+nodata
print(f"\n  {hit}/{tot} singleton bins are genomes GTDB ALSO places in a different genus")
print(f"  {miss}/{tot} GTDB keeps in the same genus (genuine disagreement)")
print(f"  {nodata}/{tot} not in GTDB at all (cannot say)")
