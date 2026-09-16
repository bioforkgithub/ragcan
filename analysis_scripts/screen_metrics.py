#!/usr/bin/env python3
"""
RaGCAn evaluated as a SCREEN, which is what Manish has always said it is.

His framing, 2026-08-27: "it was always planned as first line of identification in a taxon."
The manuscript says "first pass" three times. Under that framing, agreement/accuracy is the
wrong summary - a screen is judged on whether a FLAG is worth following, and on how much it
misses.

  flagged   = RaGCAn proposes more than one bin for the genus
  needs it  = GTDB places the genus's genomes in more than one genus

Genomes absent from GTDB are excluded per genus (counted in agreement.tsv).
Reads phase3/agreement.tsv, which is produced by phase3/agreement.py.
"""
import csv, pathlib, sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
src = ROOT/'phase3'/'agreement.tsv'
if not src.is_file():
    sys.exit("run phase3/agreement.py first")

rows = list(csv.DictReader(open(src), delimiter='\t'))
TP = sum(1 for r in rows if int(r['bins'])>1 and int(r['gtdb_genera'])>1)
FP = sum(1 for r in rows if int(r['bins'])>1 and int(r['gtdb_genera'])==1)
FN = sum(1 for r in rows if int(r['bins'])==1 and int(r['gtdb_genera'])>1)
TN = sum(1 for r in rows if int(r['bins'])==1 and int(r['gtdb_genera'])==1)
n = len(rows)
prec = TP/(TP+FP)*100 if TP+FP else float('nan')
rec  = TP/(TP+FN)*100 if TP+FN else float('nan')
spec = TN/(TN+FP)*100 if TN+FP else float('nan')

print(f"RaGCAn as a first-pass screen — {n} genera\n")
print(f"                        GTDB: >1 genus    GTDB: one genus")
print(f"  RaGCAn flags it       {TP:>10}       {FP:>13}")
print(f"  RaGCAn says fine      {FN:>10}       {TN:>13}\n")
print(f"  PRECISION   {prec:>5.1f}%   a flag is worth following up this often")
print(f"  SPECIFICITY {spec:>5.1f}%   clean genera it correctly leaves alone")
print(f"  RECALL      {rec:>5.1f}%   of the genera that need attention, it catches this many")
print(f"\n  READING: a rule-IN test, not a rule-OUT test. A flag means look;")
print(f"           'no flag' is weak evidence and must be reported as such.")
print(f"\n  ⚠ Do NOT summarise this as accuracy ({(TP+TN)/n*100:.1f}%). Accuracy is dominated by")
print(f"    the {TN} genera where both simply agree there is nothing to see, and it flatters the tool.")
