#!/usr/bin/env python3
"""
Derive the full-sweep genus list: every genus with >= MIN_SPECIES validly named species.

Uses select_genomes.py's OWN filters, imported rather than reimplemented, so the list
cannot drift from what the selector will actually accept. PLAN.md section 1 gives 1,095
genera / 17,622 genomes at >= 4 species; this script re-derives that number rather than
trusting it.

  python3 scripts/list_all_genera.py            # prints the genus list, one per line
  python3 scripts/list_all_genera.py --counts   # genus <tab> n_species
"""
import sys, os, re, argparse
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import select_genomes as sg          # reuse its constants and filters

MIN_SPECIES = 4

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--min-species', type=int, default=MIN_SPECIES)
    ap.add_argument('--counts', action='store_true')
    a = ap.parse_args()

    # one pass, applying exactly select_genomes.parse()'s per-row rules
    species = defaultdict(set)
    with open(sg.SUMMARY, newline='') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) <= sg.C_FTP:
                continue
            sp = f[sg.C_NAME].split(' ', 2)
            if len(sp) < 2:
                continue
            genus, epithet = sp[0], sp[1]
            if f[sg.C_VER] != 'latest':          continue
            if f[sg.C_CAT] not in sg.CATEGORIES: continue
            if not sg.EPITHET.match(epithet):    continue
            if f[sg.C_FTP] in ('', 'na'):        continue
            if not genus[:1].isupper():          continue
            species[genus].add(epithet)

    keep = {g: len(s) for g, s in species.items() if len(s) >= a.min_species}
    for g in sorted(keep):
        print(f"{g}\t{keep[g]}" if a.counts else g)
    print(f"\n{len(keep)} genera, {sum(keep.values())} genomes at >= {a.min_species} species",
          file=sys.stderr)

if __name__ == '__main__':
    main()
