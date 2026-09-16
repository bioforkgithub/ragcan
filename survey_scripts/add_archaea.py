#!/usr/bin/env python3
"""
Append the archaeal genera to reports/genome_selection.tsv, using select_genomes.py's OWN
filters against assembly_summary_archaea.txt.

WHY IT IS SAFE TO RUN WHILE THE BACTERIAL SWEEP IS LIVE
  run_survey.py builds its work list ONCE at startup (`order = genus_order()`, line 311) and
  only READS genome_selection.tsv - it never writes it. So appending here cannot disturb the
  running driver, and the running driver will not pick the new genera up either. They are
  collected by the NEXT driver start, which the hourly cron guarantees.
  survey_status.tsv is deliberately NOT touched: the driver rewrites it atomically and a
  concurrent write could clobber an update.

  python3 scripts/add_archaea.py [--min-species 4] [--dry-run]
"""
import sys, os, csv, argparse, shutil
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import select_genomes as sg

ARCHAEA = os.path.join(sg.ROOT, 'reference', 'assembly_summary_archaea.txt')
COLS = ['genus','species','accession','organism_name','refseq_category',
        'assembly_level','protein_coding_gene_count','total_gene_count','ftp_path','faa_name']

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--min-species', type=int, default=4)
    ap.add_argument('--dry-run', action='store_true')
    a = ap.parse_args()

    picked = {}
    with open(ARCHAEA, newline='') as fh:
        for line in fh:
            if line.startswith('#'): continue
            f = line.rstrip('\n').split('\t')
            if len(f) <= sg.C_FTP: continue
            sp = f[sg.C_NAME].split(' ', 2)
            if len(sp) < 2: continue
            genus, epithet = sp[0], sp[1]
            if not genus[:1].isupper(): continue
            if f[sg.C_VER] != 'latest': continue
            if f[sg.C_CAT] not in sg.CATEGORIES: continue
            if not sg.EPITHET.match(epithet): continue
            if f[sg.C_FTP] in ('', 'na'): continue
            row = {'genus':genus,'species':'%s %s'%(genus,epithet),'accession':f[sg.C_ACC],
                   'organism_name':f[sg.C_NAME],'refseq_category':f[sg.C_CAT],
                   'assembly_level':f[sg.C_LEVEL],
                   'protein_coding_gene_count':f[sg.C_PROT_GENE] if len(f)>sg.C_PROT_GENE else '',
                   'total_gene_count':f[sg.C_TOTAL_GENE] if len(f)>sg.C_TOTAL_GENE else '',
                   'ftp_path':f[sg.C_FTP]}
            k=(genus,epithet); prev=picked.get(k)
            if prev is None or sg._rank(row) < sg._rank(prev): picked[k]=row

    per=defaultdict(list)
    for (g,_),r in picked.items(): per[g].append(r)
    keep={g:v for g,v in per.items() if len(v)>=a.min_species}
    rows=sorted((r for v in keep.values() for r in v), key=lambda r:(r['genus'],r['species']))
    for r in rows: r['faa_name']='%s__%s.faa'%(sg.safe_name(r['organism_name']),r['accession'])

    existing={r['genus'] for r in csv.DictReader(open(sg.OUT),delimiter='\t')}
    clash=sorted({r['genus'] for r in rows} & existing)
    print(f"  archaeal genera with >= {a.min_species} species : {len(keep)}")
    print(f"  genomes                                  : {len(rows)}")
    print(f"  genus-name clashes with bacteria         : {len(clash)}{' -> '+', '.join(clash) if clash else ''}")
    if clash:
        sys.exit("REFUSING: a clashing genus name would merge two different sets. Resolve first.")
    if a.dry_run:
        print("\n  --dry-run: nothing written"); return
    shutil.copy(sg.OUT, sg.OUT + '.bak-before-archaea')
    with open(sg.OUT,'a',newline='') as fh:
        w=csv.DictWriter(fh,fieldnames=COLS,delimiter='\t',lineterminator='\n')
        for r in rows: w.writerow({c:r[c] for c in COLS})
    tot=sum(1 for _ in csv.DictReader(open(sg.OUT),delimiter='\t'))
    print(f"\n  appended. genome_selection.tsv now has {tot} genomes")
    print(f"  backup: {os.path.basename(sg.OUT)}.bak-before-archaea")

if __name__ == '__main__':
    main()
