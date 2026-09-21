#!/usr/bin/env python3
"""Download RefSeq proteomes for the four eukaryote sets and keep the longest protein per gene.

Selection follows PLAN.md. For each assembly: {asm}_protein.faa.gz and {asm}_feature_table.txt.gz
from the RefSeq FTP path; CDS rows of the feature table give product_accession -> GeneID and
product_length; the longest product per GeneID is kept. Writes <set>/proteomes/<Name>__<ACC>.faa
and prepare_log.tsv (accession, proteins in file, genes kept).
"""
import csv, gzip, io, pathlib, re, sys, time, urllib.request
X = pathlib.Path(__file__).resolve().parent
REF = X / 'reference'

def rows(div):
    for line in open(REF / f'as_{div}.txt'):
        if line.startswith('#'): continue
        f = line.rstrip('\n').split('\t')
        if len(f) >= 20 and f[10] == 'latest': yield f

SETS = {
 'yeasts': ('fungi', {'Saccharomyces','Kluyveromyces','Lachancea','Zygosaccharomyces','Torulaspora',
                      'Nakaseomyces','Kazachstania','Naumovozyma','Tetrapisispora','Vanderwaltozyma',
                      'Eremothecium','Candida'}),
 'caenorhabditis': ('invertebrate', {'Caenorhabditis'}),
 'drosophila': ('invertebrate', {'Drosophila','Scaptodrosophila'}),
 'primates': ('vertebrate_mammalian', {'Homo','Pan','Gorilla','Pongo','Hylobates','Nomascus',
                                      'Symphalangus','Macaca','Papio','Theropithecus','Chlorocebus'}),
}
RANK = {'reference genome': 0, 'representative genome': 1, 'na': 2}

def pick(div, genera):
    best = {}
    for f in rows(div):
        name = f[7].split()
        if name[0] not in genera or len(name) < 2 or name[1] in ('sp.', 'cf.', 'aff.'): continue
        acc = f[0]
        key = ' '.join(name[:2])
        if key == 'Homo sapiens':          # both human assemblies are kept on purpose (control E1)
            best[key + ' ' + acc] = f; continue
        cur = best.get(key)
        if cur is None or (RANK.get(f[4], 3), f[14]) < (RANK.get(cur[4], 3), cur[14]):
            best[key] = f
    return sorted(best.values(), key=lambda f: f[7])

def fetch(url, tries=4):
    for k in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=300) as r: return r.read()
        except Exception as e:
            if k == tries - 1: raise
            time.sleep(10 * (k + 1))

def main():
    log = open(X / 'prepare_log.tsv', 'a')
    for set_name, (div, genera) in SETS.items():
        out = X / set_name / 'proteomes'; out.mkdir(parents=True, exist_ok=True)
        chosen = pick(div, genera)
        print(f'{set_name}: {len(chosen)} assemblies', flush=True)
        for f in chosen:
            acc, ftp = f[0], f[19].replace('ftp://', 'https://').rstrip('/')
            base = ftp.rsplit('/', 1)[1]
            label = re.sub(r'[^A-Za-z0-9]+', '_', ' '.join(f[7].split()[:3])).strip('_')
            dest = out / f'{label}__{acc}.faa'
            if dest.exists() and dest.stat().st_size > 0:
                continue
            ft = gzip.decompress(fetch(f'{ftp}/{base}_feature_table.txt.gz')).decode()
            longest = {}
            for r in csv.reader(io.StringIO(ft), delimiter='\t'):
                if not r or r[0] != 'CDS' or len(r) < 19: continue
                prot, gene, plen = r[10], r[15], r[18]
                if not prot or not gene or not plen.isdigit(): continue
                if gene not in longest or int(plen) > longest[gene][1]:
                    longest[gene] = (prot, int(plen))
            keep = {p for p, _ in longest.values()}
            faa = gzip.decompress(fetch(f'{ftp}/{base}_protein.faa.gz')).decode()
            n_in = faa.count('>'); written = 0
            tmp = dest.with_suffix('.tmp')
            with open(tmp, 'w') as fh:
                on = False
                for line in faa.splitlines():
                    if line.startswith('>'):
                        on = line[1:].split()[0] in keep
                        written += on
                    if on: fh.write(line + '\n')
            tmp.rename(dest)
            log.write(f'{set_name}\t{acc}\t{f[7]}\t{n_in}\t{len(longest)}\t{written}\n'); log.flush()
            print(f'  {f[7][:45]:<45} {acc}  proteins {n_in:>7} -> genes kept {written:>6}', flush=True)
main()
