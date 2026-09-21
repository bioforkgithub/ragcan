#!/usr/bin/env python3
"""Score the eukaryote sets against PLAN.md. NCBI genus = first word of each file name.
Usage: python3 analyse.py [set ...]   (default: every set with finished results)"""
import csv, itertools, pathlib, importlib.util, sys, collections
import numpy as np
X = pathlib.Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('pg', X.parent / 'pygemini_work/PY-GEMINI/pygemini.py')
pg = importlib.util.module_from_spec(spec); spec.loader.exec_module(pg)

def load(s):
    r = list(csv.reader(open(X / s / 'run/reports/aai_matrix.tsv'), delimiter='\t'))
    ids = r[0][1:]
    M = np.array([[float(v) if v not in ('', 'NA', 'nan') else np.nan for v in row[1:]] for row in r[1:]])
    bins = list(csv.DictReader(open(X / s / 'run/reports/proposed_bins.tsv'), delimiter='\t'))
    return ids, M, bins
genus = lambda g: g.split('_')[0]
name = lambda g: g.split('__')[0].replace('_', ' ')
def ari(a, b):
    n = len(a); c2 = lambda x: x * (x - 1) // 2
    sij = sum(c2(v) for v in collections.Counter(zip(a, b)).values())
    sa = sum(c2(v) for v in collections.Counter(a).values()); sb = sum(c2(v) for v in collections.Counter(b).values())
    e = sa * sb / c2(n); mx = (sa + sb) / 2
    return float('nan') if mx == e else (sij - e) / (mx - e)

def sweep(s, ids, M):
    lab = [genus(g) for g in ids]
    _, trace = pg.bin_genomes(M, list(range(len(ids))), -np.inf)
    merges = [(h, a[0], b[0]) for h, a, b in trace]
    truth = frozenset(frozenset(i for i, l in enumerate(lab) if l == x) for x in set(lab))
    out = []
    for T in [50 + 0.5 * k for k in range(100)]:
        par = list(range(len(ids)))
        def f(x):
            while par[x] != x: par[x] = par[par[x]]; x = par[x]
            return x
        for h, i, j in merges:
            if h < T: break
            par[f(i)] = f(j)
        p = [f(i) for i in range(len(ids))]
        sets = frozenset(frozenset(i for i in range(len(ids)) if p[i] == k) for k in set(p))
        out.append((T, len(sets), ari(p, lab), sets == truth))
    (X / s / 'sweep.tsv').write_text('threshold\tbins\tari_vs_ncbi_genus\texact\n' +
        ''.join(f'{t}\t{b}\t{a:.4f}\t{int(e)}\n' for t, b, a, e in out))
    ok = [r[0] for r in out if r[3]]
    best = max(out, key=lambda r: -1 if r[2] != r[2] else r[2])
    print(f'  sweep vs NCBI genus ({len(set(lab))} genera, {len(ids)} genomes): exact at '
          f'{(str(ok[0]) + "-" + str(ok[-1]) + "%") if ok else "no threshold"}; best ARI {best[2]:.3f} at {best[0]}% ({best[1]} bins)')

sets = sys.argv[1:] or [p.name for p in X.iterdir() if (p / 'run/reports/proposed_bins.tsv').is_file()]
for s in sets:
    ids, M, bins = load(s); ix = {g: i for i, g in enumerate(ids)}; A = lambda a, b: M[ix[a], ix[b]]
    print(f'\n=== {s}: {len(ids)} genomes, {len(bins)} bin(s) at 65%')
    for b in bins:
        print(f'  bin {b["bin"]} (min AAI {b["min_aai"]}): ' + ', '.join(name(m) for m in b['members'].split()))
    off = M[~np.eye(len(ids), dtype=bool)]
    print(f'  AAI range between genomes: {np.nanmin(off):.2f}-{np.nanmax(off):.2f}')
    if s == 'primates':
        hs = [g for g in ids if g.startswith('Homo_sapiens')]
        v = A(*hs); print(f'  E1 two Homo sapiens assemblies: {v:.2f} -> {"HOLDS" if v >= 99.5 else "FALSIFIED"}')
        pt = [g for g in ids if g.startswith('Pan_troglodytes')][0]; pp = [g for g in ids if g.startswith('Pan_paniscus')][0]
        h = hs[0]
        c1 = A(pt, pp) > A(pt, h) and A(pt, pp) > A(pp, h)
        nh = [g for g in ids if not g.startswith('Homo')]
        top = max(nh, key=lambda g: A(h, g))
        print(f'  E2 Pan-Pan {A(pt,pp):.2f}; Pan-Homo {A(pt,h):.2f} / {A(pp,h):.2f}; Homo nearest non-Homo: {name(top)} {A(h,top):.2f} '
              f'-> {"HOLDS" if c1 and top.startswith("Pan") else "FALSIFIED"}')
        for g in sorted(nh, key=lambda g: -A(h, g)): print(f'      Homo (GRCh38) vs {name(g):<32} {A(h, g):.2f}')
    if s == 'yeasts':
        ng = [g for g in ids if g.startswith('Nakaseomyces_glabratus')][0]
        sc = [g for g in ids if g.startswith('Saccharomyces_cerevisiae')][0]
        ca = [g for g in ids if g.startswith('Candida_albicans')][0]
        print(f'  E3 N. glabratus vs S. cerevisiae {A(ng,sc):.2f}, vs C. albicans {A(ng,ca):.2f} -> {"HOLDS" if A(ng,sc) > A(ng,ca) else "FALSIFIED"}')
        cand = {g for g in ids if g.startswith('Candida')}
        bad = [b['bin'] for b in bins if cand & set(b['members'].split()) and set(b['members'].split()) - cand]
        print(f'  E4 Candida share a 65% bin with a Saccharomycetaceae genome: {"yes, bin " + ",".join(bad) + " -> FALSIFIED" if bad else "no -> HOLDS"}')
        print('      Candida to nearest non-Candida: ' + '; '.join(f'{name(c)} {max(A(c,g) for g in ids if g not in cand):.2f}' for c in sorted(cand)))
    if s == 'drosophila':
        sl = [g for g in ids if g.startswith('Scaptodrosophila')][0]
        dr = [g for g in ids if g.startswith('Drosophila')]
        hi = max(A(sl, g) for g in dr); lo = min(A(a, b) for a, b in itertools.combinations(dr, 2))
        print(f'  E5 Scaptodrosophila max to Drosophila {hi:.2f}; min Drosophila-Drosophila {lo:.2f} -> {"HOLDS" if hi < lo else "FALSIFIED"}')
    if s == 'caenorhabditis':
        print(f'  E6 one bin at 65%: {"HOLDS" if len(bins) == 1 else "FALSIFIED"}')
    sweep(s, ids, M)
