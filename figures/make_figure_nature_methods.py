#!/usr/bin/env python3
"""
Nature Methods version. Third file; the NAR short version and the long version are both kept.

BUILT TO NATURE PORTFOLIO FIGURE RULES, WHICH DIFFER FROM NAR'S:
  - 180 mm maximum width (NAR double column is 178 mm)
  - 5-7 pt sans-serif labels. NAR/OUP wanted LARGE text for a graphical abstract shown as a
    thumbnail; Nature wants small text because the figure is reproduced at print size.
  - >= 300 dpi, labels left as live text (vector PDF), never flattened
  - NO TITLE INSIDE THE FIGURE. Nature puts that in the caption. Panels are lettered a, b.

Most Nature Portfolio journals do not accept graphical abstracts at all (only Nature Chemical
Biology and Nature Chemistry publish them). Nature Methods' own policy could not be confirmed.
So this is built as a two-panel SUMMARY FIGURE, which is the form that journal does use.
The caption is printed to stdout, to go in the manuscript rather than in the image.

Every number is read from a result table. None typed in.
"""
import csv, pathlib, math
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

def _data_root():
    """Locate the directory holding the result tables.

    Order: --root on the command line, then $RAGCAN_DATA_ROOT, then walk up from this
    file looking for LPSN_Genus_Survey/. No path is hardcoded, so the script runs on any
    machine that has the result tree.
    """
    import argparse, os
    ap = argparse.ArgumentParser(add_help=True, description=__doc__)
    ap.add_argument('--root', metavar='DIR',
                    help='directory containing LPSN_Genus_Survey/ and the result folders')
    known, _ = ap.parse_known_args()

    def _check(d, how):
        d = pathlib.Path(d).resolve()
        if not (d / 'LPSN_Genus_Survey').is_dir():
            raise SystemExit(f'{how} is {d}, but it contains no LPSN_Genus_Survey/ '
                             'and so is not the result tree.')
        return d

    if known.root:
        return _check(known.root, '--root')
    env = os.environ.get('RAGCAN_DATA_ROOT')
    if env:
        return _check(env, '$RAGCAN_DATA_ROOT')
    here = pathlib.Path(__file__).resolve()
    for d in (here.parent, *here.parents):
        if (d / 'LPSN_Genus_Survey').is_dir():
            return d
    raise SystemExit(
        'Could not find the result tables.\n'
        'Pass --root /path/to/project, or set RAGCAN_DATA_ROOT, where that directory\n'
        'contains LPSN_Genus_Survey/ and the per-analysis result folders.')


R = _data_root()
OUT = pathlib.Path(__file__).resolve().parent
GREY, BLUE, ORANGE = '#4d4d4d', '#2c6fa8', '#c2703a'
INK, MUTED = '#1a1a1a', '#5f5f5f'

# ---- panel a: the Pseudomonas worked example
pr = R/'Stutzerimonas_Validation'/'results'/'pseudomonas_389'/'reports'
rows = list(csv.reader(open(pr/'aai_matrix.tsv'), delimiter='\t')); ids = rows[0][1:]
M = np.array([[float(x) if x not in ('', 'NA', 'nan') else np.nan for x in r[1:]]
              for r in rows[1:]], dtype=float)
bins = [r['members'].split() for r in csv.DictReader(open(pr/'proposed_bins.tsv'), delimiter='\t')]
bins.sort(key=len, reverse=True)
pos = {g: i for i, g in enumerate(ids)}
IA = [pos[g] for g in bins[0]]; IB = [pos[g] for g in bins[1]]
NBIG, NSMALL = len(IA), len(IB)
W_IN = np.concatenate([M[np.ix_(IA, IA)][np.triu_indices(NBIG, 1)],
                       M[np.ix_(IB, IB)][np.triu_indices(NSMALL, 1)]])
W_BT = M[np.ix_(IA, IB)].ravel()
W_IN = W_IN[np.isfinite(W_IN)]; W_BT = W_BT[np.isfinite(W_BT)]

# ---- panel b: survey performance, with Wilson 95% intervals
ag = list(csv.DictReader(open(R/'LPSN_Genus_Survey'/'phase3'/'agreement.tsv'), delimiter='\t'))
for r in ag:
    r['b'] = int(r['bins']); r['g'] = int(r['gtdb_genera'])
both_c = sum(1 for r in ag if r['b'] == 1 and r['g'] == 1)
both_s = sum(1 for r in ag if r['b'] > 1 and r['g'] > 1)
lump   = sum(1 for r in ag if r['b'] == 1 and r['g'] > 1)
fp     = sum(1 for r in ag if r['b'] > 1 and r['g'] == 1)

def wilson(k, n, z=1.96):
    p = k/n; d = 1 + z*z/n
    c = (p + z*z/(2*n))/d
    h = z*math.sqrt(p*(1-p)/n + z*z/(4*n*n))/d
    return 100*p, 100*(c-h), 100*(c+h)

METRICS = [('Flag corroborated\nby GTDB',      both_s, both_s+fp,   BLUE),
           ('Coherent genera\nleft undivided', both_c, both_c+fp,   BLUE),
           ('Divided genera\ndetected',        both_s, both_s+lump, ORANGE)]
print("values, each from LPSN_Genus_Survey/phase3/agreement.tsv:")
for lab, k, n, _ in METRICS:
    p, lo, hi = wilson(k, n)
    print(f"  {lab.replace(chr(10),' '):<34} {k}/{n} = {p:.1f}% [{lo:.1f}-{hi:.1f}]")
print(f"  Pseudomonas bins (proposed_bins.tsv) {NBIG} + {NSMALL}")

# ---- figure
MM = 25.4
plt.rcParams.update({'figure.dpi': 300, 'savefig.bbox': 'tight',
                     'font.family': 'Nimbus Sans', 'font.size': 7,
                     'axes.linewidth': 0.5, 'xtick.major.width': 0.5,
                     'ytick.major.width': 0.5})
fig = plt.figure(figsize=(180/MM, 62/MM))

# panel a
axa = fig.add_axes([0.055, 0.20, 0.40, 0.66])
edges = np.arange(60, 101, 0.5)
for data, colour in ((W_IN, BLUE), (W_BT, ORANGE)):
    h, _ = np.histogram(data, bins=edges)
    axa.fill_between(edges[:-1], 0, h/h.max(), step='post', color=colour, alpha=0.85, lw=0)
axa.plot([65, 65], [0, 1.03], color=GREY, ls=(0, (2.5, 1.8)), lw=0.8)
axa.text(65.6, 0.42, '65%', fontsize=6, color=GREY, ha='left')
axa.text(60.3, 1.12, f'between bins\n$n$ = {W_BT.size:,}', fontsize=6, color=ORANGE,
         ha='left', va='center', linespacing=1.3)
axa.text(89, 1.12, f'within a bin\n$n$ = {W_IN.size:,}', fontsize=6, color=BLUE,
         ha='center', va='center', linespacing=1.3)
axa.set_xlim(60, 100); axa.set_ylim(0, 1.30)
axa.set_yticks([]); axa.set_xticks([60, 70, 80, 90, 100])
axa.set_xlabel('Amino-acid identity (%)', fontsize=7, color=INK, labelpad=2)
axa.tick_params(labelsize=6.5, colors=INK, length=2, pad=1.5)
for s in ('top', 'right', 'left'): axa.spines[s].set_visible(False)
axa.spines['bottom'].set_color(INK)

# panel b
axb = fig.add_axes([0.62, 0.20, 0.33, 0.66])
ys = np.arange(len(METRICS))[::-1]
for y, (lab, k, n, colour) in zip(ys, METRICS):
    p, lo, hi = wilson(k, n)
    axb.barh(y, p, height=0.52, color=colour, alpha=0.85, lw=0)
    axb.plot([lo, hi], [y, y], color=INK, lw=0.7, solid_capstyle='butt')
    for e in (lo, hi):
        axb.plot([e, e], [y-0.10, y+0.10], color=INK, lw=0.7)
    axb.text(hi+2.5, y, f'{p:.1f}%', va='center', fontsize=6.5, color=INK)
axb.set_yticks(ys); axb.set_yticklabels([m[0] for m in METRICS], fontsize=6.5, color=INK,
                                        linespacing=1.3)
axb.set_xlim(0, 118); axb.set_xticks([0, 25, 50, 75, 100])
axb.set_xlabel('Genera (%), Wilson 95% CI', fontsize=7, color=INK, labelpad=2)
axb.tick_params(labelsize=6.5, colors=INK, length=2, pad=1.5)
for s in ('top', 'right', 'left'): axb.spines[s].set_visible(False)
axb.spines['bottom'].set_color(INK)

# panel letters, Nature house style
fig.text(0.012, 0.95, 'a', fontsize=8, fontweight='bold', color=INK, va='top')
fig.text(0.545, 0.95, 'b', fontsize=8, fontweight='bold', color=INK, va='top')

fig.savefig(OUT/'RaGCAn_figure_nature_methods.pdf')
fig.savefig(OUT/'RaGCAn_figure_nature_methods.png', dpi=600)
plt.close(fig)

print(f"\nwrote {OUT/'RaGCAn_figure_nature_methods.pdf'}  (vector, live text)")
print(f"wrote {OUT/'RaGCAn_figure_nature_methods.png'}  (600 dpi)")
print("\n--- CAPTION, for the manuscript and NOT inside the image ---")
print(f"""Fig. 1 | A core-genome coherence screen for named prokaryotic genera.
a, Distribution of pairwise amino-acid identity (AAI) among {NBIG+NSMALL} Pseudomonas genomes,
separated into pairs falling within a bin proposed by RaGCAn (blue) and pairs falling between the
two proposed bins (orange). No taxonomic labels were supplied. Every between-bin pair lies below
the 65% default threshold (dashed line); the screen returned {NBIG} genomes and {NSMALL}, and all
{NSMALL} have since been reassigned to Halopseudomonas. b, Agreement with GTDB across
{len(ag):,} genera. Bars are point estimates, whiskers Wilson 95% confidence intervals. The screen
detects a minority of the genera GTDB divides, so a genus it does not flag has not been checked.""")
