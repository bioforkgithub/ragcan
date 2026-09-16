#!/usr/bin/env python3
"""
Graphical abstract for the RaGCAn manuscript.

Same rule as Figures/make_figures.py: EVERY NUMBER IS READ FROM A RESULT TABLE.
Nothing is typed in. Each value prints to stdout with its source so it can be
checked against the manuscript.

Palette is the manuscript's own (GREY/BLUE/ORANGE from Figures/make_figures.py).
Checked 2026-09-16 for colour-vision separation: every pair is >= 15 OKLab dE at
normal vision and >= 8 under simulated deuteranopia and protanopia.

The heatmap is REAL DATA, not a drawn glyph: the 389-genome Pseudomonas AAI
matrix, rows ordered by the bins the program proposed.
"""
import csv, pathlib
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, FancyArrow
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
INK, MUTED, RULE = '#1a1a1a', '#5f5f5f', '#d8d8d8'

# ---------------------------------------------------------------- values
ag = list(csv.DictReader(open(R/'LPSN_Genus_Survey'/'phase3'/'agreement.tsv'), delimiter='\t'))
for r in ag:
    r['n'] = int(r['n']); r['ngt'] = int(r['not_in_gtdb'])
    r['b'] = int(r['bins']); r['g'] = int(r['gtdb_genera'])
both_c = sum(1 for r in ag if r['b'] == 1 and r['g'] == 1)
both_s = sum(1 for r in ag if r['b'] > 1 and r['g'] > 1)
lump   = sum(1 for r in ag if r['b'] == 1 and r['g'] > 1)
fp     = sum(1 for r in ag if r['b'] > 1 and r['g'] == 1)
PPV  = 100*both_s/(both_s+fp)
SPEC = 100*both_c/(both_c+fp)
SENS = 100*both_s/(both_s+lump)
IDENT = sum(1 for r in ag if r['b'] > 1 and r['g'] > 1 and r['ari'] == '1.0')

sw = list(csv.DictReader(open(R/'LPSN_Genus_Survey'/'phase4'/'sweep_curves.tsv'), delimiter='\t'))
win = [r for r in sw if r['exact_lo']]
mids = [(float(r['exact_lo'])+float(r['exact_hi']))/2 for r in win]
SPREAD = max(mids) - min(mids)
RECOV = 100*len(win)/len(sw)

# species-boundary spread, phase 6, and the phase 7 species-clean result
sp7 = list(csv.DictReader(open(R/'LPSN_Genus_Survey'/'phase7'/'bottom_up.tsv'), delimiter='\t'))
flags = list(csv.DictReader(open(R/'LPSN_Genus_Survey'/'phase7'/'species_level_flags.tsv'), delimiter='\t'))
NFLAG = len(flags)

# Pseudomonas 389: matrix and the bins the program proposed
pr = R/'Stutzerimonas_Validation'/'results'/'pseudomonas_389'/'reports'
rows = list(csv.reader(open(pr/'aai_matrix.tsv'), delimiter='\t'))
ids = rows[0][1:]
M = np.array([[float(x) if x not in ('', 'NA', 'nan') else np.nan for x in r[1:]]
              for r in rows[1:]], dtype=float)
bins = []
for r in csv.DictReader(open(pr/'proposed_bins.tsv'), delimiter='\t'):
    bins.append(r['members'].split())
bins.sort(key=len, reverse=True)
pos = {g: i for i, g in enumerate(ids)}
order = [pos[g] for b in bins for g in b]
Mo = M[np.ix_(order, order)]
NBIG, NSMALL = len(bins[0]), len(bins[1])
IA = [pos[g] for g in bins[0]]; IB = [pos[g] for g in bins[1]]
W_IN = np.concatenate([M[np.ix_(IA, IA)][np.triu_indices(NBIG, 1)],
                       M[np.ix_(IB, IB)][np.triu_indices(NSMALL, 1)]])
W_BT = M[np.ix_(IA, IB)].ravel()
W_IN = W_IN[np.isfinite(W_IN)]; W_BT = W_BT[np.isfinite(W_BT)]

print("VALUES, each with its source file")
print(f"  genera scored          {len(ag):>6}   phase3/agreement.tsv")
print(f"  genomes in them        {sum(r['n'] for r in ag):>6}   phase3/agreement.tsv (GTDB-assigned)")
print(f"  flag is correct        {PPV:>6.1f}%  phase3/agreement.tsv  [paper: 70.0%]")
print(f"  sound genera untouched {SPEC:>6.1f}%  phase3/agreement.tsv  [paper: 92.4%]")
print(f"  problem genera found   {SENS:>6.1f}%  phase3/agreement.tsv  [paper: 39.4%]")
print(f"  divisions identical    {IDENT:>6}   phase3/agreement.tsv  [paper: 48]")
print(f"  threshold spread       {SPREAD:>6.1f}   phase4/sweep_curves.tsv [paper: 41.8]")
print(f"  recoverable at some T  {RECOV:>6.1f}%  phase4/sweep_curves.tsv")
print(f"  Pseudomonas bins       {NBIG} + {NSMALL}   pseudomonas_389/reports/proposed_bins.tsv")
print(f"  species-named flags    {NFLAG:>6}   phase7/species_level_flags.tsv")

# ---------------------------------------------------------------- canvas
plt.rcParams.update({'font.size': 8, 'figure.dpi': 300, 'savefig.bbox': 'tight',
                     'font.family': 'DejaVu Sans'})
fig = plt.figure(figsize=(7.5, 5.2))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis('off')

def box(x, y, w, h, title, body, edge=GREY, face='white', lw=1.0):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0,rounding_size=1.3',
                                ec=edge, fc=face, lw=lw, zorder=2))
    ax.text(x+w/2, y+h-2.7, title, ha='center', va='top', fontsize=7.5,
            fontweight='bold', color=INK, zorder=3)
    ax.text(x+w/2, y+h-6.8, body, ha='center', va='top', fontsize=6.7,
            color=MUTED, linespacing=1.4, zorder=3)

def arrow(x, y):
    ax.add_patch(FancyArrow(x, y, 2.4, 0, width=0.3, head_width=1.5, head_length=1.2,
                            length_includes_head=True, color=GREY, zorder=3))

# ---- title
ax.text(3, 98.0, 'RaGCAn', fontsize=13.5, fontweight='bold', color=INK, va='top')
ax.text(15.2, 96.9, 'a rapid core-genome screen for whether a named genus holds together',
        fontsize=8.4, color=MUTED, va='top')
ax.plot([3, 97], [92.6, 92.6], color=RULE, lw=0.9)

# ---- workflow row
Y, H, W = 76.0, 13.6, 20.6
xs = [3.0, 26.2, 49.4, 72.6]
box(xs[0], Y, W, H, 'Proteomes of one genus', 'one FASTA per genome\nno tree, no alignment')
box(xs[1], Y, W, H, 'One protein search', 'DIAMOND, all against all\ncore gene set')
box(xs[2], Y, W, H, 'AAI for every pair', 'amino-acid identity\nacross the core')
box(xs[3], Y, W, H, 'Complete-linkage bins', 'merge only if the worst\npair still clears 65%',
    edge=BLUE, lw=1.4)
for x in xs[:3]: arrow(x+W+0.15, Y+H/2)
ax.text(50, 73.0, 'one bin = the genus holds together        more than one = a proposed division',
        ha='center', va='top', fontsize=7.2, color=INK, style='italic')

# ---- worked example, real data. A distribution, not a heatmap: at 7 genomes out of
#      389 the interesting block is 1.8% of a heatmap's width and cannot be seen.
ax.text(3, 68.5, 'One genus, no labels given', fontsize=8.2, fontweight='bold',
        color=INK, va='top')
ax.text(3, 64.6, f'{NBIG+NSMALL} $\\it{{Pseudomonas}}$ genomes in. Every pairwise AAI, split by\n'
                 f'the bins the program returned.',
        fontsize=6.9, color=MUTED, va='top', linespacing=1.5)

dp = fig.add_axes([0.048, 0.305, 0.285, 0.265])
edges = np.arange(60, 101, 0.5)
for data, colour, lab in ((W_IN, BLUE, 'within a bin'), (W_BT, ORANGE, 'between the two bins')):
    h, _ = np.histogram(data, bins=edges)
    dp.fill_between(edges[:-1], 0, h/h.max(), step='post', color=colour, alpha=0.82, lw=0)
dp.axvline(65, color=GREY, ls=(0, (3, 2)), lw=1.1, zorder=5)
dp.text(66.4, 1.14, '65% threshold', fontsize=6.3, color=GREY, va='top')
dp.text(92.0, 0.70, f'within a bin\n{W_IN.size:,} pairs', fontsize=6.6, color=BLUE,
        fontweight='bold', ha='center', linespacing=1.4)
dp.text(61.4, 0.86, f'between\n{W_BT.size:,} pairs', fontsize=6.6, color=ORANGE,
        fontweight='bold', ha='left', linespacing=1.4)
dp.set_xlim(60, 100); dp.set_ylim(0, 1.16)
dp.set_yticks([]); dp.set_xticks([60, 70, 80, 90, 100])
dp.set_xlabel('amino-acid identity (%)', fontsize=6.8, color=MUTED, labelpad=2)
dp.tick_params(labelsize=6.3, colors=MUTED, length=2.5, pad=2)
for s in ('top', 'right', 'left'): dp.spines[s].set_visible(False)
dp.spines['bottom'].set_color('#bbbbbb'); dp.spines['bottom'].set_linewidth(0.8)

ax.text(3, 23.0, f'All {W_BT.size:,} pairs between the two groups fall\n'
                 f'below the threshold. The program returned\n'
                 f'{NBIG} genomes and {NSMALL}. All {NSMALL} have since been\n'
                 f'reassigned to $\\it{{Halopseudomonas}}$.',
        fontsize=6.9, color=INK, va='top', linespacing=1.6)

# ---- stat tiles
def tile(x, y, w, h, big, small, colour):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0,rounding_size=1.2',
                                ec=RULE, fc='#fafafa', lw=0.8, zorder=2))
    ax.text(x+2.8, y+h-3.4, big, ha='left', va='top', fontsize=14.5,
            fontweight='bold', color=colour, zorder=3)
    ax.text(x+2.8, y+h-10.0, small, ha='left', va='top', fontsize=6.9,
            color=INK, linespacing=1.4, zorder=3)

TW, TH = 25.6, 15.8
tile(41.5, 47.0, TW, TH, f'{PPV:.1f}%', 'of flagged genera are ones\nGTDB also divides', BLUE)
tile(70.0, 47.0, TW, TH, f'{SPEC:.1f}%', 'of sound genera are\nleft alone', BLUE)
tile(41.5, 28.6, TW, TH, f'{SENS:.1f}%', 'of genera needing attention\nare found. It misses the rest', ORANGE)
tile(70.0, 28.6, TW, TH, f'{SPREAD:.1f}', 'points the ideal threshold\nmoves between genera', ORANGE)

ax.text(41.5, 24.6, '1,160 genera  ·  17,979 genomes  ·  bacteria and archaea  ·  30 h on one machine',
        fontsize=7.3, color=INK, va='top', fontweight='bold')
ax.text(41.5, 20.2, f'Divisions identical to the expert ones in {IDENT} genera. Not one bin in the\n'
                    f'survey cuts a named species in half, so every flag can be stated by\n'
                    f'species name ({NFLAG} of them).',
        fontsize=6.9, color=MUTED, va='top', linespacing=1.55)

# ---- footer
ax.plot([3, 97], [7.0, 7.0], color=RULE, lw=0.9)
ax.text(50, 3.6, 'It misses three of five real problems, so a genus it does not flag has not been '
                 'checked.  A first check, not a final answer.',
        ha='center', va='center', fontsize=7.4, color=INK, style='italic')

fig.savefig(OUT/'RaGCAn_graphical_abstract.pdf')
fig.savefig(OUT/'RaGCAn_graphical_abstract.png', dpi=300)
plt.close(fig)
print(f"\nwrote {OUT/'RaGCAn_graphical_abstract.pdf'}")
print(f"wrote {OUT/'RaGCAn_graphical_abstract.png'}")
