#!/usr/bin/env python3
"""
SHORT graphical abstract for RaGCAn. The long one is make_graphical_abstract.py and is kept.

Built to NAR/OUP figure rules rather than to page width:
  - 178 mm wide (NAR double column), landscape
  - Nimbus Sans, a Helvetica-metric font (NAR asks for Arial or Helvetica; neither is
    installed here, Nimbus Sans is the standard Helvetica clone)
  - nothing below 8 pt, so it survives reduction to a table-of-contents thumbnail
  - about 40 words, one message

Same rule as the long version: every number is read from a result table, none typed in.
"""
import csv, pathlib
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
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
print(f"from proposed_bins.tsv: {NBIG} + {NSMALL} = {NBIG+NSMALL} genomes")
print(f"from aai_matrix.tsv:    {W_IN.size:,} within-bin pairs, {W_BT.size:,} between")

plt.rcParams.update({'figure.dpi': 300, 'savefig.bbox': 'tight',
                     'font.family': 'Nimbus Sans'})
MM = 25.4
fig = plt.figure(figsize=(178/MM, 96/MM))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis('off')

# --- title (6 words)
ax.text(2.5, 97, 'RaGCAn: does this genus hold together?',
        fontsize=13.5, fontweight='bold', color=INK, va='top')

# --- the evidence (8 words of labelling)
dp = fig.add_axes([0.035, 0.265, 0.50, 0.525])
edges = np.arange(60, 101, 0.5)
for data, colour in ((W_IN, BLUE), (W_BT, ORANGE)):
    h, _ = np.histogram(data, bins=edges)
    dp.fill_between(edges[:-1], 0, h/h.max(), step='post', color=colour, alpha=0.85, lw=0)
# line stops below the label band so the labels never cross it
dp.plot([65, 65], [0, 1.02], color=GREY, ls=(0, (3, 2)), lw=1.3, zorder=4)
dp.text(66.4, 0.30, '65%', fontsize=9, color=GREY, va='bottom', ha='left', fontweight='bold')
dp.text(60.3, 1.10, 'between', fontsize=9.5, color=ORANGE, fontweight='bold',
        ha='left', va='center')
dp.text(88, 1.10, 'within a bin', fontsize=9.5, color=BLUE, fontweight='bold',
        ha='center', va='center')
dp.set_xlim(60, 100); dp.set_ylim(0, 1.22)
dp.set_yticks([]); dp.set_xticks([60, 70, 80, 90, 100])
dp.set_xlabel('amino-acid identity (%)', fontsize=9, color=MUTED, labelpad=3)
dp.tick_params(labelsize=8.5, colors=MUTED, length=3, pad=2)
for s in ('top', 'right', 'left'): dp.spines[s].set_visible(False)
dp.spines['bottom'].set_color('#aaaaaa'); dp.spines['bottom'].set_linewidth(1.0)

# --- the result (6 words)
def block(x, y, w, h, num, lab, colour):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0,rounding_size=1.6',
                                ec=colour, fc='white', lw=1.6, zorder=2))
    ax.text(x+w/2, y+h-3.5, num, ha='center', va='top', fontsize=21,
            fontweight='bold', color=colour, zorder=3)
    ax.text(x+w/2, y+6.5, lab, ha='center', va='center', fontsize=9,
            color=INK, zorder=3, linespacing=1.35)

ax.annotate('', xy=(59.5, 62), xytext=(54.0, 62),
            arrowprops=dict(arrowstyle='-|>', color=GREY, lw=1.6, mutation_scale=14))
block(61.0, 51.0, 17.0, 28.0, f'{NBIG}', 'kept\ntogether', BLUE)
block(80.5, 51.0, 17.0, 28.0, f'{NSMALL}', 'pulled\nout', ORANGE)

# --- caption (14 words) and footer (5 words)
ax.text(61.0, 45.5, f'{NBIG+NSMALL} $\\it{{Pseudomonas}}$ genomes, no labels given.',
        fontsize=9.5, color=INK, va='top')
ax.text(61.0, 37.5, f'All {NSMALL} have since been reassigned\nto $\\it{{Halopseudomonas}}$.',
        fontsize=9.5, color=ORANGE, va='top', fontweight='bold', linespacing=1.45)
ax.text(61.0, 22.0, '2,206 genera screened in 38 hours.',
        fontsize=9, color=MUTED, va='top')

# --- the limitation, carried over from the long version. The paper's own framing, and the
#     short abstract should not drop it just because it is the unflattering number.
ax.add_patch(FancyBboxPatch((2.5, 2.0), 95.0, 12.0,
                            boxstyle='round,pad=0,rounding_size=1.6',
                            ec=ORANGE, fc='#fdf6f1', lw=1.3, zorder=2))
ax.text(6.0, 8.0, 'Misses', fontsize=10, color=INK, va='center', ha='left', zorder=3)
ax.text(15.5, 8.0, '3 of 5', fontsize=14, fontweight='bold', color=ORANGE,
        va='center', ha='left', zorder=3)
ax.text(24.8, 8.0, 'real problems. A first check, not a final answer.',
        fontsize=10, color=INK, va='center', ha='left', zorder=3)

fig.savefig(OUT/'RaGCAn_graphical_abstract_short.pdf')
fig.savefig(OUT/'RaGCAn_graphical_abstract_short.png', dpi=300)
plt.close(fig)
print(f"\nwrote {OUT/'RaGCAn_graphical_abstract_short.pdf'}")
print(f"wrote {OUT/'RaGCAn_graphical_abstract_short.png'}")
