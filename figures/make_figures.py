#!/usr/bin/env python3
"""Figures for the RaGCAn manuscript. Every value is read from the result tables, none typed in."""
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
ag = list(csv.DictReader(open(R/'LPSN_Genus_Survey'/'phase3'/'agreement.tsv'), delimiter='\t'))
pocp = {r['genus']: int(r['bins']) for r in
        csv.DictReader(open(R/'POCP_run'/'reports'/'pocp_status.tsv'), delimiter='\t')
        if r['state'] == 'done' and r['bins']}
for r in ag:
    r['b'] = int(r['bins']); r['g'] = int(r['gtdb_genera'])
    r['a'] = float(r['ari']) if r['ari'] not in ('nan', '') else None

plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False,
                     'figure.dpi': 300, 'savefig.bbox': 'tight'})
GREY, BLUE, ORANGE = '#4d4d4d', '#2c6fa8', '#c2703a'

# ---- FIG 1: the four outcomes -------------------------------------------------
both_c = sum(1 for r in ag if r['b'] == 1 and r['g'] == 1)
both_s = sum(1 for r in ag if r['b'] > 1 and r['g'] > 1)
lump   = sum(1 for r in ag if r['b'] == 1 and r['g'] > 1)
split  = sum(1 for r in ag if r['b'] > 1 and r['g'] == 1)
fig, ax = plt.subplots(figsize=(6.2, 2.6))
labels = ['Both agree:\ncoherent', 'Both agree:\ndivided',
          'RaGCAn one bin,\nGTDB divides', 'RaGCAn divides,\nGTDB one genus']
vals = [both_c, both_s, lump, split]
cols = [GREY, BLUE, ORANGE, ORANGE]
b = ax.barh(range(4), vals, color=cols, height=0.62)
for i, v in enumerate(vals):
    ax.text(v + 12, i, f'{v}  ({v/len(ag)*100:.1f}%)', va='center', fontsize=8.5)
ax.set_yticks(range(4)); ax.set_yticklabels(labels, fontsize=8.5)
ax.invert_yaxis(); ax.set_xlim(0, max(vals)*1.22)
ax.set_xlabel(f'Genera (n = {len(ag)})')
ax.axvline(0, color='k', lw=0.8)
fig.savefig(OUT/'Figure1_outcomes.pdf'); fig.savefig(OUT/'Figure1_outcomes.png'); plt.close(fig)

# ---- FIG 2: ARI distribution where both divide --------------------------------
a = sorted(r['a'] for r in ag if r['b'] > 1 and r['g'] > 1 and r['a'] is not None)
fig, ax = plt.subplots(figsize=(6.2, 2.9))
ax.hist(a, bins=np.arange(-0.05, 1.06, 0.05), color=BLUE, edgecolor='white', lw=0.6)
ax.axvline(0, color=GREY, ls=':', lw=1.2)
ax.axvline(float(np.median(a)), color=ORANGE, ls='--', lw=1.4)
ax.text(0.01, ax.get_ylim()[1]*0.92, 'chance', color=GREY, fontsize=8, rotation=90, va='top')
ax.text(np.median(a)+0.012, ax.get_ylim()[1]*0.97, f'median {np.median(a):.2f}',
        color=ORANGE, fontsize=8, va='top')
n1 = sum(1 for x in a if x >= 0.999)
ax.annotate(f'{n1} genera partitioned\nidentically to GTDB', xy=(0.99, n1*0.92), xytext=(0.40, n1*0.55),
            fontsize=8.5, ha='center', arrowprops=dict(arrowstyle='->', lw=0.9, color='k'))
ax.set_xlabel('Adjusted Rand index vs GTDB'); ax.set_ylabel('Genera')
ax.set_title(f'Agreement where both divide the genus (n = {len(a)})', fontsize=9.5, loc='left')
fig.savefig(OUT/'Figure2_ari.pdf'); fig.savefig(OUT/'Figure2_ari.png'); plt.close(fig)

# ---- FIG 3: RaGCAn vs POCP -------------------------------------------------
# Horizontal layout: labels sit beside their bars, values are printed, and the overlap
# between each pair is visible directly rather than needing a caption to point it out.
def wilson(k, n, z=1.959963985):
    p = k/n; d = 1+z*z/n; c = (p+z*z/(2*n))/d
    h = z*math.sqrt(p*(1-p)/n + z*z/(4*n*n))/d
    return p*100, max(0, c-h)*100, min(1, c+h)*100
common = [r for r in ag if r['genus'] in pocp]
def stats(get):
    TP = sum(1 for r in common if get(r) > 1 and r['g'] > 1)
    FP = sum(1 for r in common if get(r) > 1 and r['g'] == 1)
    FN = sum(1 for r in common if get(r) == 1 and r['g'] > 1)
    TN = sum(1 for r in common if get(r) == 1 and r['g'] == 1)
    return [wilson(TP, TP+FP), wilson(TN, TN+FP), wilson(TP, TP+FN)]
S = {'RaGCAn': stats(lambda r: r['b']), 'POCP': stats(lambda r: pocp[r['genus']])}

names = ['Precision\nwhen it flags a genus,\nGTDB agrees',
         'Specificity\nsound genera\nleft alone',
         'Recall\nof genera needing\nattention, found']
fig, axes = plt.subplots(3, 1, figsize=(6.4, 4.4), sharex=True)
for k, ax in enumerate(axes):
    for row, (lab, col) in enumerate((('RaGCAn', BLUE), ('POCP', ORANGE))):
        m, lo, hi = S[lab][k]
        y = 0.72 - row*0.44
        ax.plot([lo, hi], [y, y], color=col, lw=2.4, solid_capstyle='butt', zorder=2)
        for x in (lo, hi):
            ax.plot([x, x], [y-0.10, y+0.10], color=col, lw=1.6, zorder=2)
        ax.plot(m, y, 'o', color=col, ms=7, zorder=3,
                markeredgecolor='white', markeredgewidth=0.9)
        # name and value both on the RIGHT, so nothing can collide with the row label
        ax.text(hi + 1.6, y, f'{lab}  {m:.1f}%', va='center', ha='left',
                fontsize=8.5, color=col, fontweight='bold')
    # shade the region where the two intervals overlap
    a, b = S['RaGCAn'][k], S['POCP'][k]
    ov_lo, ov_hi = max(a[1], b[1]), min(a[2], b[2])
    if ov_hi > ov_lo:
        ax.axvspan(ov_lo, ov_hi, color='#999999', alpha=0.16, zorder=0, lw=0)
    ax.set_ylim(0.05, 1.05); ax.set_yticks([])
    ax.set_ylabel(names[k], rotation=0, ha='right', va='center', fontsize=8.5, labelpad=8)
    for sp in ('left', 'right', 'top'):
        ax.spines[sp].set_visible(False)
    ax.tick_params(axis='y', length=0)
axes[-1].set_xlim(28, 118)
axes[-1].set_xlabel('% (points are estimates, bars are Wilson 95% confidence intervals)', fontsize=8.5)
axes[0].set_title(f'RaGCAn and POCP scored against GTDB on the same {len(common)} genera',
                  fontsize=9.5, loc='left', pad=10)
fig.text(0.5, -0.035, 'Grey shading marks where the two intervals overlap. '
         'They overlap on all three measures:\nthe difference between the methods is not '
         'statistically distinguishable on these data.',
         ha='center', fontsize=8.2, style='italic', color=GREY)
fig.subplots_adjust(hspace=0.28)
fig.savefig(OUT/'Figure3_pocp.pdf'); fig.savefig(OUT/'Figure3_pocp.png'); plt.close(fig)

print(f"  Fig 3: n={len(common)}  RaGCAn {S['RaGCAn'][0][0]:.1f}%  POCP {S['POCP'][0][0]:.1f}%")
