"""
ABE OG pairwise LFC correlation grid — lower triangle, hexbin density.

Follows extended_figure11_CLEANED.ipynb cell 2 (the KO-screen panel) but on the
ABE OG screen with the sensor/input filters the user specified:
  - targeting guides only
  - sensor target base editing >= 20%
  - Input_median >= 100
"""
import os
import warnings

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.stats

warnings.filterwarnings('ignore')

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Helvetica', 'Arial', 'DejaVu Sans']

BASE = '/Users/kexindong/Documents/GitHub/B-ALL-in-vivo-base-editing/'
OUT = BASE + 'figures/abe-og-correlation/'
os.makedirs(OUT, exist_ok=True)

# ── Filters (as requested) ────────────────────────────────────────────────
MIN_EDIT = 20     # % sensor target base editing
MIN_INPUT = 100   # median input counts

# Tissue panels, in display order. 'd15' is the in vitro arm.
SAMPLES = ['d15', 'spleen', 'bm', 'men']
S_DICT = {'d15': 'LFC In Vitro', 'spleen': 'LFC Spleen',
          'bm': 'LFC Bone Marrow', 'men': 'LFC Meninges'}

df = pd.read_csv(BASE + 'LFC-FDR/ABE_OG_LFC_FDR.csv', low_memory=False)
n_total = len(df)
t = df[df['classification'] == 'targeting guide']
n_targ = len(t)

# Editing filter: a guide is kept if it reaches MIN_EDIT in ANY of the plotted
# tissues (same convention as the UMAP notebook), so every panel shares one
# guide set — a pairwise grid has to plot the same points on both axes.
# Note this means a point may appear in a panel whose tissue it did not edit
# well in; requiring ALL four instead would give n=662.
edit_cols = [f'target_base_edit_perc_{s}' for s in SAMPLES]
t = t[(t[edit_cols] >= MIN_EDIT).any(axis=1)]
n_edit = len(t)

t = t[t['Input_median'] >= MIN_INPUT]
n_input = len(t)

lfc_cols = [f'LFC_median_{s}' for s in SAMPLES]
t = t.dropna(subset=lfc_cols)
n_final = len(t)

print(f'total rows          : {n_total}')
print(f'targeting guides    : {n_targ}')
print(f'editing >={MIN_EDIT}% (any) : {n_edit}')
print(f'Input_median >={MIN_INPUT}   : {n_input}')
print(f'final (LFC complete): {n_final}')

# ── Plot: 4x4, keep lower triangle only ───────────────────────────────────
n = len(SAMPLES)
fig, ax = plt.subplots(n, n, figsize=(10, 10), sharex=True, sharey=True)

for row, samp in enumerate(SAMPLES):
    for col, samp2 in enumerate(SAMPLES):
        a = ax[row][col]
        x = t[f'LFC_median_{samp2}']
        y = t[f'LFC_median_{samp}']
        a.hexbin(x, y, bins='log', gridsize=30, linewidth=0,
                 extent=(-8, 8, -8, 8))
        a.plot([0, 0], [-7.5, 6], linestyle='dashed', color='black', linewidth=1)
        a.plot([-8, 5], [0, 0], linestyle='dashed', color='black', linewidth=1)
        a.plot([-8, 6], [-8, 6], linestyle='dashed', color='black', linewidth=1)
        ax[row][0].set_ylabel(S_DICT[samp], fontsize=13)
        ax[n - 1][col].set_xlabel(S_DICT[samp2], fontsize=13)
        a.spines[['top', 'right']].set_visible(False)
        a.tick_params(axis='both', which='major', labelsize=13)

        spear = np.round(scipy.stats.spearmanr(x, y).statistic, 2)
        pear = np.round(scipy.stats.pearsonr(x, y).statistic, 2)
        a.set_title('$R_{p} = $' + f'{pear}' + ' | $R_{s} = $' + f'{spear}',
                    fontsize=13)

# drop the diagonal and upper triangle
for row in range(n):
    for col in range(n):
        if col >= row:
            fig.delaxes(ax[row][col])

fig.tight_layout()
fig.savefig(OUT + 'ABE_OG_LFC_correlation_grid.pdf', bbox_inches='tight')
fig.savefig(OUT + 'ABE_OG_LFC_correlation_grid.png', dpi=200,
            bbox_inches='tight')
print('\nwrote', OUT + 'ABE_OG_LFC_correlation_grid.pdf')

# record the exact correlations for the legend / source data
rows = []
for row, samp in enumerate(SAMPLES):
    for col, samp2 in enumerate(SAMPLES):
        if col >= row:
            continue
        x = t[f'LFC_median_{samp2}']
        y = t[f'LFC_median_{samp}']
        rows.append({
            'x': S_DICT[samp2].replace('LFC ', ''),
            'y': S_DICT[samp].replace('LFC ', ''),
            'pearson_r': round(scipy.stats.pearsonr(x, y).statistic, 3),
            'spearman_rho': round(scipy.stats.spearmanr(x, y).statistic, 3),
            'n_gRNA': len(t)})
stats = pd.DataFrame(rows)
stats.to_csv(OUT + 'ABE_OG_LFC_correlation_stats.csv', index=False)
print(stats.to_string(index=False))

with open(OUT + 'filter_report.txt', 'w') as fh:
    fh.write(f'ABE OG LFC correlation grid — filtering\n'
             f'total rows in table      : {n_total}\n'
             f'targeting guides         : {n_targ}\n'
             f'sensor editing >={MIN_EDIT}% (any of 4 tissues): {n_edit}\n'
             f'Input_median >={MIN_INPUT}      : {n_input}\n'
             f'final (complete LFC)     : {n_final}\n')
plt.close(fig)


# ── Variant: extent scaled to the ABE data instead of the KO screen's -8..8 ──
lo = np.floor(t[lfc_cols].min().min())
hi = np.ceil(t[lfc_cols].max().max())
lim = max(abs(lo), abs(hi))
EXT = (-lim, lim, -lim, lim)
print(f'\ntight extent: +/-{lim:.0f}')

fig, ax = plt.subplots(n, n, figsize=(10, 10), sharex=True, sharey=True)
for row, samp in enumerate(SAMPLES):
    for col, samp2 in enumerate(SAMPLES):
        a = ax[row][col]
        x = t[f'LFC_median_{samp2}']
        y = t[f'LFC_median_{samp}']
        a.hexbin(x, y, bins='log', gridsize=30, linewidth=0, extent=EXT)
        a.axvline(0, linestyle='dashed', color='black', linewidth=1)
        a.axhline(0, linestyle='dashed', color='black', linewidth=1)
        a.plot([-lim, lim], [-lim, lim], linestyle='dashed', color='black',
               linewidth=1)
        a.set_xlim(-lim, lim); a.set_ylim(-lim, lim)
        ax[row][0].set_ylabel(S_DICT[samp], fontsize=13)
        ax[n - 1][col].set_xlabel(S_DICT[samp2], fontsize=13)
        a.spines[['top', 'right']].set_visible(False)
        a.tick_params(axis='both', which='major', labelsize=13)
        spear = np.round(scipy.stats.spearmanr(x, y).statistic, 2)
        pear = np.round(scipy.stats.pearsonr(x, y).statistic, 2)
        a.set_title('$R_{p} = $' + f'{pear}' + ' | $R_{s} = $' + f'{spear}',
                    fontsize=13)
for row in range(n):
    for col in range(n):
        if col >= row:
            fig.delaxes(ax[row][col])
fig.tight_layout()
fig.savefig(OUT + 'ABE_OG_LFC_correlation_grid_tightaxes.pdf',
            bbox_inches='tight')
fig.savefig(OUT + 'ABE_OG_LFC_correlation_grid_tightaxes.png', dpi=200,
            bbox_inches='tight')
plt.close(fig)
print('wrote tight-axes variant')
