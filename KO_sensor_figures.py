"""
KO (Cas9-NG) sensor indel analysis -- figure panels.

Panels
------
a  Sensor indel efficiency by tissue (bar + replicate dots)
b  Distribution of target_indel_perc, split by COSMIC role
c  KO vs ABE original LFC correlation, paneled by COSMIC role, per tissue
d  KO sensor indel % vs KO LFC, by COSMIC role

Data
----
KO_sensor_indel_table.csv   this analysis (CRISPResso indel quantification)
KO_LFC_FDR.csv              published KO LFC/FDR            (SG repo)
ABE_LFC_FDR_df.csv          published ABE original LFC/FDR  (SG repo)
Census_all*.csv             COSMIC Cancer Gene Census       (SG repo)
"""

import os
import warnings

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats
import seaborn as sns
from matplotlib.patches import Patch

warnings.filterwarnings('ignore')

# ── Global style: final print size, Nature look ───────────────────────────
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['text.usetex'] = False
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Helvetica']
mpl.rcParams.update({'mathtext.fontset': 'custom', 'mathtext.rm': 'Helvetica',
                     'mathtext.it': 'Arial:italic', 'mathtext.bf': 'Arial:bold'})
mpl.rcParams.update({'font.size': 7, 'axes.labelsize': 7, 'axes.titlesize': 7,
                     'xtick.labelsize': 6, 'ytick.labelsize': 6, 'legend.fontsize': 6,
                     'axes.linewidth': .6, 'lines.linewidth': .75, 'patch.linewidth': .5,
                     'xtick.major.width': .6, 'ytick.major.width': .6,
                     'xtick.major.size': 3, 'ytick.major.size': 3,
                     'errorbar.capsize': 2, 'legend.frameon': False, 'figure.dpi': 200})
MM = 1 / 25.4

blues = sns.color_palette('Blues').as_hex()
greens = sns.color_palette('Greens').as_hex()
reds = sns.color_palette('Reds').as_hex()
yellows = sns.color_palette('Oranges').as_hex()
purples = sns.color_palette('Purples').as_hex()
greys = sns.color_palette('Greys').as_hex()

BASE = '/Users/kexindong/Documents/GitHub/B-ALL-in-vivo-base-editing/'
SG_BASE = os.environ.get('SG_BALL_ROOT',
                         '/Users/kexindong/Documents/GitHub/SG-B-ALL-in-vivo-base-editing/')
OUT = BASE + 'figures/KO_sensor/'
os.makedirs(OUT, exist_ok=True)

# fixed semantic colors (sci-plot conventions)
COSMIC_COLORS = {'Oncogene': 'tab:red', 'Oncogene/TSG': 'tab:purple',
                 'TSG': 'tab:blue', 'Undefined': 'tab:grey'}
COSMIC_ORDER = ['Oncogene', 'Oncogene/TSG', 'TSG', 'Undefined']
TISSUE_COLORS = {'d15': purples[5], 'spleen': blues[3],
                 'bonemarrow': greens[3], 'meninges': reds[3]}
TISSUE_LABEL = {'d15': 'In vitro', 'spleen': 'Spleen',
                'bonemarrow': 'Bone marrow', 'meninges': 'Meninges'}
TISSUES = ['d15', 'spleen', 'bonemarrow', 'meninges']

# replicate columns per tissue (see KO_sample_tissue_map.csv)
TISSUE_REPS = {'d15': ['d15_rep1', 'd15_rep2', 'd15_rep3'],
               'spleen': ['S1', 'S2', 'S3_5'],
               'bonemarrow': ['B1', 'B2', 'B3', 'B4', 'B5'],
               'meninges': ['M1', 'M2', 'M3_5']}

MIN_EDIT = 20.0     # indel % cutoff, mirrors the ABE arm's target_base_edit_perc >= 20
MIN_READS = 100     # Input_median_KO cutoff -- panel d only (as in
                    # extended_figure11 cell 2, the KO-only tissue grid).
                    # Panel c deliberately omits it: cell 3 does not use it.


def italic(g):
    return f"$\\mathit{{{g}}}$"


def corr_text(x, y):
    """Both correlations, as the sci-plot convention requires."""
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 3:
        return 'n/a'
    r = scipy.stats.pearsonr(x[ok], y[ok]).statistic
    rho = scipy.stats.spearmanr(x[ok], y[ok]).statistic
    return f"$\\rho$ = {rho:.2f} and R = {r:.2f}"


# ══════════════════════════════════════════════════════════════════════════
# Load and annotate
# ══════════════════════════════════════════════════════════════════════════
def load():
    ind = pd.read_csv(BASE + 'KO_sensor_indel_table.csv')
    ko = pd.read_csv(SG_BASE + 'KO_LFC_FDR.csv')
    abe = pd.read_csv(SG_BASE + 'ABE_LFC_FDR_df.csv')

    # ── COSMIC role, joined on the HUMAN symbol (census is human-only) ────
    cosmic = pd.read_csv(
        SG_BASE + 'source_data/Census_allSun Nov 17 02_26_47 2024.csv').fillna('Undefined')
    t_dict = {'TSG': 'TSG', 'TSG, fusion': 'TSG', 'Undefined': 'Undefined',
              'fusion': 'Undefined', 'oncogene': 'Oncogene',
              'oncogene, TSG': 'Oncogene/TSG',
              'oncogene, TSG, fusion': 'Oncogene/TSG',
              'oncogene, fusion': 'Oncogene'}
    role_by_symbol = {}
    for _, r in cosmic.iterrows():
        role_by_symbol.setdefault(r['Gene Symbol'], t_dict[r['Role in Cancer']])
    ko['Role in Cancer'] = [role_by_symbol.get(g, 'Undefined') for g in ko['gene_name_h']]

    # KO LFC columns renamed so they do not collide with the ABE ones
    ren = {'sgrna': 'gRNA_id'}
    for t in ['d15', 'spleen', 'bonemarrow', 'meninges']:
        ren[f'LFC_median_{t}'] = f'LFC_median_{t}_KO'
    ko = ko.rename(columns=ren)

    # Both tables carry Input_median and they are NOT the same (KO >=100:
    # 530 guides, ABE >=100: 779). Name them explicitly so no panel can
    # silently filter on the wrong arm.
    ko = ko.rename(columns={'Input_median': 'Input_median_KO'})

    m = ko.merge(ind, on='gRNA_id', how='left')
    m = m.merge(abe[['gRNA_id', 'classification', 'target_base_edit_perc',
                     'Input_median']
                    + [f'LFC_median_{t}' for t in TISSUES]]
                .rename(columns={'Input_median': 'Input_median_ABE'}),
                on='gRNA_id', how='left', suffixes=('', '_ABE'))
    return m


# ══════════════════════════════════════════════════════════════════════════
# a  Sensor indel efficiency by tissue
# ══════════════════════════════════════════════════════════════════════════
def panel_a(m):
    fig, ax = plt.subplots(figsize=(5, 3.2))

    tg = m[m['classification'] == 'targeting guide']
    xs = np.arange(len(TISSUES))

    for i, t in enumerate(TISSUES):
        reps = [f'indel_perc_{r}' for r in TISSUE_REPS[t]
                if f'indel_perc_{r}' in tg.columns]
        rep_medians = [tg[c].median() for c in reps]
        mean = np.mean(rep_medians)
        sd = np.std(rep_medians, ddof=1) if len(rep_medians) > 1 else 0.0

        ax.bar(i, mean, color=TISSUE_COLORS[t], edgecolor='black',
               linewidth=1, zorder=1, width=.65)
        ax.errorbar(i, mean, yerr=sd, color='black', lw=1, capsize=3, zorder=2)
        jitter = np.linspace(-.13, .13, len(rep_medians))
        ax.scatter(i + jitter, rep_medians, color='black', s=28,
                   linewidth=.5, edgecolors='black', zorder=3, alpha=.8)

    ax.set_xticks(xs)
    ax.set_xticklabels([TISSUE_LABEL[t] for t in TISSUES])
    ax.set_ylabel('Median sensor indel (%)', fontsize=13)
    ax.set_ylim(0, 75)
    ax.set_yticks([0, 20, 40, 60])
    ax.spines[['top', 'right']].set_visible(False)
    ax.tick_params(axis='both', which='major', labelsize=13)
    ax.set_title('Cas9-NG sensor editing', fontsize=13)

    fig.tight_layout()
    fig.savefig(OUT + 'KO_a_indel_by_tissue.pdf', bbox_inches='tight')
    fig.savefig(OUT + 'KO_a_indel_by_tissue.png', dpi=150, bbox_inches='tight')
    plt.close(fig)


# ══════════════════════════════════════════════════════════════════════════
# b  Indel distribution split by COSMIC role
# ══════════════════════════════════════════════════════════════════════════
def panel_b(m):
    fig, ax = plt.subplots(figsize=(6, 3.4))

    tg = m[(m['classification'] == 'targeting guide')
           & m['target_indel_perc'].notna()]

    data = [tg[tg['Role in Cancer'] == r]['target_indel_perc'].values
            for r in COSMIC_ORDER]

    bp = ax.boxplot(data, positions=np.arange(len(COSMIC_ORDER)),
                    widths=.6, showfliers=False, patch_artist=True,
                    medianprops=dict(color='black', lw=1),
                    whiskerprops=dict(lw=.6), capprops=dict(lw=.6),
                    boxprops=dict(lw=.6))
    for patch, r in zip(bp['boxes'], COSMIC_ORDER):
        patch.set_facecolor(COSMIC_COLORS[r])
        patch.set_alpha(.55)
        patch.set_edgecolor('black')

    rng = np.random.default_rng(0)
    for i, (r, d) in enumerate(zip(COSMIC_ORDER, data)):
        if len(d) == 0:
            continue
        x = i + rng.uniform(-.18, .18, len(d))
        ax.scatter(x, d, s=6, color=COSMIC_COLORS[r], alpha=.35,
                   linewidths=0, zorder=3)

    ax.axhline(MIN_EDIT, color='grey', ls='--', lw=.8, zorder=0)
    ax.set_xticks(np.arange(len(COSMIC_ORDER)))
    ax.set_xticklabels([f"{r}\n(n={len(d)})" for r, d in zip(COSMIC_ORDER, data)])
    ax.set_ylabel('Sensor indel (%)', fontsize=13)
    ax.set_ylim(-3, 105)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.spines[['top', 'right']].set_visible(False)
    ax.tick_params(axis='both', which='major', labelsize=13)
    ax.set_title('Sensor editing by COSMIC role', fontsize=13)

    fig.tight_layout()
    fig.savefig(OUT + 'KO_b_indel_by_role.pdf', bbox_inches='tight')
    fig.savefig(OUT + 'KO_b_indel_by_role.png', dpi=150, bbox_inches='tight')
    plt.close(fig)


# ══════════════════════════════════════════════════════════════════════════
# c  KO vs ABE original LFC, paneled by COSMIC role (rows) x tissue (cols)
# ══════════════════════════════════════════════════════════════════════════
def panel_c(m, filter_edit=True):
    """
    Matches the published Figure 3b / extended_figure11 style exactly:
    log-scaled hexbin density, equal -8..8 axes with the y = x diagonal,
    '$R_p$ | $R_s$' titles, screen-mode font sizes. Rows add the COSMIC
    role split; the top row reproduces the published panel (All guides).
    """
    rows = ['All'] + COSMIC_ORDER
    names = ['In Vitro', 'Spleen', 'Bone', 'Meninges']

    # draw-big screen mode, to match the existing figure
    fig, axs = plt.subplots(len(rows), len(TISSUES), figsize=(10, 2.9 * len(rows)),
                            sharex=True, sharey=True)

    # Published baseline (extended_figure11 cell 3): targeting guides with
    # ABE target_base_edit_perc >= 20, merged with the KO table. NO
    # Input_median cutoff -- that belongs to cell 2 (the KO-only tissue
    # grid), not to this panel. Applying it here gave n = 711 instead of
    # 795 and shifted every correlation.
    sub = m[(m['classification'] == 'targeting guide')
            & (m['target_base_edit_perc'] >= MIN_EDIT)]
    if filter_edit:
        # the new, symmetric filter: the KO arm is held to the same
        # >=20% editing bar as the ABE arm
        sub = sub[sub['target_indel_perc'] >= MIN_EDIT]

    for i, role in enumerate(rows):
        d_role = sub if role == 'All' else sub[sub['Role in Cancer'] == role]
        for j, t in enumerate(TISSUES):
            ax = axs[i][j]
            x = d_role[f'LFC_median_{t}'].to_numpy(dtype=float)       # ABE
            y = d_role[f'LFC_median_{t}_KO'].to_numpy(dtype=float)    # KO
            ok = np.isfinite(x) & np.isfinite(y)
            x, y = x[ok], y[ok]

            ax.hexbin(x, y, bins='log', gridsize=30, linewidth=0,
                      extent=(-8, 8, -8, 8))
            ax.plot([0, 0], [-8, 8], linestyle='dashed', color='black', linewidth=1)
            ax.plot([-8, 8], [0, 0], linestyle='dashed', color='black', linewidth=1)
            ax.plot([-8, 8], [-8, 8], linestyle='dashed', color='black', linewidth=1)
            ax.spines[['top', 'right']].set_visible(False)
            ax.tick_params(axis='both', which='major', labelsize=13)
            ax.set_xticks([-5, 0, 5])
            ax.set_yticks([-5, 0, 5])

            spear = np.round(scipy.stats.spearmanr(x, y).statistic, 2)
            pear = np.round(scipy.stats.pearsonr(x, y).statistic, 2)
            head = names[j] if i == 0 else ''
            ax.set_title((f'{head}\n' if head else '')
                         + '$R_{p} = $' + f'{pear}'
                         + ' | $R_{s} = $' + f'{spear}'
                         + f'\nn = {len(x)}', fontsize=13)

            if j == 0:
                lbl = 'All guides' if role == 'All' else role
                ax.set_ylabel(f'{lbl}\nLFC CRISPR KO', fontsize=13)
            if i == len(rows) - 1:
                ax.set_xlabel('LFC ABE', fontsize=13)

    tag = 'edit_filtered' if filter_edit else 'unfiltered'
    fig.tight_layout()
    fig.savefig(OUT + f'KO_c_LFC_corr_by_role_{tag}.pdf',
                bbox_inches='tight', transparent=True)
    fig.savefig(OUT + f'KO_c_LFC_corr_by_role_{tag}.png',
                dpi=150, bbox_inches='tight')
    plt.close(fig)


# ══════════════════════════════════════════════════════════════════════════
# d  Sensor indel % vs KO LFC, by role
# ══════════════════════════════════════════════════════════════════════════
def panel_d(m):
    names = ['In Vitro', 'Spleen', 'Bone', 'Meninges']
    fig, axs = plt.subplots(1, len(TISSUES), figsize=(10, 3),
                            sharex=True, sharey=True)

    sub = m[(m['classification'] == 'targeting guide')
            & (m['Input_median_KO'] >= MIN_READS)]

    for j, t in enumerate(TISSUES):
        ax = axs[j]
        x = sub[f'indel_perc_{t}'].to_numpy(dtype=float)
        y = sub[f'LFC_median_{t}_KO'].to_numpy(dtype=float)
        ok = np.isfinite(x) & np.isfinite(y)

        ax.hexbin(x[ok], y[ok], bins='log', gridsize=30, linewidth=0,
                  extent=(0, 100, -8, 8))
        ax.plot([0, 100], [0, 0], linestyle='dashed', color='black', linewidth=1)
        ax.plot([MIN_EDIT, MIN_EDIT], [-8, 8], linestyle='dashed',
                color='black', linewidth=1)

        spear = np.round(scipy.stats.spearmanr(x[ok], y[ok]).statistic, 2)
        pear = np.round(scipy.stats.pearsonr(x[ok], y[ok]).statistic, 2)
        ax.set_title(f'{names[j]}\n' + '$R_{p} = $' + f'{pear}'
                     + ' | $R_{s} = $' + f'{spear}', fontsize=13)

        ax.set_xlim(0, 100)
        ax.set_ylim(-8, 8)
        ax.set_xticks([0, 50, 100])
        ax.set_yticks([-5, 0, 5])
        ax.spines[['top', 'right']].set_visible(False)
        ax.tick_params(axis='both', which='major', labelsize=13)
        ax.set_xlabel('Sensor indel (%)', fontsize=13)
        if j == 0:
            ax.set_ylabel('LFC CRISPR KO', fontsize=13)

    fig.tight_layout()
    fig.savefig(OUT + 'KO_d_indel_vs_LFC.pdf', bbox_inches='tight')
    fig.savefig(OUT + 'KO_d_indel_vs_LFC.png', dpi=150, bbox_inches='tight')
    plt.close(fig)


if __name__ == '__main__':
    m = load()
    print("merged rows:", len(m))
    print("\nRole in Cancer (targeting guides):")
    print(m[m['classification'] == 'targeting guide']['Role in Cancer']
          .value_counts().to_string())
    print("\nguides with sensor indel:", int(m['target_indel_perc'].notna().sum()))

    panel_a(m)
    panel_b(m)
    panel_c(m, filter_edit=True)
    panel_c(m, filter_edit=False)
    panel_d(m)
    print("\nwrote panels to", OUT)
