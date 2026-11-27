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



# ══════════════════════════════════════════════════════════════════════════
# a2  Sensor indel efficiency, lollipop style matching the ABE editing panel
#     (horizontal bars, one open dot per replicate at the bar tip)
#     Green = Cas9-NG indel, as the counterpart of the purple/blue ABE panel.
# ══════════════════════════════════════════════════════════════════════════
#     The Input samples (278012/13/14) are the low-depth ones flagged
#     usable=NO in KO_sample_tissue_map.csv -- 7-9% barcode match, and two of
#     the three give a per-guide median of 0. Plotting them would read as
#     "no editing in Input" when it is really "not enough reads to call it",
#     so the row is left out. The ABE panel's Plasmid row has no KO
#     counterpart in this run either.
LOLLI_ROWS = [('d5', 'D5', ['d5_rep1', 'd5_rep2', 'd5_rep3']),
              ('d15', 'D15', ['d15_rep1', 'd15_rep2', 'd15_rep3']),
              ('spleen', 'Spleen', TISSUE_REPS['spleen']),
              ('bonemarrow', 'Bone Marrow', TISSUE_REPS['bonemarrow']),
              ('meninges', 'Meninges', TISSUE_REPS['meninges'])]


def panel_a2(m):
    """Horizontal lollipop bars in green, styled like the ABE 'Average
    Editing %' panel: flat bar to the mean, one black-edged dot per
    replicate, open axes, category names on the left."""
    fig, ax = plt.subplots(figsize=(5.4, 3.6))

    tg = m[m['classification'] == 'targeting guide']
    bar_color = greens[3]
    ys = np.arange(len(LOLLI_ROWS))[::-1]    # first row at the top

    for y, (_key, _label, reps) in zip(ys, LOLLI_ROWS):
        cols = [f'indel_perc_{r}' for r in reps if f'indel_perc_{r}' in tg.columns]
        rep_medians = [tg[c].median() for c in cols]
        mean = np.mean(rep_medians)

        ax.barh(y, mean, color=bar_color, edgecolor='black',
                linewidth=1, height=.55, zorder=1)
        # light fill + dark rim so a dot stays visible both on the bar and
        # off its end -- a same-green dot disappears inside the bar.
        ax.scatter(rep_medians, [y] * len(rep_medians), s=55,
                   facecolors=greens[1], edgecolors='black',
                   linewidth=1, zorder=3)

    ax.set_yticks(ys)
    ax.set_yticklabels([lab for _k, lab, _r in LOLLI_ROWS])
    ax.set_xlabel('Average Indel %', fontsize=13)
    ax.set_xlim(0, 70)
    ax.set_xticks([0, 10, 20, 30, 40, 50, 60, 70])
    ax.set_ylim(-.6, len(LOLLI_ROWS) - .4)
    ax.spines[['top', 'right']].set_visible(False)
    ax.tick_params(axis='both', which='major', labelsize=13)

    legend = [Patch(facecolor=bar_color, edgecolor='black',
                    label='Target Indel (Cas9-NG)')]
    ax.legend(handles=legend, loc='upper center', bbox_to_anchor=(.5, -.26),
              fontsize=11, frameon=False, handlelength=1.2, handleheight=1.2)

    fig.tight_layout()
    fig.savefig(OUT + 'KO_a2_indel_lollipop_green.pdf', bbox_inches='tight')
    fig.savefig(OUT + 'KO_a2_indel_lollipop_green.png', dpi=200,
                bbox_inches='tight')
    plt.close(fig)


# ══════════════════════════════════════════════════════════════════════════
# e  UpSet: how ABE and KO hit calls overlap across tissues
# f  Heatmaps of the concordant sets (UP in both / DOWN in both)
#
# Requested by a co-author. Read the caveat in `concordant_sets()` before
# interpreting: at the standard thresholds the "UP in both" set is EMPTY,
# and that is a real result about the screens, not a plotting bug.
# ══════════════════════════════════════════════════════════════════════════
HIT_TISSUES = ['spleen', 'bonemarrow', 'meninges']
LFC_CUTOFF = 1.0
FDR_CUTOFF = 0.1


def hit_frame(min_reads=MIN_READS, tissues=None):
    """ABE x KO merged on gRNA_id, with one boolean column per
    (arm, tissue, direction) hit call.

    Filters, and why each one is here:
      * targeting guides only, ABE target_base_edit_perc >= 20
        -- the ABE arm must actually have made its edit, as in the
        published Fig 3b merge.
      * Input_median >= min_reads on BOTH arms. This one matters a lot:
        every KO guide that scores as 'enriched' has an input median of
        11-32 reads, so its large positive LFC is a small-denominator
        artifact. Dropping this cutoff does not reveal biology, it
        reveals noise -- see concordant_sets().
    """
    tissues = HIT_TISSUES if tissues is None else list(tissues)
    abe = pd.read_csv(SG_BASE + 'ABE_LFC_FDR_df.csv')
    ko = pd.read_csv(SG_BASE + 'KO_LFC_FDR.csv')

    abe = abe[(abe['classification'] == 'targeting guide')
              & (abe['target_base_edit_perc'] >= MIN_EDIT)]
    m = abe.merge(ko, on='gRNA_id', suffixes=('_ABE', '_KO'))
    m = m[(m['Input_median_ABE'] >= min_reads)
          & (m['Input_median_KO'] >= min_reads)]

    for arm in ('ABE', 'KO'):
        for t in tissues:
            lfc = m[f'LFC_median_{t}_{arm}']
            m[f'{arm}_{t}_up'] = (lfc >= LFC_CUTOFF) & \
                                 (m[f'FDR_high_{t}_{arm}'] < FDR_CUTOFF)
            m[f'{arm}_{t}_dn'] = (lfc <= -LFC_CUTOFF) & \
                                 (m[f'FDR_low_{t}_{arm}'] < FDR_CUTOFF)
    # "a hit anywhere in vivo", the summary used by the UpSet and heatmaps
    for arm in ('ABE', 'KO'):
        for d in ('up', 'dn'):
            m[f'{arm}_{d}'] = m[[f'{arm}_{t}_{d}' for t in tissues]].any(axis=1)
    m.attrs['tissues'] = list(tissues)
    return m


def concordant_sets(m):
    """The four quadrants the co-author asked about, plus their sizes."""
    return {'UP in both':   m[m['ABE_up'] & m['KO_up']],
            'DOWN in both': m[m['ABE_dn'] & m['KO_dn']],
            'ABE up only':  m[m['ABE_up'] & ~m['KO_up']],
            'KO down only': m[m['KO_dn'] & ~m['ABE_dn']]}


def panel_e_upset(m, tag=''):
    """Hand-drawn UpSet over the six hit sets (ABE/KO x up/down pooled over
    tissues, plus the per-tissue breakdown). Drawn by hand rather than with
    upsetplot so it follows the same open-axes / fixed-color conventions as
    every other panel in this figure."""
    sets = {'ABE up': m['ABE_up'], 'ABE down': m['ABE_dn'],
            'KO up': m['KO_up'], 'KO down': m['KO_dn']}
    names = list(sets)
    mat = np.column_stack([sets[n].values for n in names])

    # every observed non-empty combination, largest first
    combos = {}
    for row in mat:
        key = tuple(row)
        if not any(key):
            continue
        combos[key] = combos.get(key, 0) + 1
    combos = sorted(combos.items(), key=lambda kv: -kv[1])

    n_c = len(combos)
    fig = plt.figure(figsize=(max(5.0, .62 * n_c + 2.4), 4.6))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.15, 3.2], height_ratios=[2.1, 1.3],
                          wspace=.05, hspace=.08)
    ax_bar = fig.add_subplot(gs[0, 1])
    ax_dot = fig.add_subplot(gs[1, 1], sharex=ax_bar)
    ax_set = fig.add_subplot(gs[1, 0], sharey=ax_dot)

    xs = np.arange(n_c)
    counts = [c for _k, c in combos]

    # ── intersection sizes ────────────────────────────────────────────────
    ax_bar.bar(xs, counts, color=greys[4], edgecolor='black',
               linewidth=.8, width=.62)
    for x, c in zip(xs, counts):
        ax_bar.text(x, c + max(counts) * .025, str(c), ha='center',
                    va='bottom', fontsize=9)
    ax_bar.set_ylabel('Guides in\nintersection', fontsize=11)
    ax_bar.set_ylim(0, max(counts) * 1.18)
    ax_bar.spines[['top', 'right']].set_visible(False)
    ax_bar.tick_params(axis='x', labelbottom=False, bottom=False)
    ax_bar.tick_params(axis='y', labelsize=10)

    # ── membership dots ───────────────────────────────────────────────────
    ys = np.arange(len(names))[::-1]
    for xi, (key, _c) in enumerate(combos):
        on = [y for y, k in zip(ys, key) if k]
        ax_dot.scatter([xi] * len(names), ys, s=60,
                       color=greys[1], zorder=1)
        if on:
            ax_dot.scatter([xi] * len(on), on, s=60, color='black', zorder=2)
            if len(on) > 1:
                ax_dot.plot([xi, xi], [min(on), max(on)], color='black',
                            lw=1.4, zorder=2)
    ax_dot.set_yticks(ys)
    ax_dot.set_yticklabels([''] * len(names))
    ax_dot.set_ylim(-.6, len(names) - .4)
    ax_dot.set_xlim(-.7, n_c - .3)
    ax_dot.set_xticks([])
    for s in ('top', 'right', 'bottom', 'left'):
        ax_dot.spines[s].set_visible(False)
    ax_dot.tick_params(left=False)

    # ── set sizes on the left ─────────────────────────────────────────────
    # Set names + total size as text only. A barh here shares the y axis with
    # the dot matrix, which overrides any x-limit set to reserve label space,
    # so the bars ran under their own labels. The totals are in the text.
    tot = [int(sets[n].sum()) for n in names]
    for y, v, n in zip(ys, tot, names):
        ax_set.text(1.0, y, f'{n}  ({v})', ha='right', va='center',
                    fontsize=10, transform=ax_set.get_yaxis_transform())
    ax_set.set_xlim(0, 1)
    ax_set.axis('off')

    tis = m.attrs.get('tissues', HIT_TISSUES)
    where = ('' if len(tis) > 1 else f' in {TISSUE_LABEL[tis[0]].lower()}')
    ax_bar.set_title(
        f'ABE vs Cas9-NG hit overlap{where}  |  {len(m)} guides, '
        f'|LFC| $\\geq$ {LFC_CUTOFF:g}, FDR < {FDR_CUTOFF:g}', fontsize=11)

    fig.savefig(OUT + f'KO_e_upset_hits{tag}.pdf', bbox_inches='tight')
    fig.savefig(OUT + f'KO_e_upset_hits{tag}.png', dpi=200,
                bbox_inches='tight')
    plt.close(fig)
    return combos


def _cosmic_labels(sub, cosmic_map):
    """Two-line label per guide, with the COSMIC role appended:
       'TP53 p.R248Q [TSG]'
       'Trp53 p.R245Q'
    Colored by role when the axis is drawn."""
    labels, roles = [], []
    for _i, v in sub.iterrows():
        gh = v['gene_name_h_ABE']
        ph = v['HGVSp_h_ABE'] if isinstance(v['HGVSp_h_ABE'], str) else 'Non-coding'
        gm = v['gene_name_m_corrected_ABE']
        pm = v['HGVSp_m_ABE'] if isinstance(v['HGVSp_m_ABE'], str) else ''
        role = cosmic_map.get(gh, 'Undefined')
        labels.append(f'{gh} {ph} [{role}]\n{gm} {pm}'.rstrip())
        roles.append(role)
    return labels, roles


def _heatmap(sub, title, fname, cosmic_map, sort_col, tissues=None):
    if len(sub) == 0:
        print(f'  [skip] {title}: no guides in this set')
        return
    tissues = HIT_TISSUES if tissues is None else list(tissues)
    sub = sub.sort_values(by=sort_col, ascending=False)
    cols = [f'LFC_median_{t}_{a}' for a in ('ABE', 'KO') for t in tissues]
    dat = sub[cols]

    fig, ax = plt.subplots(figsize=(0.62 * len(cols) + 5.0,
                                    max(2.2, .52 * len(sub) + 1.3)))
    g = sns.heatmap(dat, cmap='seismic', vmin=-5, vmax=5, annot=True, fmt='.1f',
                    annot_kws={'size': 9}, linewidth=1, linecolor='white',
                    cbar_kws={'label': 'median LFC', 'shrink': .45,
                              'aspect': 12, 'pad': .06}, ax=ax)

    labels, roles = _cosmic_labels(sub, cosmic_map)
    g.set_yticklabels(labels, rotation=0, fontsize=9)
    for tick, role in zip(g.get_yticklabels(), roles):
        tick.set_color(COSMIC_COLORS.get(role, 'tab:grey'))
    # Tissue only on the ticks (rotated -- 'Bone marrow' collides at 0 deg);
    # the ABE | KO split is named once per half, under the axis.
    g.set_xticklabels([TISSUE_LABEL[t] for _a in ('ABE', 'KO')
                       for t in tissues],
                      rotation=45, ha='right', rotation_mode='anchor',
                      fontsize=9)
    for k, arm in enumerate(('ABE', 'KO')):
        ax.text((k + .5) / 2, -.028 - 1.55 / max(len(sub), 3), arm, ha='center', va='top', fontsize=11,
                fontweight='bold', transform=ax.transAxes, clip_on=False)
    ax.axvline(len(tissues), color='black', lw=2)   # ABE | KO divider
    ax.set_title(f'{title}  (n = {len(sub)})', fontsize=11)

    handles = [Patch(facecolor=COSMIC_COLORS[r], label=r)
               for r in COSMIC_ORDER if r in set(roles)]
    if handles:
        ax.legend(handles=handles, loc='upper left', bbox_to_anchor=(1.42, 1),
                  fontsize=9, frameon=False, title='COSMIC role',
                  title_fontsize=9)

    fig.tight_layout()
    fig.savefig(OUT + fname + '.pdf', bbox_inches='tight')
    fig.savefig(OUT + fname + '.png', dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f'  wrote {fname}  (n = {len(sub)})')


def panel_f_concordant(m, cosmic_map, tag='', tissues=None):
    tissues = m.attrs.get('tissues', HIT_TISSUES) if tissues is None else tissues
    where = '' if len(tissues) > 1 else f' ({TISSUE_LABEL[tissues[0]]})'
    sets = concordant_sets(m)
    sort_t = tissues[0]
    _heatmap(sets['UP in both'], f'Enriched in both ABE and Cas9-NG{where}',
             f'KO_f1_up_in_both{tag}', cosmic_map,
             f'LFC_median_{sort_t}_ABE', tissues)
    _heatmap(sets['DOWN in both'].head(25),
             f'Depleted in both ABE and Cas9-NG{where}',
             f'KO_f2_down_in_both{tag}', cosmic_map,
             f'LFC_median_{sort_t}_KO', tissues)
    return sets


# ══════════════════════════════════════════════════════════════════════════
# g  Most-divergent heatmap (Sam's ranking, restricted to the discordant
#    spleen sets)
#
# Ranking is Sam's original (extended_figure11_CLEANED.ipynb cell 4):
#     delta_spleen = LFC_median_spleen_ABE - LFC_median_spleen_KO
#     sort descending, take the top N.
# Sam's version selected on LFC_median_spleen >= 1 alone, with no FDR and
# no KO-side condition, so it could rank a guide whose KO arm was merely
# noisy. Here the candidate pool is instead one of the two discordant
# quadrants defined in hit_frame()/concordant_sets():
#     'ABE up / KO down'  = ABE_spleen_up & KO_spleen_dn
#     'ABE down / KO up'  = ABE_spleen_dn & KO_spleen_up
# Both arms must therefore clear |LFC| >= LFC_CUTOFF with the matching
# directional FDR < FDR_CUTOFF, and both must pass Input_median >= MIN_READS.
# The ordering within the panel is still purely Sam's delta.
# ══════════════════════════════════════════════════════════════════════════
def divergent_sets(m):
    """The two discordant quadrants, spleen, with Sam's delta attached."""
    d = m['LFC_median_spleen_ABE'] - m['LFC_median_spleen_KO']
    m = m.assign(delta_spleen=d)
    return {'ABE up / KO down': m[m['ABE_spleen_up'] & m['KO_spleen_dn']],
            'ABE down / KO up': m[m['ABE_spleen_dn'] & m['KO_spleen_up']]}


def _sam_heatmap(sub, fname, cosmic_map):
    """Sam's heatmap (extended_figure11_CLEANED.ipynb cell 4).

    Departures from his code, both requested:
      * y-tick labels coloured by COSMIC role
      * larger fonts on the mutation labels and the colourbar ticks
    Colour scale is Sam's: seismic, auto-scaled to the panel's own range.
    """
    if len(sub) == 0:
        print(f'  [skip] {fname}: no guides in this set')
        return

    cols = ['LFC_median_spleen_ABE', 'LFC_median_spleen_KO']

    fig, ax = plt.subplots(figsize=(5, 9))
    # Sam's call: no vmin/vmax, so seaborn auto-scales to the panel's own
    # min/max. Zero is then off-centre, which is what his figure shows.
    g = sns.heatmap(sub[cols], cmap='seismic',
                    annot=True, annot_kws={'size': 14},
                    linewidth=1, edgecolor='white', ax=ax)

    name, roles = [], []
    for _i, val in sub.iterrows():
        a = val['gene_name_h_ABE']
        b = val['HGVSp_h_ABE']
        c = val['gene_name_m_corrected_ABE']
        d = val['HGVSp_m_ABE']
        if type(b) != str:
            b = 'Non-coding'
        if type(d) != str:
            d = 'Non-coding'
        name.append(f'{a} {b}\n{c} {d}')
        roles.append(cosmic_map.get(a, 'Undefined'))

    g.set_yticklabels(name, rotation=0, fontsize=14)
    for tick, role in zip(g.get_yticklabels(), roles):
        tick.set_color(COSMIC_COLORS.get(role, 'tab:grey'))
    g.set_xticklabels(['ABE', 'KO'], rotation=0, fontsize=16)
    g.collections[0].colorbar.ax.tick_params(labelsize=14)
    fig.tight_layout()
    fig.savefig(OUT + fname + '.pdf', transparent=True, bbox_inches='tight')
    fig.savefig(OUT + fname + '.png', dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f'  wrote {fname}  (n = {len(sub)})')


def panel_g_divergent(m, cosmic_map, top_n=15, tag='_spleen'):
    """Sam's most-divergent heatmap, pool restricted to the discordant sets."""
    sets = divergent_sets(m)

    ud = sets['ABE up / KO down'].sort_values('delta_spleen', ascending=False)
    _sam_heatmap(ud.head(top_n), f'KO_g1_divergent_up_dn{tag}', cosmic_map)

    du = sets['ABE down / KO up'].sort_values('delta_spleen', ascending=True)
    _sam_heatmap(du.head(top_n), f'KO_g2_divergent_dn_up{tag}', cosmic_map)
    return sets


def panel_f_up_relaxed(cosmic_map, tag='', tissues=None):
    """Why 'UP in both' is empty, made explicit.

    At the standard thresholds no guide is enriched in both arms. The
    reason is not that the cutoffs are too strict in general -- it is that
    EVERY KO guide calling as enriched is a low-abundance guide. Dropping
    the Input_median floor to 0 lets them in, and this panel shows what
    they look like: input medians of ~11-32 reads, where a handful of
    counts in the numerator produces a large positive LFC. Show this panel
    next to the depleted one so the asymmetry is visible, and do not read
    these as hits.
    """
    tissues = HIT_TISSUES if tissues is None else list(tissues)
    where = '' if len(tissues) > 1 else f' ({TISSUE_LABEL[tissues[0]]})'
    sort_t = tissues[0]
    m = hit_frame(min_reads=0, tissues=tissues)
    sub = m[m['ABE_up'] & m['KO_up']]
    if len(sub) == 0:
        print('  [skip] UP in both: empty even with no read floor')
        return sub
    sub = sub.sort_values(by=f'LFC_median_{sort_t}_KO',
                          ascending=False).head(25)
    _heatmap(sub,
             f'Enriched in both{where} -- NO read-depth floor (artifact)',
             f'KO_f1_up_in_both_relaxed{tag}', cosmic_map,
             f'LFC_median_{sort_t}_KO', tissues)
    print('  input medians of these guides: '
          f"ABE {sub['Input_median_ABE'].median():.0f}, "
          f"KO {sub['Input_median_KO'].median():.0f} "
          '(cutoff for a real hit is 100)')
    return sub


def cosmic_role_map():
    """gene symbol (HUMAN) -> one of COSMIC_ORDER.

    Same recipe as extended_figure9_CLEANED cell 5: fillna('Undefined')
    BEFORE the lookup, and join on the human symbol.
    """
    cos = pd.read_csv(SG_BASE +
                      'source_data/Census_allSun Nov 17 02_26_47 2024.csv'
                      ).fillna('Undefined')
    t_dict = {'TSG': 'TSG', 'TSG, fusion': 'TSG', 'Undefined': 'Undefined',
              'fusion': 'Undefined', 'oncogene': 'Oncogene',
              'oncogene, TSG': 'Oncogene/TSG',
              'oncogene, TSG, fusion': 'Oncogene/TSG',
              'oncogene, fusion': 'Oncogene'}
    return {r['Gene Symbol']: t_dict[r['Role in Cancer']]
            for _i, r in cos.iterrows()}



# ══════════════════════════════════════════════════════════════════════════
# c_single  The published Fig 3b layout as a standalone one-row figure.
#
# panel_c draws this as the top row of a 5-row (All + 4 COSMIC roles)
# figure. This writes the same row on its own, so it can be dropped into
# the manuscript without cropping a multi-row PDF.
# ══════════════════════════════════════════════════════════════════════════
def panel_c_single(m, filter_edit=True):
    names = ['In Vitro', 'Spleen', 'Bone', 'Meninges']

    # Taller than panel_c's row so each panel reads close to square: the
    # data range is 8 x 8, so the axes box should be too. set_aspect(1)
    # below enforces it regardless of the final figure size.
    fig, axs = plt.subplots(1, len(TISSUES), figsize=(10, 4.2),
                            sharex=True, sharey=True)

    sub = m[(m['classification'] == 'targeting guide')
            & (m['target_base_edit_perc'] >= MIN_EDIT)]
    if filter_edit:
        sub = sub[sub['target_indel_perc'] >= MIN_EDIT]

    for j, t in enumerate(TISSUES):
        ax = axs[j]
        x = sub[f'LFC_median_{t}'].to_numpy(dtype=float)       # ABE
        y = sub[f'LFC_median_{t}_KO'].to_numpy(dtype=float)    # KO
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
        ax.set_xlim(-8, 8)
        ax.set_ylim(-8, 8)
        ax.set_aspect('equal', adjustable='box')

        spear = np.round(scipy.stats.spearmanr(x, y).statistic, 2)
        pear = np.round(scipy.stats.pearsonr(x, y).statistic, 2)
        ax.set_title(f'{names[j]}\n'
                     + '$R_{p} = $' + f'{pear}'
                     + ' | $R_{s} = $' + f'{spear}'
                     + f'\nn = {len(x)}', fontsize=13)
        ax.set_xlabel('LFC ABE', fontsize=13)
        if j == 0:
            ax.set_ylabel('LFC CRISPR KO', fontsize=13)

    tag = 'edit_filtered' if filter_edit else 'unfiltered'
    fig.tight_layout()
    fig.savefig(OUT + f'KO_c_single_LFC_corr_{tag}.pdf',
                bbox_inches='tight', transparent=True)
    fig.savefig(OUT + f'KO_c_single_LFC_corr_{tag}.png',
                dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f'  wrote KO_c_single_LFC_corr_{tag}  (n = {len(x)})')


if __name__ == '__main__':
    m = load()
    print("merged rows:", len(m))
    print("\nRole in Cancer (targeting guides):")
    print(m[m['classification'] == 'targeting guide']['Role in Cancer']
          .value_counts().to_string())
    print("\nguides with sensor indel:", int(m['target_indel_perc'].notna().sum()))

    panel_a(m)
    panel_a2(m)
    panel_b(m)
    panel_c(m, filter_edit=True)
    panel_c(m, filter_edit=False)
    panel_c_single(m, filter_edit=True)
    panel_c_single(m, filter_edit=False)
    panel_d(m)

    # ── co-author request: UpSet + concordant-set heatmaps ────────────────
    cosmic_map = cosmic_role_map()
    h = hit_frame()
    print("\nhit frame:", len(h), "guides (targeting, ABE edit >= "
          f"{MIN_EDIT:g}%, Input_median >= {MIN_READS} on both arms)")
    for arm in ('ABE', 'KO'):
        for d in ('up', 'dn'):
            print(f"  {arm}_{d}: {int(h[f'{arm}_{d}'].sum())}")
    panel_e_upset(h)
    panel_f_concordant(h, cosmic_map)
    panel_f_up_relaxed(cosmic_map)

    # ── same analysis, spleen only ───────────────────────────────────────
    # A guide counts as a hit on its spleen LFC/FDR alone, so the sets are
    # not just the pooled ones re-coloured: pooling over three tissues
    # calls anything significant in ANY of them.
    hs = hit_frame(tissues=['spleen'])
    print("\nspleen only:", len(hs), "guides")
    for arm in ('ABE', 'KO'):
        for d in ('up', 'dn'):
            print(f"  {arm}_{d}: {int(hs[f'{arm}_{d}'].sum())}")
    panel_e_upset(hs, tag='_spleen')
    panel_f_concordant(hs, cosmic_map, tag='_spleen')
    panel_f_up_relaxed(cosmic_map, tag='_spleen', tissues=['spleen'])

    # ── most-divergent heatmap, Sam's ranking inside the discordant sets ──
    dsets = panel_g_divergent(hs, cosmic_map)
    for k, v in dsets.items():
        print(f"  divergent '{k}': {len(v)} guides")

    print("\nwrote panels to", OUT)
