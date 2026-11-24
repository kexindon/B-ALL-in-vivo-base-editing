"""
Proposed new analyses / visualizations for the sensor-calibration paragraph.
B-ALL in vivo base editing. Draft panels for review — NOT final figures.

Outputs one main figure + 4 extended figures as single-page PDFs, then a
compacted review PDF.
"""
import os, warnings
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns
import scipy.stats
from adjustText import adjust_text
warnings.filterwarnings('ignore')

# ── Global style: final print size, Nature look ───────────────────────────
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['text.usetex'] = False
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Helvetica', 'Arial', 'DejaVu Sans']
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
teals = sns.color_palette('BuPu').as_hex()

BASE = '/Users/kexindong/Documents/GitHub/B-ALL-in-vivo-base-editing/'
OUT = os.path.dirname(os.path.abspath(__file__))
os.makedirs(OUT, exist_ok=True)

# Fixed semantic colors (sci-plot conventions)
TISSUE_COLORS = {'Plasmid': 'black', 'Input': yellows[3], 'D5': purples[2],
                 'D15': purples[5], 'Spleen': blues[3], 'Bone marrow': greens[3],
                 'Meninges': reds[3]}
EDITOR_COLORS = {'ABE': 'tab:purple', 'CBE': 'tab:blue'}
EDIT_TYPE_COLORS = {'target': blues[4], 'correct': teals[2]}

TISSUES = ['d5', 'd15', 'spleen', 'bm', 'men']
TISSUE_LABEL = {'d5': 'D5', 'd15': 'D15', 'spleen': 'Spleen',
                'bm': 'Bone marrow', 'men': 'Meninges'}
INVIVO = ['spleen', 'bm', 'men']
REPS = {'bm': ['bonemarrow1', 'bonemarrow2', 'bonemarrow3', 'bonemarrow4', 'bonemarrow5'],
        'spleen': ['spleen1', 'spleen2', 'spleen3', 'spleen4'],
        'men': ['meninges1', 'meninges2', 'meninges3', 'meninges4', 'meninges5'],
        'd15': ['d15_rep1', 'd15_rep2', 'd15_rep3'],
        'd5': ['d5_rep1', 'd5_rep2', 'd5_rep3']}

print('Loading…')
ABE = pd.read_csv(BASE + 'LFC-FDR/ABE_OG_LFC_FDR.csv', low_memory=False)
CBE = pd.read_csv(BASE + 'LFC-FDR/CBE_OG_LFC_FDR.csv', low_memory=False)
ABE_C = pd.read_csv(BASE + 'counts/og/OG_ABE_COUNTS.txt', sep='\t')
CBE_C = pd.read_csv(BASE + 'counts/og/OG_CBE_COUNTS.txt', sep='\t')
DFS = {'ABE': ABE, 'CBE': CBE}
CNT = {'ABE': ABE_C, 'CBE': CBE_C}

MIN_READS = 100
FDR_CUT = 0.1


SHORT = {'d5': 'D5', 'd15': 'D15', 'spleen': 'Spl', 'bm': 'BM', 'men': 'Men'}


def targ(d):
    return d[d['classification'] == 'targeting guide']


def ctrls(d):
    return d[d['classification'].isin(['safe-targeting control', 'non-targeting control'])]


def open_ax(ax):
    ax.spines[['top', 'right']].set_visible(False)
    ax.tick_params(axis='both', which='major', labelsize=6,
                   direction='out', length=3, width=.6)


def panel_letter(ax, letter, dx=-0.085, dy=1.06):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=8,
            fontweight='bold', va='bottom', ha='left')


def gini(x):
    x = np.sort(np.asarray(x, float))
    x = x[~np.isnan(x)]
    n = len(x)
    if n == 0 or x.sum() == 0:
        return np.nan
    return (((2 * np.arange(1, n + 1) - n - 1) * x).sum()) / (n * x.sum())


def eff_n(x):
    """Shannon effective number of guides."""
    p = np.asarray(x, float)
    p = p[p > 0]
    p = p / p.sum()
    return float(np.exp(-(p * np.log(p)).sum()))


def sample_group(col):
    c = col.lower()
    if 'plasmid' in c:
        return 'Plasmid'
    if 'input' in c:
        return 'Input'
    if c.startswith('d5'):
        return 'D5'
    if c.startswith('d15'):
        return 'D15'
    if 'spleen' in c:
        return 'Spleen'
    if 'bonemarrow' in c:
        return 'Bone marrow'
    if 'meninges' in c:
        return 'Meninges'
    return None


# ══════════════════════════════════════════════════════════════════════════
# Shared computations
# ══════════════════════════════════════════════════════════════════════════
THRS = np.arange(0, 95, 5)


def sweep_table(editor, metric='target_base_edit_perc'):
    """Signal-to-noise + replicate concordance + n retained, vs editing threshold."""
    d = DFS[editor]
    t, c = targ(d), ctrls(d)
    rows = []
    for s in TISSUES:
        rc = [r for r in REPS[s] if r in d.columns]
        null_sd = c[f'LFC_median_{s}'].std()
        for thr in THRS:
            sel = t[(t[f'{metric}_{s}'] >= thr) & (t[f'sensor_reads_{s}'] >= MIN_READS)]
            if len(sel) < 20:
                continue
            snr = sel[f'LFC_median_{s}'].std() / null_sd
            sub = sel[rc].dropna()
            if len(sub) > 10 and len(rc) > 1:
                cm = sub.corr(method='spearman').values
                rep_r = np.nanmean(cm[np.triu_indices_from(cm, 1)])
            else:
                rep_r = np.nan
            fdrcol = f'FDR_{s}_stouffers'
            frac_sig = (sel[fdrcol] < FDR_CUT).mean() * 100 if fdrcol in sel else np.nan
            rows.append(dict(editor=editor, tissue=s, thr=thr, n=len(sel),
                             snr=snr, rep_r=rep_r, frac_sig=frac_sig,
                             frac_retained=len(sel) / len(t) * 100))
    return pd.DataFrame(rows)


print('Computing threshold sweeps…')
SWEEP = pd.concat([sweep_table('ABE'), sweep_table('CBE')], ignore_index=True)


def bottleneck_table():
    rows = []
    for editor, d in CNT.items():
        n_lib = len(d)
        for col in d.columns[2:]:
            grp = sample_group(col)
            if grp is None:
                continue
            v = d[col].astype(float)
            if v.sum() == 0:
                continue
            cpm = v / v.sum() * 1e6
            rows.append(dict(editor=editor, sample=col, group=grp,
                             gini=gini(cpm),
                             skew=np.percentile(cpm, 90) / max(np.percentile(cpm, 10), 1e-9),
                             eff_frac=eff_n(cpm) / n_lib * 100,
                             n_lib=n_lib))
    return pd.DataFrame(rows)


BOTT = bottleneck_table()
print('Bottleneck table:', BOTT.shape)


# ══════════════════════════════════════════════════════════════════════════
# MAIN FIGURE — Sensor calibration as a tunable analytical lever
# ══════════════════════════════════════════════════════════════════════════
def main_figure(return_fig=False):
    fig = plt.figure(figsize=(183 * MM, 170 * MM))
    gs = fig.add_gridspec(3, 3, hspace=.95, wspace=.5,
                          left=.08, right=.97, top=.89, bottom=.07)

    # ── a: editing-outcome composition, stacked, per tissue, both editors ──
    ax = fig.add_subplot(gs[0, 0])
    groups, lefts = [], []
    for editor in ['ABE', 'CBE']:
        t = targ(DFS[editor])
        for s in ['d15', 'bm', 'men']:
            ok = t[f'sensor_reads_{s}'] >= MIN_READS
            tb = t.loc[ok, f'target_base_edit_perc_{s}']
            cp = t.loc[ok, f'corr_perc_{s}']
            groups.append(dict(lab=f'{editor}\n{SHORT[s]}',
                               correct=cp.median(),
                               byst=max(tb.median() - cp.median(), 0),
                               unedited=max(100 - tb.median(), 0)))
    gd = pd.DataFrame(groups)
    x = np.arange(len(gd))
    ax.bar(x, gd['correct'], color=EDIT_TYPE_COLORS['correct'], edgecolor='black',
           linewidth=.5, label='Perfect edit')
    ax.bar(x, gd['byst'], bottom=gd['correct'], color=EDIT_TYPE_COLORS['target'],
           edgecolor='black', linewidth=.5, label='Target + bystander')
    ax.bar(x, gd['unedited'], bottom=gd['correct'] + gd['byst'], color=greys[1],
           edgecolor='black', linewidth=.5, label='Unedited')
    ax.set_xticks(x)
    ax.set_xticklabels(gd['lab'], fontsize=5)
    ax.set_ylabel('Sensor reads (%, median)')
    ax.set_ylim(0, 100)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_yticklabels(['0', '25', '50', '75', '100%'])
    open_ax(ax)
    ax.legend(loc='lower center', bbox_to_anchor=(.5, 1.03), ncol=3,
              fontsize=4.5, frameon=False, handlelength=.9, handleheight=.9,
              columnspacing=.6, handletextpad=.4)
    panel_letter(ax, 'a')

    # ── b: editing purity (perfect / target) — ABE vs CBE ──────────────────
    ax = fig.add_subplot(gs[0, 1])
    rows = []
    for editor in ['ABE', 'CBE']:
        t = targ(DFS[editor])
        for s in INVIVO:
            ok = (t[f'sensor_reads_{s}'] >= MIN_READS) & (t[f'target_base_edit_perc_{s}'] >= 20)
            pur = (t.loc[ok, f'corr_perc_{s}'] / t.loc[ok, f'target_base_edit_perc_{s}'] * 100)
            for v in pur.dropna():
                rows.append(dict(editor=editor, tissue=TISSUE_LABEL[s], purity=v))
    pur_df = pd.DataFrame(rows)
    sns.boxplot(data=pur_df, x='tissue', y='purity', hue='editor', ax=ax,
                palette=EDITOR_COLORS, showfliers=False, linewidth=.6,
                width=.65, legend=False)
    ax.set_ylabel('Editing purity\n(perfect / target, %)')
    ax.set_xlabel('')
    ax.set_ylim(0, 100)
    ax.tick_params(axis='x', labelsize=5)
    open_ax(ax)
    h = [Patch(facecolor=EDITOR_COLORS[e], edgecolor='black', lw=.5, label=e)
         for e in ['ABE', 'CBE']]
    ax.legend(handles=h, loc='lower center', bbox_to_anchor=(.5, 1.02),
              ncol=2, fontsize=6, frameon=False, handlelength=1, handleheight=1)
    panel_letter(ax, 'b')

    # ── c: per-guide editing reproducibility across contexts (ABE) ─────────
    ax = fig.add_subplot(gs[0, 2])
    t = targ(ABE)
    cols = [f'target_base_edit_perc_{s}' for s in TISSUES]
    ok = (t[[f'sensor_reads_{s}' for s in TISSUES]] >= MIN_READS).all(axis=1)
    sub = t.loc[ok, cols].dropna()
    cm = sub.corr(method='spearman')
    sns.heatmap(cm, ax=ax, cmap='Reds', vmin=0, vmax=1, annot=True, fmt='.2f',
                annot_kws={'size': 5}, linewidth=1, linecolor='white',
                cbar_kws={'label': 'Spearman ρ', 'shrink': .7})
    ax.set_xticklabels([TISSUE_LABEL[s] for s in TISSUES], fontsize=5, rotation=45,
                       ha='right', rotation_mode='anchor')
    ax.set_yticklabels([TISSUE_LABEL[s] for s in TISSUES], fontsize=5, rotation=0)
    ax.set_title(f'ABE per-gRNA target editing\n(n = {len(sub)} gRNAs)', fontsize=6)
    cb = ax.collections[0].colorbar
    cb.ax.tick_params(labelsize=5)
    cb.set_label('Spearman ρ', size=6)
    panel_letter(ax, 'c', dx=-0.18)

    # ── d: THE ANALOG LEVER — SNR vs editing threshold ─────────────────────
    ax = fig.add_subplot(gs[1, 0])
    sw = SWEEP[(SWEEP.editor == 'ABE')]
    for s in TISSUES:
        ss = sw[sw.tissue == s]
        ax.plot(ss.thr, ss.snr, '-o', color=TISSUE_COLORS[TISSUE_LABEL[s]],
                ms=2.5, lw=1, mec='black', mew=.3, label=TISSUE_LABEL[s])
    ax.axhline(1, color='black', ls='--', lw=.8)
    ax.axvline(20, color=greys[3], ls=':', lw=.8)
    ax.text(21, ax.get_ylim()[1] * .98, '20% cutoff', fontsize=5,
            color=greys[4], va='top')
    ax.set_xlabel('Target base editing cutoff (%)')
    ax.set_ylabel('Signal-to-noise\n(SD targeting / SD control)')
    open_ax(ax)
    ax.legend(loc='lower center', bbox_to_anchor=(.5, 1.02), ncol=3,
              fontsize=5, frameon=False, handlelength=1)
    panel_letter(ax, 'd')

    # ── e: replicate concordance vs threshold — ABE vs CBE ─────────────────
    ax = fig.add_subplot(gs[1, 1])
    for editor, ls in [('ABE', '-'), ('CBE', '--')]:
        sw = SWEEP[SWEEP.editor == editor]
        for s in ['bm', 'spleen', 'men']:
            ss = sw[sw.tissue == s]
            ax.plot(ss.thr, ss.rep_r, ls, color=TISSUE_COLORS[TISSUE_LABEL[s]],
                    lw=1, marker='o' if editor == 'ABE' else '^', ms=2, mec='black',
                    mew=.3)
    ax.set_xlabel('Target base editing cutoff (%)')
    ax.set_ylabel('Mean inter-replicate ρ\n(gRNA LFC)')
    ax.set_ylim(0, .8)
    open_ax(ax)
    h = [plt.Line2D([], [], color='black', ls='-', marker='o', ms=2.5, lw=1, label='ABE'),
         plt.Line2D([], [], color='black', ls='--', marker='^', ms=2.5, lw=1, label='CBE')]
    h += [Patch(facecolor=TISSUE_COLORS[TISSUE_LABEL[s]], edgecolor='none',
                label=TISSUE_LABEL[s]) for s in ['spleen', 'bm', 'men']]
    ax.legend(handles=h, loc='lower center', bbox_to_anchor=(.5, 1.02), ncol=3,
              fontsize=5, frameon=False, handlelength=1.2, columnspacing=.8)
    panel_letter(ax, 'e')

    # ── f: the cost of the lever — guides retained vs threshold ────────────
    ax = fig.add_subplot(gs[1, 2])
    for editor, ls, mk in [('ABE', '-', 'o'), ('CBE', '--', '^')]:
        sw = SWEEP[(SWEEP.editor == editor) & (SWEEP.tissue == 'bm')]
        ax.plot(sw.thr, sw.frac_retained, ls, marker=mk, color=EDITOR_COLORS[editor],
                ms=2.5, lw=1, mec='black', mew=.3, label=editor)
    ax.axvline(20, color=greys[3], ls=':', lw=.8)
    ax.set_xlabel('Target base editing cutoff (%)')
    ax.set_ylabel('gRNAs retained (%)')
    ax.set_ylim(0, 100)
    ax.set_title('Bone marrow', fontsize=6)
    open_ax(ax)
    ax.legend(loc='upper right', fontsize=6, frameon=False, handlelength=1.2)
    panel_letter(ax, 'f')

    # ── g: calibration changes the hit list — LFC vs editing, BM ABE ───────
    ax = fig.add_subplot(gs[2, 0])
    t = targ(ABE)
    ok = t['sensor_reads_bm'] >= MIN_READS
    tt = t[ok]
    sig = tt['FDR_bm_stouffers'] < FDR_CUT
    ax.scatter(tt.loc[~sig, 'target_base_edit_perc_bm'], tt.loc[~sig, 'LFC_median_bm'],
               s=5, color='lightgrey', linewidths=0, alpha=.6, label='n.s.')
    enr = sig & (tt['LFC_median_bm'] > 0)
    dep = sig & (tt['LFC_median_bm'] < 0)
    ax.scatter(tt.loc[enr, 'target_base_edit_perc_bm'], tt.loc[enr, 'LFC_median_bm'],
               s=9, color=reds[2], edgecolors='black', linewidths=.2, label='Enriched')
    ax.scatter(tt.loc[dep, 'target_base_edit_perc_bm'], tt.loc[dep, 'LFC_median_bm'],
               s=9, color=blues[2], edgecolors='black', linewidths=.2, label='Depleted')
    ax.axhline(0, color='black', lw=.8, ls='--')
    ax.axvspan(0, 20, color=greys[1], alpha=.45, lw=0)
    ax.text(10, ax.get_ylim()[1] * .9, 'excluded', fontsize=5, ha='center',
            color=greys[4])
    ax.set_xlabel('Target base editing (%), bone marrow')
    ax.set_ylabel('gRNA LFC (median)')
    open_ax(ax)
    ax.legend(loc='lower center', bbox_to_anchor=(.5, 1.02), ncol=3, fontsize=5,
              frameon=False, handlelength=1)
    panel_letter(ax, 'g')

    # ── h: effective library coverage — the bottleneck ─────────────────────
    ax = fig.add_subplot(gs[2, 1])
    order = ['Plasmid', 'Input', 'D5', 'D15', 'Spleen', 'Bone marrow', 'Meninges']
    bb = BOTT[BOTT.group.isin(order)].copy()
    bb['group'] = pd.Categorical(bb['group'], order, ordered=True)
    xs = np.arange(len(order))
    w = .36
    for i, editor in enumerate(['ABE', 'CBE']):
        sub = bb[bb.editor == editor]
        m = sub.groupby('group', observed=True)['eff_frac'].mean().reindex(order)
        sd = sub.groupby('group', observed=True)['eff_frac'].std().reindex(order)
        pos = xs + (i - .5) * w
        ax.bar(pos, m, w, yerr=sd, color=EDITOR_COLORS[editor], edgecolor='black',
               linewidth=.5, error_kw=dict(lw=.7, capsize=1.5, ecolor='black'),
               label=editor, zorder=1)
        for g, p in zip(order, pos):
            v = sub[sub.group == g]['eff_frac'].values
            if len(v):
                ax.scatter(np.full(len(v), p) + np.random.uniform(-.08, .08, len(v)),
                           v, s=3, color='black', alpha=.55, linewidths=0, zorder=3)
    ax.set_xticks(xs)
    ax.set_xticklabels(order, fontsize=5, rotation=45, ha='right',
                       rotation_mode='anchor')
    ax.set_ylabel('Effective library coverage\n(% of library)')
    ax.set_ylim(0, 100)
    open_ax(ax)
    ax.legend(loc='lower center', bbox_to_anchor=(.5, 1.02), ncol=2, fontsize=6,
              frameon=False, handlelength=1, handleheight=1)
    panel_letter(ax, 'h')

    # ── i: library size vs recoverable signal ──────────────────────────────
    ax = fig.add_subplot(gs[2, 2])
    _texts_i = []
    for editor, mk in [('ABE', 'o'), ('CBE', '^')]:
        sw = SWEEP[(SWEEP.editor == editor) & (SWEEP.thr == 20)]
        bb = BOTT[BOTT.editor == editor]
        for s in INVIVO:
            lab = TISSUE_LABEL[s]
            grp = {'spleen': 'Spleen', 'bm': 'Bone marrow', 'men': 'Meninges'}[s]
            ef = bb[bb.group == grp]['eff_frac'].mean()
            rr = sw[sw.tissue == s]['rep_r'].values
            if len(rr) == 0:
                continue
            ax.scatter(ef, rr[0], s=38, marker=mk, color=TISSUE_COLORS[lab],
                       edgecolors='black', linewidths=.5, zorder=3)
            _texts_i.append(ax.text(ef, rr[0], f'{editor[0]}·{lab.split()[0]}',
                                    fontsize=5))
    ax.set_xlabel('Effective library coverage (%)')
    ax.set_ylabel('Mean inter-replicate ρ\n(at 20% cutoff)')
    ax.set_xlim(0, 100)
    ax.set_ylim(0, .8)
    open_ax(ax)
    adjust_text(_texts_i, ax=ax, expand=(1.3, 1.6),
                arrowprops=dict(arrowstyle='-', color='black', lw=.4))
    panel_letter(ax, 'i')

    fig.suptitle('PROPOSED MAIN FIGURE — Sensor-informed calibration as a tunable analytical lever',
                 fontsize=8, fontweight='bold', y=1.0)
    fig.savefig(f'{OUT}/main_figure.pdf', bbox_inches='tight')
    fig.savefig(f'{OUT}/main_figure.png', dpi=150, bbox_inches='tight')
    print('  main_figure done')
    if return_fig:
        return fig
    plt.close(fig)


# ══════════════════════════════════════════════════════════════════════════
# EXTENDED 1 — Full threshold-sweep atlas
# ══════════════════════════════════════════════════════════════════════════
def ext1(return_fig=False):
    fig, axs = plt.subplots(2, 3, figsize=(183 * MM, 105 * MM))
    for r, editor in enumerate(['ABE', 'CBE']):
        sw = SWEEP[SWEEP.editor == editor]
        for c, (met, lab) in enumerate([
                ('snr', 'Signal-to-noise\n(SD targeting / SD control)'),
                ('frac_sig', f'gRNAs with FDR < {FDR_CUT} (%)'),
                ('frac_retained', 'gRNAs retained (%)')]):
            ax = axs[r, c]
            for s in TISSUES:
                ss = sw[sw.tissue == s]
                ax.plot(ss.thr, ss[met], '-o', color=TISSUE_COLORS[TISSUE_LABEL[s]],
                        ms=2, lw=.9, mec='black', mew=.25, label=TISSUE_LABEL[s])
            ax.axvline(20, color=greys[3], ls=':', lw=.8)
            if met == 'snr':
                ax.axhline(1, color='black', ls='--', lw=.8)
            ax.set_xlabel('Target base editing cutoff (%)')
            ax.set_ylabel(lab)
            ax.set_title(f'{editor}', fontsize=7)
            open_ax(ax)
            if r == 0 and c == 0:
                ax.legend(loc='lower center', bbox_to_anchor=(.5, 1.12), ncol=5,
                          fontsize=5, frameon=False, handlelength=1,
                          columnspacing=.7)
    for ax, L in zip(axs.ravel(), 'abcdef'):
        panel_letter(ax, L, dx=-0.22, dy=1.04)
    fig.suptitle('EXTENDED 1 — Threshold sweep atlas: the sensor cutoff as a continuous lever',
                 fontsize=8, fontweight='bold', y=1.0)
    fig.tight_layout(rect=[0, 0, 1, .95])
    fig.savefig(f'{OUT}/ext1.pdf', bbox_inches='tight')
    fig.savefig(f'{OUT}/ext1.png', dpi=150, bbox_inches='tight')
    print('  ext1 done')
    if return_fig:
        return fig
    plt.close(fig)


# ══════════════════════════════════════════════════════════════════════════
# EXTENDED 2 — Perfect vs target editing: metric choice
# ══════════════════════════════════════════════════════════════════════════
def ext2(return_fig=False):
    fig, axs = plt.subplots(2, 3, figsize=(183 * MM, 110 * MM))

    # a/d: target vs perfect scatter
    for r, editor in enumerate(['ABE', 'CBE']):
        t = targ(DFS[editor])
        ok = t['sensor_reads_d15'] >= MIN_READS
        x = t.loc[ok, 'target_base_edit_perc_d15']
        y = t.loc[ok, 'corr_perc_d15']
        ax = axs[r, 0]
        hb = ax.hexbin(x, y, gridsize=34, bins='log', cmap='inferno', mincnt=1,
                       linewidths=0)
        ax.plot([0, 100], [0, 100], '--', color='black', lw=.8)
        rs = scipy.stats.spearmanr(x, y).statistic
        rp = scipy.stats.pearsonr(x, y).statistic
        ax.set_title(f'{editor} D15 — ρ = {rs:.2f} and R = {rp:.2f}', fontsize=6)
        ax.set_xlabel('Target base editing (%)')
        ax.set_ylabel('Perfect base editing (%)')
        ax.set_xlim(0, 100)
        ax.set_ylim(0, 100)
        open_ax(ax)
        cb = fig.colorbar(hb, ax=ax, shrink=.75, pad=.02)
        cb.ax.tick_params(labelsize=5)
        cb.set_label('gRNAs (log)', size=5)

        # b/e: purity distribution
        ax = axs[r, 1]
        for s in INVIVO:
            ok2 = (t[f'sensor_reads_{s}'] >= MIN_READS) & (t[f'target_base_edit_perc_{s}'] >= 20)
            pur = (t.loc[ok2, f'corr_perc_{s}'] / t.loc[ok2, f'target_base_edit_perc_{s}'] * 100)
            ax.hist(pur.dropna(), bins=np.linspace(0, 100, 26), histtype='step',
                    lw=1, color=TISSUE_COLORS[TISSUE_LABEL[s]],
                    label=TISSUE_LABEL[s], density=True)
        ax.set_xlabel('Editing purity (perfect / target, %)')
        ax.set_ylabel('Density')
        ax.set_title(f'{editor}', fontsize=6)
        open_ax(ax)
        if r == 0:
            ax.legend(fontsize=5, frameon=False, handlelength=1.2)

        # c/f: SNR under the two metrics
        ax = axs[r, 2]
        sw_t = sweep_table(editor, 'target_base_edit_perc')
        sw_c = sweep_table(editor, 'corr_perc')
        for met, sw, ls, mk in [('Target', sw_t, '-', 'o'), ('Perfect', sw_c, '--', 's')]:
            ss = sw[sw.tissue == 'bm']
            ax.plot(ss.thr, ss.snr, ls, marker=mk, ms=2.5, lw=1,
                    color=EDIT_TYPE_COLORS['target' if met == 'Target' else 'correct'],
                    mec='black', mew=.3, label=f'{met} editing cutoff')
        ax.axhline(1, color='black', ls='--', lw=.8)
        ax.set_xlabel('Editing cutoff (%)')
        ax.set_ylabel('Signal-to-noise')
        ax.set_title(f'{editor} bone marrow', fontsize=6)
        open_ax(ax)
        ax.legend(fontsize=5, frameon=False, handlelength=1.4)

    for ax, L in zip(axs.ravel(), 'abcdef'):
        panel_letter(ax, L, dx=-0.24, dy=1.05)
    fig.suptitle('EXTENDED 2 — Target vs perfect base editing: which metric should gate the data?',
                 fontsize=8, fontweight='bold', y=1.0)
    fig.tight_layout(rect=[0, 0, 1, .95])
    fig.savefig(f'{OUT}/ext2.pdf', bbox_inches='tight')
    fig.savefig(f'{OUT}/ext2.png', dpi=150, bbox_inches='tight')
    print('  ext2 done')
    if return_fig:
        return fig
    plt.close(fig)


# ══════════════════════════════════════════════════════════════════════════
# EXTENDED 3 — Bottlenecking: library size vs recoverable biology
# ══════════════════════════════════════════════════════════════════════════
def ext3(return_fig=False):
    fig, axs = plt.subplots(2, 3, figsize=(183 * MM, 110 * MM))
    order = ['Plasmid', 'Input', 'D5', 'D15', 'Spleen', 'Bone marrow', 'Meninges']

    # a: Lorenz curves
    ax = axs[0, 0]
    for editor, ls in [('ABE', '-'), ('CBE', '--')]:
        d = CNT[editor]
        for grp, col in [('Plasmid', None), ('Bone marrow', None), ('Meninges', None)]:
            cands = [c for c in d.columns[2:] if sample_group(c) == grp]
            if not cands:
                continue
            v = np.sort(d[cands[0]].astype(float).values)
            cum = np.cumsum(v) / v.sum() * 100
            ax.plot(np.arange(1, len(v) + 1) / len(v) * 100, cum, ls, lw=1,
                    color=TISSUE_COLORS[grp],
                    label=f'{editor} {grp}' if editor == 'ABE' or grp != 'Plasmid' else None)
    ax.plot([0, 100], [0, 100], ':', color='black', lw=.8)
    ax.set_xlabel('gRNAs, ranked (%)')
    ax.set_ylabel('Cumulative reads (%)')
    ax.set_title('Lorenz curves (solid ABE, dashed CBE)', fontsize=6)
    open_ax(ax)
    ax.legend(fontsize=4.5, frameon=False, handlelength=1.4, loc='upper left')

    # b: Gini by group
    ax = axs[0, 1]
    bb = BOTT[BOTT.group.isin(order)].copy()
    sns.boxplot(data=bb, x='group', y='gini', hue='editor', order=order,
                palette=EDITOR_COLORS, ax=ax, showfliers=False, linewidth=.6,
                width=.7, legend=False)
    sns.stripplot(data=bb, x='group', y='gini', hue='editor', order=order,
                  dodge=True, ax=ax, color='black', size=1.8, alpha=.6,
                  linewidth=0, legend=False)
    ax.set_xticklabels(order, fontsize=5, rotation=45, ha='right',
                       rotation_mode='anchor')
    ax.set_xlabel('')
    ax.set_ylabel('Gini coefficient\n(count inequality)')
    open_ax(ax)
    h = [Patch(facecolor=EDITOR_COLORS[e], edgecolor='black', lw=.5, label=e)
         for e in ['ABE', 'CBE']]
    ax.legend(handles=h, fontsize=6, frameon=False, loc='upper left',
              handlelength=1, handleheight=1)

    # c: effective coverage, absolute guides
    ax = axs[0, 2]
    for editor in ['ABE', 'CBE']:
        sub = BOTT[(BOTT.editor == editor) & (BOTT.group.isin(order))]
        m = sub.groupby('group', observed=True).apply(
            lambda x: (x['eff_frac'] / 100 * x['n_lib']).mean()).reindex(order)
        ax.plot(np.arange(len(order)), m, '-o', color=EDITOR_COLORS[editor],
                ms=3, lw=1, mec='black', mew=.3, label=editor)
    ax.set_yscale('log')
    ax.set_xticks(np.arange(len(order)))
    ax.set_xticklabels(order, fontsize=5, rotation=45, ha='right',
                       rotation_mode='anchor')
    ax.set_ylabel('Effective no. of gRNAs\n(Shannon, log scale)')
    open_ax(ax)
    ax.legend(fontsize=6, frameon=False, handlelength=1.2)

    # d: skew vs replicate concordance
    ax = axs[1, 0]
    _t3 = []
    for editor, mk in [('ABE', 'o'), ('CBE', '^')]:
        d = DFS[editor]
        t = targ(d)
        for s in ['d15', 'spleen', 'bm', 'men']:
            grp = {'d15': 'D15', 'spleen': 'Spleen', 'bm': 'Bone marrow',
                   'men': 'Meninges'}[s]
            sk = BOTT[(BOTT.editor == editor) & (BOTT.group == grp)]['skew'].mean()
            rc = [r for r in REPS[s] if r in d.columns]
            sel = t[(t[f'target_base_edit_perc_{s}'] >= 20) &
                    (t[f'sensor_reads_{s}'] >= MIN_READS)][rc].dropna()
            if len(sel) < 10:
                continue
            cm = sel.corr(method='spearman').values
            rr = np.nanmean(cm[np.triu_indices_from(cm, 1)])
            ax.scatter(sk, rr, s=36, marker=mk, color=TISSUE_COLORS[grp],
                       edgecolors='black', linewidths=.5, zorder=3)
            _t3.append(ax.text(sk, rr, f'{editor[0]}·{grp.split()[0]}', fontsize=5))
    ax.set_xscale('log')
    ax.set_xlabel('90/10 skew ratio (log)')
    ax.set_ylabel('Mean inter-replicate ρ')
    ax.set_ylim(0, .9)
    open_ax(ax)
    adjust_text(_t3, ax=ax, expand=(1.3, 1.6),
                arrowprops=dict(arrowstyle='-', color='black', lw=.4))

    # e: control-guide null width per tissue
    ax = axs[1, 1]
    rows = []
    for editor in ['ABE', 'CBE']:
        c = ctrls(DFS[editor])
        for s in TISSUES:
            rows.append(dict(editor=editor, tissue=TISSUE_LABEL[s],
                             sd=c[f'LFC_median_{s}'].std()))
    nd = pd.DataFrame(rows)
    xs = np.arange(len(TISSUES))
    for i, editor in enumerate(['ABE', 'CBE']):
        sub = nd[nd.editor == editor].set_index('tissue').reindex(
            [TISSUE_LABEL[s] for s in TISSUES])
        ax.bar(xs + (i - .5) * .36, sub['sd'], .36, color=EDITOR_COLORS[editor],
               edgecolor='black', linewidth=.5, label=editor)
    ax.set_xticks(xs)
    ax.set_xticklabels([TISSUE_LABEL[s] for s in TISSUES], fontsize=5,
                       rotation=45, ha='right', rotation_mode='anchor')
    ax.set_ylabel('Control-gRNA LFC SD\n(noise floor)')
    open_ax(ax)
    ax.legend(fontsize=6, frameon=False, handlelength=1, handleheight=1)

    # f: dropout — zero/low-count guides
    ax = axs[1, 2]
    for editor, mk in [('ABE', 'o'), ('CBE', '^')]:
        d = CNT[editor]
        vals = []
        for grp in order:
            cands = [c for c in d.columns[2:] if sample_group(c) == grp]
            if not cands:
                vals.append(np.nan)
                continue
            fr = []
            for c in cands:
                v = d[c].astype(float)
                cpm = v / v.sum() * 1e6
                fr.append((cpm < 10).mean() * 100)
            vals.append(np.mean(fr))
        ax.plot(np.arange(len(order)), vals, '-', marker=mk,
                color=EDITOR_COLORS[editor], ms=3, lw=1, mec='black', mew=.3,
                label=editor)
    ax.set_xticks(np.arange(len(order)))
    ax.set_xticklabels(order, fontsize=5, rotation=45, ha='right',
                       rotation_mode='anchor')
    ax.set_ylabel('gRNAs < 10 CPM (%)')
    open_ax(ax)
    ax.legend(fontsize=6, frameon=False, handlelength=1.2)

    for ax, L in zip(axs.ravel(), 'abcdef'):
        panel_letter(ax, L, dx=-0.24, dy=1.05)
    fig.suptitle('EXTENDED 3 — Genetic bottlenecking sets the ceiling on recoverable signal',
                 fontsize=8, fontweight='bold', y=1.0)
    fig.tight_layout(rect=[0, 0, 1, .95])
    fig.savefig(f'{OUT}/ext3.pdf', bbox_inches='tight')
    fig.savefig(f'{OUT}/ext3.png', dpi=150, bbox_inches='tight')
    print('  ext3 done')
    if return_fig:
        return fig
    plt.close(fig)


# ══════════════════════════════════════════════════════════════════════════
# EXTENDED 4 — Editing depth is not a confounder of the hit list
# ══════════════════════════════════════════════════════════════════════════
def ext4(return_fig=False):
    fig, axs = plt.subplots(2, 3, figsize=(183 * MM, 110 * MM))

    # a: |LFC| by editing bin, ABE BM, with control band
    for r, editor in enumerate(['ABE', 'CBE']):
        d = DFS[editor]
        t, c = targ(d), ctrls(d)
        ax = axs[r, 0]
        rows = []
        bins = [0, 20, 40, 60, 80, 100]
        labs = ['0–20', '20–40', '40–60', '60–80', '80–100']
        for s in INVIVO:
            ok = t[f'sensor_reads_{s}'] >= MIN_READS
            tt = t[ok].copy()
            tt['bin'] = pd.cut(tt[f'target_base_edit_perc_{s}'], bins,
                               labels=labs, include_lowest=True)
            for b in labs:
                sub = tt[tt['bin'] == b]
                if len(sub) < 5:
                    continue
                rows.append(dict(tissue=TISSUE_LABEL[s], bin=b,
                                 val=sub[f'LFC_median_{s}'].abs().median(),
                                 n=len(sub)))
        rd = pd.DataFrame(rows)
        for s in INVIVO:
            sub = rd[rd.tissue == TISSUE_LABEL[s]]
            ax.plot(sub['bin'], sub['val'], '-o', ms=2.5, lw=1,
                    color=TISSUE_COLORS[TISSUE_LABEL[s]], mec='black', mew=.3,
                    label=TISSUE_LABEL[s])
        nullv = np.mean([c[f'LFC_median_{s}'].abs().median() for s in INVIVO])
        ax.axhline(nullv, color='black', ls='--', lw=.8)
        ax.text(.02, nullv, ' control median', fontsize=5, va='bottom',
                transform=ax.get_yaxis_transform())
        ax.set_xlabel('Target base editing bin (%)')
        ax.set_ylabel('Median absolute LFC')
        ax.set_title(editor, fontsize=6)
        ax.tick_params(axis='x', labelsize=5, rotation=45)
        open_ax(ax)
        if r == 0:
            ax.legend(fontsize=5, frameon=False, handlelength=1.2)

        # b: fraction significant by bin
        ax = axs[r, 1]
        rows = []
        for s in INVIVO:
            ok = t[f'sensor_reads_{s}'] >= MIN_READS
            tt = t[ok].copy()
            tt['bin'] = pd.cut(tt[f'target_base_edit_perc_{s}'], bins,
                               labels=labs, include_lowest=True)
            for b in labs:
                sub = tt[tt['bin'] == b]
                if len(sub) < 5:
                    continue
                rows.append(dict(tissue=TISSUE_LABEL[s], bin=b,
                                 val=(sub[f'FDR_{s}_stouffers'] < FDR_CUT).mean() * 100))
        rd = pd.DataFrame(rows)
        for s in INVIVO:
            sub = rd[rd.tissue == TISSUE_LABEL[s]]
            ax.plot(sub['bin'], sub['val'], '-o', ms=2.5, lw=1,
                    color=TISSUE_COLORS[TISSUE_LABEL[s]], mec='black', mew=.3)
        ax.set_xlabel('Target base editing bin (%)')
        ax.set_ylabel(f'gRNAs with FDR < {FDR_CUT} (%)')
        ax.set_title(editor, fontsize=6)
        ax.tick_params(axis='x', labelsize=5, rotation=45)
        open_ax(ax)

        # c: sensor read depth vs |LFC| — depth is not driving effect size
        ax = axs[r, 2]
        s = 'bm'
        ok = t[f'sensor_reads_{s}'] > 0
        x = t.loc[ok, f'sensor_reads_{s}'].clip(lower=1)
        y = t.loc[ok, f'LFC_median_{s}'].abs()
        ax.scatter(x, y, s=4, color=blues[3], alpha=.3, linewidths=0)
        rs = scipy.stats.spearmanr(x, y, nan_policy='omit').statistic
        ax.set_xscale('log')
        ax.axvline(MIN_READS, color=greys[4], ls=':', lw=.8)
        ax.text(MIN_READS * 1.15, ax.get_ylim()[1] * .92,
                f'{MIN_READS} reads', fontsize=5, color=greys[4])
        ax.set_xlabel('Sensor reads, bone marrow (log)')
        ax.set_ylabel('Absolute LFC (median)')
        ax.set_title(f'{editor} — ρ = {rs:.2f}', fontsize=6)
        open_ax(ax)

    for ax, L in zip(axs.ravel(), 'abcdef'):
        panel_letter(ax, L, dx=-0.24, dy=1.05)
    fig.suptitle('EXTENDED 4 — Editing depth and read depth do not manufacture the fitness signal',
                 fontsize=8, fontweight='bold', y=1.0)
    fig.tight_layout(rect=[0, 0, 1, .95])
    fig.savefig(f'{OUT}/ext4.pdf', bbox_inches='tight')
    fig.savefig(f'{OUT}/ext4.png', dpi=150, bbox_inches='tight')
    print('  ext4 done')
    if return_fig:
        return fig
    plt.close(fig)


def build_pdf():
    """Render cover + all figures into one multi-page review PDF."""
    from matplotlib.backends.backend_pdf import PdfPages
    import make_cover_page as mcp
    figs = []
    figs.append(mcp.cover_fig())
    for fn in (main_figure, ext1, ext2, ext3, ext4):
        figs.append(fn(return_fig=True))
    path = f'{OUT}/sensor_calibration_proposal.pdf'
    with PdfPages(path) as pdf:
        for f in figs:
            pdf.savefig(f, bbox_inches='tight')
            plt.close(f)
    print('wrote', path)


if __name__ == '__main__':
    build_pdf()
    print('All figures written to', OUT)
