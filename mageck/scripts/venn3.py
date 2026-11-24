"""3-way Venn of targeting-guide hits per tissue: Empirical / MAGeCK / MAGeCK-NC, ABE OG.

    python venn3.py              # all guides
    python venn3.py --input100   # median raw input count >= 100
"""
import sys
import matplotlib as mpl, matplotlib.pyplot as plt, seaborn as sns
from matplotlib_venn import venn3
from calls import *
INPUT100='--input100' in sys.argv
TAG='_input100' if INPUT100 else ''

mpl.rcParams['pdf.fonttype']=42; mpl.rcParams['ps.fonttype']=42
plt.rcParams['font.family']='sans-serif'; plt.rcParams['font.sans-serif']=['Helvetica','Arial','DejaVu Sans']
COLORS=('#457b9d','#f4a261',sns.color_palette('Greens').as_hex()[3])   # Empirical, MAGeCK, MAGeCK-NC

T=table(INPUT100)
fig,axes=plt.subplots(2,2,figsize=(9,9))
for ax,(e,title) in zip(axes.flat,zip(EMP,TITLES)):
    tg=T[e][T[e].classification=='targeting guide']
    sets=[set(tg.gRNA_id[tg[h]]) for h in ('emp_hit','default_hit','nc_hit')]
    names=('Empirical','MAGeCK','MAGeCK-NC')
    v=venn3(sets,set_labels=[f'{n} ({len(x)})' for n,x in zip(names,sets)],set_colors=COLORS,alpha=.55,ax=ax)
    for t,x in zip(v.set_labels,sets):
        if t: t.set_fontsize(12); t.set_visible(len(x)>0)   # an empty set has no circle to label
    empty=[n for n,x in zip(names,sets) if not x]
    if empty: ax.text(.5,-.04,', '.join(empty)+': no hits',transform=ax.transAxes,ha='center',va='top',fontsize=11,color='dimgrey')
    for t in v.subset_labels:
        if t: t.set_fontsize(12)
    ax.set_title(f'{title}\nABE',fontsize=14)
sub='median input count ≥100, ' if INPUT100 else ''
fig.suptitle(f'ABE OG ({sub}no editing cut-off)',fontsize=14,y=1.0)
fig.tight_layout(h_pad=4)
for ext in ['pdf','png']: fig.savefig(f'{REPO}mageck/figures/venn_comparison_ABE_OG_3way{TAG}.{ext}',bbox_inches='tight',dpi=300)
