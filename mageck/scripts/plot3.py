"""3-row volcano: Empirical / MAGeCK (median norm) / MAGeCK-NC (control norm), ABE OG.
Each row is coloured by its own hit calls: enriched red, depleted blue, rest grey; controls white.

    python plot3.py              # all guides
    python plot3.py --input100   # median raw input count >= 100
"""
import sys
import numpy as np, matplotlib as mpl, matplotlib.pyplot as plt, seaborn as sns
from calls import *
INPUT100='--input100' in sys.argv
TAG='_input100' if INPUT100 else ''

mpl.rcParams['pdf.fonttype']=42; mpl.rcParams['ps.fonttype']=42
plt.rcParams['font.family']='sans-serif'; plt.rcParams['font.sans-serif']=['Helvetica','Arial','DejaVu Sans']
reds=sns.color_palette('Reds').as_hex(); blues=sns.color_palette('Blues').as_hex()
# same red/blue as the original Ext Fig 4 comparison panels
DIR={'Enriched':'#e63946','Depleted':'#457b9d','Not significant':'lightgrey'}

nl=lambda s:-np.log10(s.clip(lower=1e-300))
ROWS=[('emp','emp_LFC','emp_FDR','emp_hit','Empirical\n-log10(FDR)','Median LFC (vs. Input)'),
      ('mag','default_LFC','default_fdr_dir','default_hit','MAGeCK\n-log10(FDR)','MAGeCK LFC'),
      ('nc','nc_LFC','nc_fdr_dir','nc_hit','MAGeCK-NC\n-log10(FDR)','MAGeCK-NC LFC')]

T=table(INPUT100)
fig,axes=plt.subplots(3,4,figsize=(20,14),sharey='row')
for k,(e,title) in enumerate(zip(EMP,TITLES)):
    d=T[e]; tg=d[d.classification=='targeting guide']; ct=d[d.classification.isin(CTRL)]
    for i,(_,x,y,hit,yl,xl) in enumerate(ROWS):
        ax=axes[i][k]
        groups={'Not significant':tg[~tg[hit]],'Enriched':tg[tg[hit]&(tg[x]>0)],'Depleted':tg[tg[hit]&(tg[x]<0)]}
        for name,s in groups.items():
            ns=name=='Not significant'
            ax.scatter(s[x],nl(s[y]),color=DIR[name],edgecolor='none' if ns else 'black',linewidth=.3,
                       s=20 if ns else 50,alpha=.5,zorder=2 if ns else 4,label=f'{name} (n={len(s)})')
        ax.scatter(ct[x],nl(ct[y]),color='white',edgecolor='black',linewidth=.5,s=20,zorder=5,
                   label=f'Controls ({int(ct[hit].sum())}/{len(ct)} called)')
        ax.axhline(1,color='grey',ls='--',lw=.8)
        for v in (-1,1): ax.axvline(v,color='grey',ls='--',lw=.8)
        ax.set_title(title,fontsize=16); ax.spines[['top','right']].set_visible(False)
        ax.tick_params(axis='both',which='major',labelsize=14,direction='out')
        ax.set_xlabel(xl,fontsize=14)
        if k==0: ax.set_ylabel(yl,fontsize=14)
        ax.legend(fontsize=10,frameon=False,loc='center' if i==2 else 'upper left')
sub='median input count ≥100, ' if INPUT100 else ''
fig.suptitle(f'ABE OG ({sub}no editing cut-off)',fontsize=16,y=1.0)
fig.tight_layout()
for ext in ['pdf','png']: fig.savefig(f'{REPO}mageck/figures/volcano_comparison_ABE_OG_3way{TAG}.{ext}',bbox_inches='tight',dpi=300)
