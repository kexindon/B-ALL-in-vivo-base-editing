import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt
import sys
from calls import *
INPUT100='--input100' in sys.argv
TAG='_input100' if INPUT100 else ''
mpl.rcParams['pdf.fonttype']=42; plt.rcParams['font.family']='sans-serif'
plt.rcParams['font.sans-serif']=['Helvetica','Arial','DejaVu Sans']

# fixed category -> color (color follows category, not rank); original 4 colours kept
CATS=[('Neither','#cccccc'),
      ('Empirical only','#457b9d'),
      ('MAGeCK only','#f4a261'),
      ('MAGeCK-NC only','#2a9d8f'),
      ('Empirical + MAGeCK','#e63946'),
      ('Empirical + MAGeCK-NC','#7b2cbf'),
      ('MAGeCK + MAGeCK-NC','#8d6e63'),
      ('All three','#111111')]
COL=dict(CATS)
def cat(r):
    k=(r.emp_hit,r.default_hit,r.nc_hit)
    return {(0,0,0):'Neither',(1,0,0):'Empirical only',(0,1,0):'MAGeCK only',(0,0,1):'MAGeCK-NC only',
            (1,1,0):'Empirical + MAGeCK',(1,0,1):'Empirical + MAGeCK-NC',(0,1,1):'MAGeCK + MAGeCK-NC',(1,1,1):'All three'}[tuple(map(int,k))]
nl=lambda s:-np.log10(s.clip(lower=1e-300))
ROWS=[('emp_LFC','emp_FDR','Empirical\n-log10(FDR)','Median LFC (vs. Input)'),
      ('default_LFC','default_fdr_dir','MAGeCK (median norm.)\n-log10(one-sided FDR)','MAGeCK LFC'),
      ('nc_LFC','nc_fdr_dir','MAGeCK-NC (control norm.)\n-log10(one-sided FDR)','MAGeCK-NC LFC')]
T=table(INPUT100)
fig,axes=plt.subplots(3,4,figsize=(20,14),sharey='row')
for k,(e,title) in enumerate(zip(EMP,TITLES)):
    d=T[e]; tg=d[d.classification=='targeting guide'].copy(); ct=d[d.classification.isin(CTRL)]
    tg['cat']=tg.apply(cat,axis=1)
    for i,(x,y,yl,xl) in enumerate(ROWS):
        ax=axes[i][k]
        for name,c in CATS:
            s=tg[tg.cat==name]
            if not len(s): continue
            ax.scatter(s[x],nl(s[y]),color=c,edgecolor='none' if name=='Neither' else 'black',linewidth=.3,
                       s=20 if name=='Neither' else 35,alpha=.35 if name=='Neither' else .85,
                       zorder=2 if name=='Neither' else 3+[n for n,_ in CATS].index(name),label=f'{name} (n={len(s)})')
        ax.scatter(ct[x],nl(ct[y]),color='white',edgecolor='black',linewidth=.5,s=20,zorder=20,
                   label=f'Controls (FP {int(ct[["emp_hit","default_hit","nc_hit"][i]].sum())}/{len(ct)})')
        ax.axhline(1,color='red',ls='--',lw=1,alpha=.7)
        for v in (-1,1): ax.axvline(v,color='grey',ls=':',lw=.8)
        ax.set_title(title,fontsize=16); ax.spines[['top','right']].set_visible(False)
        ax.tick_params(labelsize=13); ax.set_xlabel(xl,fontsize=13)
        if k==0: ax.set_ylabel(yl,fontsize=13)
        ax.legend(fontsize=8.5,frameon=False,loc='center' if i==2 else 'upper left')
fig.suptitle('ABE OG (no editing cut-off'+(', median input count ≥100' if INPUT100 else '')+') — Empirical vs MAGeCK vs MAGeCK-NC\n'
             'Hit: Empirical = Fisher FDR<0.1 & |median LFC|>1;  MAGeCK arms = one-sided raw P<0.01 & BH FDR<0.1 & |LFC|>1',fontsize=15,y=1.0)
fig.tight_layout()
for ext in ['png','pdf']: fig.savefig(f'{REPO}mageck/figures/volcano_comparison_ABE_OG_3way{TAG}.{ext}',bbox_inches='tight',dpi=300)
