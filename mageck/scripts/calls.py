"""Hit calls for ABE OG: Empirical vs MAGeCK (median norm) vs MAGeCK-NC (control norm).

    python calls.py              # all 1,070 guides
    python calls.py --input100   # median raw input count >= 100 (paper's coverage filter)
"""
import os, sys
import pandas as pd, numpy as np
from statsmodels.stats.multitest import multipletests

REPO=os.path.abspath(os.path.join(os.path.dirname(__file__),'..','..'))+'/'
MAG_DIR=REPO+'mageck/og/ABE'
EMP=['d15','spleen','bm','men']; MAG=['d15','spleen','bonemarrow','meninges']; TITLES=['D15','Spleen','Bone Marrow','Meninges']
FDR, LFC, P = 0.1, 1.0, 0.01
CTRL=['non-targeting control','safe-targeting control']
INPUTS=['input_rep1','input_rep2','input_rep3']
MIN_INPUT=100
bh=lambda p: multipletests(p,method='fdr_bh')[1]

def load_mag(mode, cond, input100=False):
    sub=f'input100/{mode}' if input100 else mode
    m=pd.read_csv(f'{MAG_DIR}/{sub}/ABE_{cond}_input.sgrna_summary.txt',sep='\t').rename(columns={'sgrna':'gRNA_id'})
    # one-sided tests, BH within each tail (paper: raw P<0.01, BH FDR<0.1, |LFC|>1)
    m['fdr_high']=bh(m['p.high']); m['fdr_low']=bh(m['p.low'])
    enr=(m.LFC> LFC)&(m['p.high']<P)&(m.fdr_high<FDR)
    dep=(m.LFC<-LFC)&(m['p.low'] <P)&(m.fdr_low <FDR)
    m['hit']=enr|dep
    # one-sided FDR in the direction of the LFC, for the y-axis
    m['fdr_dir']=np.where(m.LFC>=0,m.fdr_high,m.fdr_low)
    return m[['gRNA_id','LFC','fdr_dir','FDR','hit']]

def input100_ids():
    c=pd.read_csv(REPO+'counts/og/OG_ABE_COUNTS.txt',sep='\t')
    return set(c.loc[c[INPUTS].median(1)>=MIN_INPUT,'sgRNA'])

def table(input100=False):
    emp=pd.read_csv(REPO+'LFC-FDR/ABE_OG_LFC_FDR.csv',low_memory=False)
    emp=emp.loc[:,~emp.columns.duplicated()]
    if input100:
        # drop low-coverage guides, then redo BH on the published Fisher p-values over the kept set
        emp=emp[emp.gRNA_id.isin(input100_ids())].copy()
        for e in EMP:
            emp[f'FDR_{e}_fishers']=np.minimum(bh(emp[f'p_high_unadjusted_{e}_fishers']),
                                               bh(emp[f'p_low_unadjusted_{e}_fishers']))
    out={}
    for e,m in zip(EMP,MAG):
        d=emp[['gRNA_id','classification',f'LFC_median_{e}',f'FDR_{e}_fishers']].rename(columns={f'LFC_median_{e}':'emp_LFC',f'FDR_{e}_fishers':'emp_FDR'})
        d['emp_hit']=(d.emp_FDR<FDR)&(d.emp_LFC.abs()>LFC)
        for mode in ['default','nc']:
            x=load_mag(mode,m,input100).add_prefix(f'{mode}_').rename(columns={f'{mode}_gRNA_id':'gRNA_id'})
            d=d.merge(x,on='gRNA_id',how='inner')
        out[e]=d
    return out

if __name__=='__main__':
    input100='--input100' in sys.argv
    rows=[]
    for (e,d),t in zip(table(input100).items(),TITLES):
        tg=d[d.classification=='targeting guide']; c=d[d.classification.isin(CTRL)]
        r={'tissue':t,'n targeting':len(tg)}
        for k,lab in [('emp_hit','Empirical'),('default_hit','MAGeCK'),('nc_hit','MAGeCK-NC')]:
            r[f'{lab} hits']=int(tg[k].sum()); r[f'{lab} ctrl FP']=f'{int(c[k].sum())}/{len(c)}'
        r['E&MAG']=int((tg.emp_hit&tg.default_hit).sum()); r['E&NC']=int((tg.emp_hit&tg.nc_hit).sum())
        r['enr/dep E']=f"{int((tg.emp_hit&(tg.emp_LFC>0)).sum())}/{int((tg.emp_hit&(tg.emp_LFC<0)).sum())}"
        r['enr/dep MAG']=f"{int((tg.default_hit&(tg.default_LFC>0)).sum())}/{int((tg.default_hit&(tg.default_LFC<0)).sum())}"
        r['enr/dep NC']=f"{int((tg.nc_hit&(tg.nc_LFC>0)).sum())}/{int((tg.nc_hit&(tg.nc_LFC<0)).sum())}"
        rows.append(r)
    pd.set_option('display.width',250); print(pd.DataFrame(rows).set_index('tissue').T)
