"""Exploratory library-aware Welch comparisons of log pseudobulk CPM."""
from pathlib import Path
import numpy as np,pandas as pd,anndata as ad
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
R=Path(__file__).resolve().parents[1];T=R/'analysis/corrected_20260920';a=ad.read_h5ad(R/'data/processed/recovered_counts.h5ad');a=a[a.obs.study=='GSE216146'].copy();groups={};meta=[]
for (typ,sid,arm),ix in a.obs.groupby(['source_cell_type','sample_id','arm'],observed=True).indices.items():
 meta.append({'cell_type':typ,'sample_id':sid,'arm':arm,'n_cells':len(ix)})
 if len(ix)>=20:groups[(typ,sid,arm)]=np.asarray(a.X[ix].sum(0)).ravel()
rows=[];summ=[]
for typ in a.obs.source_cell_type.unique():
 for A,B in [('cisplatin','control'),('cisplatin_GENUS','cisplatin'),('control_GENUS','control')]:
  aa=[v for (t,s,c),v in groups.items() if t==typ and c==A];bb=[v for (t,s,c),v in groups.items() if t==typ and c==B]
  if len(aa)!=3 or len(bb)!=3:continue
  raw=np.array(aa+bb);keep=(raw>=10).sum(0)>=3;v=np.log2(raw/raw.sum(1)[:,None]*1e6+1)[:,keep];stat,p=ttest_ind(v[:3],v[3:],axis=0,equal_var=False);p=np.nan_to_num(p,nan=1);q=multipletests(p,method='fdr_bh')[1];effect=v[:3].mean(0)-v[3:].mean(0);genes=a.var_names[keep];contrast=A+'_vs_'+B;summ.append({'cell_type':typ,'contrast':contrast,'n_genes':len(genes),'n_fdr_005':sum(q<.05)})
  for g,e,s,pv,qv in zip(genes,effect,stat,p,q):rows.append({'cell_type':typ,'contrast':contrast,'gene':g,'mean_log2CPM_difference':e,'welch_t':s,'p_value':pv,'q_value':qv,'n_a':3,'n_b':3})
pd.DataFrame(rows).to_csv(T/'pseudobulk_exploratory.csv.gz',index=False);pd.DataFrame(meta).to_csv(T/'pseudobulk_library_qc.csv',index=False);pd.DataFrame(summ).to_csv(T/'pseudobulk_summary.csv',index=False)
