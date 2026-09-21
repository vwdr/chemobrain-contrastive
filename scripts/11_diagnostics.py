"""Latent diagnostics and independent-library rescue comparisons."""
from pathlib import Path
import numpy as np,pandas as pd,itertools
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score
R=Path(__file__).resolve().parents[1];RUN=R/'runs/corrected_20260920';T=R/'analysis/corrected_20260920';z=np.load(RUN/'inputs.npz');o=pd.read_csv(RUN/'cells.csv',index_col=0);tr=z['train'];te=z['test'];rng=np.random.default_rng(112);rows=[];var=[];samples=[];dep=[];ab=[]
def cka(x,y):
 x=x-x.mean(0);y=y-y.mean(0);return float(np.linalg.norm(x.T@y,'fro')**2/(np.linalg.norm(x.T@x,'fro')*np.linalg.norm(y.T@y,'fro')+1e-20))
for mode in ['full','no_hsic','no_gating','gaussian_vae']:
 for seed in [0,1,2]:
  lat=np.load(RUN/f'{mode}_{seed}_latents.npz');ix=te[z['d'][te]>0];aa=tr[z['b'][tr]==0];bb=te[z['b'][te]==0];clf=LogisticRegression(max_iter=800,class_weight='balanced').fit(lat['bg'][aa],o.source_cell_type.iloc[aa]);ar={'mode':mode,'seed':seed,'background_cell_type_balanced_accuracy':balanced_accuracy_score(o.source_cell_type.iloc[bb],clf.predict(lat['bg'][bb]))}
  for a,b in [('bg','shared'),('bg','drug'),('shared','drug')]:ar[a+'_'+b+'_cka']=cka(lat[a][ix],lat[b][ix])
  ab.append(ar)
  if mode!='full':continue
  for label,rep in [('background',lat['bg']),('shared_ungated',lat['raw_shared']),('drug_ungated',lat['raw_drug'])]:
   c=LogisticRegression(max_iter=800,class_weight='balanced').fit(rep[aa],o.source_cell_type.iloc[aa]);rows.append(dict(seed=seed,representation=label,task='cis_source_cell_type',balanced_accuracy=balanced_accuracy_score(o.source_cell_type.iloc[bb],c.predict(rep[bb]))))
   for study,b in [('doxorubicin',1),('cisplatin',0)]:
    ti=tr[z['b'][tr]==b];vi=te[z['b'][te]==b];y=(z['d']>0).astype(int);c=LogisticRegression(max_iter=800,class_weight='balanced').fit(rep[ti],y[ti]);rows.append(dict(seed=seed,representation=label,task=study+'_treatment',balanced_accuracy=balanced_accuracy_score(y[vi],c.predict(rep[vi]))))
  for a,b in [('bg','shared'),('bg','drug'),('shared','drug')]:
   actual=cka(lat[a][ix],lat[b][ix]);null=[]
   for it in range(100):
    pp=np.arange(len(ix))
    for dd in [1,2]:
     ii=np.where(z['d'][ix]==dd)[0];pp[ii]=rng.permutation(ii)
    null.append(cka(lat[a][ix],lat[b][ix][pp]))
   dep.append(dict(seed=seed,space_a=a,space_b=b,linear_cka=actual,null_mean=np.mean(null),cell_permutation_p=(1+sum(t>=actual for t in null))/101))
  for study in ['GSE216146','GSE271055']:
   ref=tr[(o.study.iloc[tr].values==study)&(z['d'][tr]==0)];dist=np.linalg.norm(lat['raw_shared']-lat['raw_shared'][ref].mean(0),axis=1);oi=o.assign(distance=dist)
   for (sid,arm),sub in oi[oi.study==study].groupby(['sample_id','arm']):samples.append(dict(seed=seed,study=study,sample_id=sid,arm=arm,n_cells=len(sub),mean_distance=sub.distance.mean(),median_distance=sub.distance.median()))
  for drug,k in [('pooled',None),('doxorubicin',1),('cisplatin',2)]:
   ii=te[z['d'][te]>0] if k is None else te[z['d'][te]==k];sh=lat['shared'][ii];dr=lat['drug'][ii]
   def frac(a,b):return float(np.var(a,axis=0).sum()/(np.var(a,axis=0).sum()+np.var(b,axis=0).sum()))
   vals=[]
   for it in range(200):jj=rng.integers(len(ii),size=len(ii));vals.append(frac(sh[jj],dr[jj]))
   var.append(dict(seed=seed,drug=drug,n_test_cells=len(ii),fraction=frac(sh,dr),cell_bootstrap_low=np.quantile(vals,.025),cell_bootstrap_high=np.quantile(vals,.975)))
for name,data in [('latent_probes',rows),('latent_dependence',dep),('variance_robustness',var),('rescue_library_projection',samples),('ablation_diagnostics',ab)]:pd.DataFrame(data).to_csv(T/f'{name}.csv',index=False)
sdf=pd.DataFrame(samples).groupby(['study','sample_id','arm']).mean(numeric_only=True).reset_index();rr=[]
for study,treat,rescue in [('GSE216146','cisplatin','cisplatin_GENUS'),('GSE271055','doxorubicin','doxorubicin_ACY1083')]:
 sub=sdf[(sdf.study==study)&sdf.arm.isin([treat,rescue])].copy();A=sub[sub.arm==rescue].mean_distance.values;B=sub[sub.arm==treat].mean_distance.values;effect=A.mean()-B.mean();row=dict(study=study,contrast=rescue+'_vs_'+treat,mean_difference=effect,n_rescue_libraries=len(A),n_treated_libraries=len(B),batch_stratified_permutation_p=None,n_label_permutations=None)
 if len(A)==3 and len(B)==3:
  sub['batch']=sub.sample_id.str.split('-').str[0];groups=[]
  for batch,g in sub.groupby('batch'):
   vals=g.mean_distance.to_numpy();n=(g.arm==rescue).sum();groups.append([(vals[list(c)].sum(),np.delete(vals,c).sum()) for c in itertools.combinations(range(len(vals)),n)])
  perms=[sum(a for a,b in c)/3-sum(b for a,b in c)/3 for c in itertools.product(*groups)];row['batch_stratified_permutation_p']=sum(abs(v)>=abs(effect)-1e-12 for v in perms)/len(perms);row['n_label_permutations']=len(perms)
 rr.append(row)
pd.DataFrame(rr).to_csv(T/'rescue_exploratory_comparisons.csv',index=False);st=[]
for s,t in itertools.combinations([0,1,2],2):
 A=np.load(RUN/f'full_{s}_latents.npz');B=np.load(RUN/f'full_{t}_latents.npz')
 for key in ['bg','shared','drug']:st.append(dict(seed_a=s,seed_b=t,space=key,test_cka=cka(A[key][te],B[key][te])))
pd.DataFrame(st).to_csv(T/'cross_seed_latent_stability.csv',index=False)
