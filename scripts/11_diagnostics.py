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
samp=pd.DataFrame(samples)
sdf=samp.groupby(['study','sample_id','arm']).mean(numeric_only=True).reset_index()
rr=[];pair_rows=[];perm_rows=[];seed_rows=[]
# GSE216146 GEO replicate pairs: cisplatin nonstim -> cisplatin + GENUS.
pairs=[
 ('replicate_1','D20-6409','D20-6410'),
 ('replicate_2','D21-2750','D21-2752'),
 ('replicate_3','D21-2751','D21-2753'),
]
cis=sdf[(sdf.study=='GSE216146')&sdf.arm.isin(['cisplatin','cisplatin_GENUS'])].copy()
lookup=dict(zip(cis.sample_id,cis.mean_distance))
diffs=[]
for rep,treat_sid,rescue_sid in pairs:
 assert treat_sid in lookup and rescue_sid in lookup,(rep,treat_sid,rescue_sid)
 d=float(lookup[rescue_sid]-lookup[treat_sid]);diffs.append(d)
 pair_rows.append(dict(replicate=rep,treated_library=treat_sid,rescue_library=rescue_sid,treated_mean_distance=float(lookup[treat_sid]),rescue_mean_distance=float(lookup[rescue_sid]),rescue_minus_treatment=d))
effect=float(np.mean(diffs));vals=[]
for signs in itertools.product([-1,1],repeat=3):
 stat=float(np.mean(np.asarray(diffs)*np.asarray(signs)));vals.append(stat)
 perm_rows.append(dict(sign_replicate_1=signs[0],sign_replicate_2=signs[1],sign_replicate_3=signs[2],mean_difference=stat,absolute_ge_observed=abs(stat)>=abs(effect)-1e-12))
p=float(sum(abs(v)>=abs(effect)-1e-12 for v in vals)/len(vals))
rr.append(dict(study='GSE216146',contrast='cisplatin_GENUS_vs_cisplatin',mean_difference=effect,n_rescue_libraries=3,n_treated_libraries=3,batch_stratified_permutation_p=p,n_label_permutations=len(vals),reference_design='paired_GEO_replicate_sign_flip',exact_two_sided_p=p,n_sign_assignments=len(vals),minimum_attainable_two_sided_p=2/len(vals),minimum_attainable_one_sided_p=1/len(vals)))
for seed in [0,1,2]:
 q=samp[(samp.seed==seed)&(samp.study=='GSE216146')]
 lk=dict(zip(q.sample_id,q.mean_distance));dd=np.asarray([float(lk[r]-lk[t]) for _,t,r in pairs]);ef=float(dd.mean());pv=[]
 for signs in itertools.product([-1,1],repeat=3):pv.append(float(np.mean(dd*np.asarray(signs))))
 seed_rows.append(dict(seed=seed,mean_difference=ef,exact_two_sided_p=sum(abs(v)>=abs(ef)-1e-12 for v in pv)/len(pv),replicate_1_difference=dd[0],replicate_2_difference=dd[1],replicate_3_difference=dd[2]))
# GSE271055 contains one pooled treated library and one pooled rescue library; descriptive only.
dox=sdf[(sdf.study=='GSE271055')&sdf.arm.isin(['doxorubicin','doxorubicin_ACY1083'])].copy();A=dox[dox.arm=='doxorubicin_ACY1083'].mean_distance.values;B=dox[dox.arm=='doxorubicin'].mean_distance.values
rr.append(dict(study='GSE271055',contrast='doxorubicin_ACY1083_vs_doxorubicin',mean_difference=float(A.mean()-B.mean()),n_rescue_libraries=len(A),n_treated_libraries=len(B),batch_stratified_permutation_p=None,n_label_permutations=None,reference_design='descriptive_one_pooled_library_per_arm',exact_two_sided_p=None,n_sign_assignments=None,minimum_attainable_two_sided_p=None,minimum_attainable_one_sided_p=None))
pd.DataFrame(rr).to_csv(T/'rescue_exploratory_comparisons.csv',index=False)
pd.DataFrame(pair_rows).to_csv(T/'rescue_replicate_pairs.csv',index=False)
pd.DataFrame(perm_rows).to_csv(T/'rescue_signflip_distribution.csv',index=False)
pd.DataFrame(seed_rows).to_csv(T/'rescue_seedwise_signflip.csv',index=False)
st=[]
for s,t in itertools.combinations([0,1,2],2):
 A=np.load(RUN/f'full_{s}_latents.npz');B=np.load(RUN/f'full_{t}_latents.npz')
 for key in ['bg','shared','drug']:st.append(dict(seed_a=s,seed_b=t,space=key,test_cka=cka(A[key][te],B[key][te])))
pd.DataFrame(st).to_csv(T/'cross_seed_latent_stability.csv',index=False)