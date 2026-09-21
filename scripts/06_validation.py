"""Audited benchmark. Held-out cells measure interpolation within deposited libraries."""
import sys,json,copy,argparse,time
from pathlib import Path
import numpy as np,pandas as pd,torch,anndata as ad,scanpy as sc
from sklearn.model_selection import train_test_split
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score,mean_squared_error
R=Path(__file__).resolve().parents[1];sys.path.insert(0,str(R));import src.models.mc_contrastive_vi as mm
O=R/'analysis/corrected_20260920';RUN=R/'runs/corrected_20260920';RUN.mkdir(parents=True,exist_ok=True);torch.set_num_threads(2)
def fast_hsic(x,y):
 if len(x)<4:return x.new_tensor(0.)
 K=mm._rbf_kernel(x);L=mm._rbf_kernel(y);K=K-K.mean(0,keepdim=True)-K.mean(1,keepdim=True)+K.mean();L=L-L.mean(0,keepdim=True)-L.mean(1,keepdim=True)+L.mean();return (K*L).sum()/(len(x)-1)**2
mm.hsic=fast_hsic
def prepare():
 a=ad.read_h5ad(R/'data/processed/recovered_counts.h5ad');obs=a.obs.copy();pure=np.where(~obs.rescue.values)[0];strata=(obs.study.astype(str)+'_'+obs.drug.astype(str)).values
 tr,te=train_test_split(pure,test_size=.2,stratify=strata[pure],random_state=1729);va,te=train_test_split(te,test_size=.5,stratify=strata[te],random_state=1729);rng=np.random.default_rng(1729)
 tr=np.concatenate([rng.choice(tr[strata[tr]==s],min(1500,sum(strata[tr]==s)),replace=False) for s in np.unique(strata[tr])]);v=a[tr].copy();sc.pp.highly_variable_genes(v,flavor='seurat_v3',n_top_genes=1500,batch_key='study');genes=v.var_names[v.var.highly_variable]
 counts=a[:,genes].X.toarray().astype('float32');library=np.asarray(a.X.sum(1)).ravel().astype('float32');np.savez_compressed(RUN/'count_inputs.npz',C=counts,library=library)
 sc.pp.normalize_total(a,target_sum=1e4);sc.pp.log1p(a);X=a[:,genes].X.toarray().astype('float32');d=obs.drug.map({'control':0,'doxorubicin':1,'cisplatin':2}).values.astype('int64');b=(obs.study=='GSE271055').values.astype('int64')
 np.savez_compressed(RUN/'inputs.npz',X=X,d=d,b=b,train=tr,val=va,test=te,genes=np.asarray(genes.tolist(),dtype='U'));obs.to_csv(RUN/'cells.csv');pd.DataFrame({'gene':genes}).to_csv(O/'benchmark_gene_universe.csv',index=False)
 split=np.full(len(obs),'not_used',dtype=object);split[tr]='train';split[va]='validation';split[te]='test';split[obs.rescue.values]='rescue_projection';obs.assign(split=split).groupby(['study','arm','sample_id','split'],observed=True).size().rename('n_cells').reset_index().to_csv(O/'benchmark_split_counts.csv',index=False);print('Prepared',X.shape,len(tr),len(va),len(te),flush=True)
class Ungated(mm.MCContrastiveVI):
 def encode(self,x,drug_idx):
  e=super().encode(x,drug_idx);e['z_shared_gated']=e['z_shared'];e['z_drug_cat']=torch.cat([mm.reparameterize(m,v) for m,v in zip(e['mu_drug_list'],e['lv_drug_list'])],1);return e
def means(m,x,d,mode):
 e=m.encode(x,d);bg=e['mu_bg'];sh=e['mu_shared'];dr=torch.cat(e['mu_drug_list'],1)
 if mode not in ('no_gating','gaussian_vae'):
  sh=sh*(d>0)[:,None];dr=torch.cat([mu*(d==k+1)[:,None] for k,mu in enumerate(e['mu_drug_list'])],1)
 return bg,sh,dr,e
def train(mode,seed,epochs):
 z=np.load(RUN/'inputs.npz');X=torch.from_numpy(z['X']);d=torch.from_numpy(z['d']);b=torch.from_numpy(z['b']);tr=z['train'];va=z['val'];te=z['test'];torch.manual_seed(seed);rng=np.random.default_rng(seed)
 cfg=mm.MCContrastiveVIConfig(n_genes=X.shape[1],n_drugs=3,n_batches=2,hidden_dim=128,dropout=.1);m=(Ungated if mode=='no_gating' else mm.MCContrastiveVI)(cfg);opt=torch.optim.AdamW(m.parameters(),lr=.001,weight_decay=1e-5);best=float('inf');bad=0;history=[];start=time.time()
 for epoch in range(epochs):
  m.train();vals=[]
  for ix in np.array_split(rng.permutation(tr),int(np.ceil(len(tr)/256))):
   x,dd,bb=X[ix],d[ix],b[ix];e=m(x,dd,bb)
   if mode=='gaussian_vae':
    e=m.encode(x,dd);e['z_shared_gated']=e['z_shared'];e['z_drug_cat']=torch.cat([mm.reparameterize(u,v) for u,v in zip(e['mu_drug_list'],e['lv_drug_list'])],1);e['recon']=m.decode(e['z_bg'],e['z_shared_gated'],e['z_drug_cat'],bb)
   if mode in ('no_gating','gaussian_vae'):
    loss=mm.gaussian_nll(x,e['recon'],m.log_sigma).mean()+sum(mm.kl_standard_normal(u,v).mean() for u,v in [(e['mu_bg'],e['lv_bg']),(e['mu_shared'],e['lv_shared'])]+list(zip(e['mu_drug_list'],e['lv_drug_list'])))
    if mode=='no_gating':loss=loss+10*fast_hsic(e['z_bg'],e['z_shared'])+10*fast_hsic(e['z_bg'],e['z_drug_cat'])+5*fast_hsic(e['z_shared'],e['z_drug_cat'])
   else:
    w=0 if mode=='no_hsic' else 1;loss=mm.mc_contrastive_loss(m,e,x,dd,1.,10*w,10*w,5*w)['loss']
   opt.zero_grad();loss.backward();torch.nn.utils.clip_grad_norm_(m.parameters(),5);opt.step();vals.append(float(loss.detach()))
  m.eval();sq=0;n=0
  with torch.no_grad():
   for ix in np.array_split(va,int(np.ceil(len(va)/512))):
    bg,sh,dr,_=means(m,X[ix],d[ix],mode);p=m.decode(bg,sh,dr,b[ix]);sq+=float(((p-X[ix])**2).sum());n+=X[ix].numel()
  vm=sq/n;history.append({'epoch':epoch,'train_objective':np.mean(vals),'val_mse':vm})
  if vm<best-1e-5:best=vm;state=copy.deepcopy(m.state_dict());be=epoch;bad=0
  else:bad+=1
  if epoch%20==0:print(mode,seed,epoch,'val_mse',round(vm,5),'seconds',round(time.time()-start),flush=True)
  if bad>=15:break
 m.load_state_dict(state);m.eval();prefix=RUN/f'{mode}_{seed}';torch.save({'state':state,'config':cfg.__dict__,'seed':seed,'mode':mode,'best_epoch':be},str(prefix)+'.pt');pd.DataFrame(history).to_csv(str(prefix)+'_history.csv',index=False)
 arrays={k:[] for k in ['bg','shared','drug','raw_shared','raw_drug','mse','nll']}
 with torch.no_grad():
  for ix in np.array_split(np.arange(len(X)),int(np.ceil(len(X)/512))):
   bg,sh,dr,e=means(m,X[ix],d[ix],mode);p=m.decode(bg,sh,dr,b[ix]);arrays['bg'].append(bg.numpy());arrays['shared'].append(sh.numpy());arrays['drug'].append(dr.numpy());arrays['raw_shared'].append(e['mu_shared'].numpy());arrays['raw_drug'].append(torch.cat(e['mu_drug_list'],1).numpy());arrays['mse'].append(((p-X[ix])**2).mean(1).numpy());arrays['nll'].append(mm.gaussian_nll(X[ix],p,m.log_sigma).numpy())
 arrays={k:np.concatenate(v) for k,v in arrays.items()};treated=te[z['d'][te]>0];vs=np.var(arrays['shared'][treated],axis=0).sum();vd=np.var(arrays['drug'][treated],axis=0).sum();np.savez_compressed(str(prefix)+'_latents.npz',**arrays)
 result={'mode':mode,'seed':seed,'best_epoch':be,'epochs_run':len(history),'test_mse':float(arrays['mse'][te].mean()),'test_plugin_nll':float(arrays['nll'][te].mean()),'shared_fraction':float(vs/(vs+vd))}
 for name,rep in [('bg',arrays['bg']),('shared_ungated',arrays['raw_shared']),('drug_gated',arrays['drug']),('drug_ungated',arrays['raw_drug'])]:
  for label,y in [('study',z['b']),('condition',z['d'])]:
   clf=LogisticRegression(max_iter=500,class_weight='balanced').fit(rep[tr],y[tr]);result[name+'_'+label+'_balanced_accuracy']=float(balanced_accuracy_score(y[te],clf.predict(rep[te])))
 Path(str(prefix)+'_metrics.json').write_text(json.dumps(result,indent=2));print(result,flush=True)
def pca():
 z=np.load(RUN/'inputs.npz');m=PCA(32,random_state=1729).fit(z['X'][z['train']]);r=m.transform(z['X']);pred=m.inverse_transform(r[z['test']]);row={'test_mse':float(mean_squared_error(z['X'][z['test']],pred))}
 for label,y in [('study',z['b']),('condition',z['d'])]:
  c=LogisticRegression(max_iter=1000,class_weight='balanced').fit(r[z['train']],y[z['train']]);row[label+'_balanced_accuracy']=float(balanced_accuracy_score(y[z['test']],c.predict(r[z['test']])))
 (RUN/'pca_metrics.json').write_text(json.dumps(row,indent=2));np.save(RUN/'pca_latents.npy',r)
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--prepare',action='store_true');ap.add_argument('--mode',default='full');ap.add_argument('--seed',type=int,default=0);ap.add_argument('--epochs',type=int,default=100);ap.add_argument('--pca',action='store_true');a=ap.parse_args()
 if a.prepare:prepare()
 elif a.pca:pca()
 else:train(a.mode,a.seed,a.epochs)
