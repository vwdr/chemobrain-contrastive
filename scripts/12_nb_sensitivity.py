"""Negative-binomial count likelihood sensitivity with matched benchmark splits."""
from pathlib import Path
import sys,json,copy,argparse,importlib.util
import numpy as np,pandas as pd,torch
import torch.nn.functional as F
R=Path(__file__).resolve().parents[1];sys.path.insert(0,str(R));from src.models.mc_contrastive_vi import MCContrastiveVI,MCContrastiveVIConfig,kl_standard_normal
spec=importlib.util.spec_from_file_location('v',R/'scripts/06_validation.py');v=importlib.util.module_from_spec(spec);spec.loader.exec_module(v);RUN=R/'runs/corrected_20260920';torch.set_num_threads(1)
def train(seed):
 z=np.load(RUN/'inputs.npz');q=np.load(RUN/'count_inputs.npz');X=torch.tensor(z['X']);C=torch.tensor(q['C']);lib=torch.tensor(q['library']);d=torch.tensor(z['d']);b=torch.tensor(z['b']);tr=z['train'];va=z['val'];te=z['test'];assert torch.allclose(C,C.round())
 torch.manual_seed(seed);rng=np.random.default_rng(seed);cfg=MCContrastiveVIConfig(n_genes=X.shape[1],n_drugs=3,n_batches=2,hidden_dim=128);m=MCContrastiveVI(cfg);theta=torch.nn.Parameter(torch.zeros(X.shape[1]));opt=torch.optim.AdamW(list(m.parameters())+[theta],lr=.001,weight_decay=1e-5);best=float('inf');hist=[];bad=0
 def nll(c,mu):
  t=F.softplus(theta)+1e-4;mu=mu.clamp_min(1e-6)
  return -(torch.lgamma(c+t)-torch.lgamma(t)-torch.lgamma(c+1)+t*(torch.log(t)-torch.log(t+mu))+c*(torch.log(mu)-torch.log(t+mu))).sum(1)
 for epoch in range(100):
  m.train()
  for ix in np.array_split(rng.permutation(tr),int(np.ceil(len(tr)/256))):
   e=m(X[ix],d[ix],b[ix]);mu=F.softmax(e['recon'],dim=1)*C[ix].sum(1,keepdim=True);loss=nll(C[ix],mu).mean()+kl_standard_normal(e['mu_bg'],e['lv_bg']).mean();mask=d[ix]>0
   if mask.any():loss=loss+kl_standard_normal(e['mu_shared'][mask],e['lv_shared'][mask]).mean()
   for k,(u,w) in enumerate(zip(e['mu_drug_list'],e['lv_drug_list'])):
    mask=d[ix]==k+1
    if mask.any():loss=loss+kl_standard_normal(u[mask],w[mask]).mean()
   loss=loss+10*v.fast_hsic(e['z_bg'],e['z_shared'])+10*v.fast_hsic(e['z_bg'],e['z_drug_cat'])+5*v.fast_hsic(e['z_shared'],e['z_drug_cat']);opt.zero_grad();loss.backward();torch.nn.utils.clip_grad_norm_(list(m.parameters())+[theta],5);opt.step()
  m.eval();vals=[]
  with torch.no_grad():
   for ix in np.array_split(va,10):
    bg,sh,dr,_=v.means(m,X[ix],d[ix],'full');mu=F.softmax(m.decode(bg,sh,dr,b[ix]),1)*C[ix].sum(1,keepdim=True);pred=torch.log1p(mu/lib[ix,None]*1e4);vals.extend(((pred-X[ix])**2).mean(1).tolist())
  vm=float(np.mean(vals));hist.append({'epoch':epoch,'val_mse':vm})
  if vm<best-1e-5:best=vm;state=copy.deepcopy(m.state_dict());disp=theta.detach().clone();be=epoch;bad=0
  else:bad+=1
  if epoch%20==0:print('NB',seed,epoch,vm,flush=True)
  if bad>=15:break
 m.load_state_dict(state);theta.data.copy_(disp);m.eval();vals=[];ls=[];shs=[];drs=[]
 with torch.no_grad():
  for ix in np.array_split(te,10):
   bg,sh,dr,_=v.means(m,X[ix],d[ix],'full');mu=F.softmax(m.decode(bg,sh,dr,b[ix]),1)*C[ix].sum(1,keepdim=True);pred=torch.log1p(mu/lib[ix,None]*1e4);vals.extend(((pred-X[ix])**2).mean(1).tolist());ls.extend(nll(C[ix],mu).tolist());shs.append(sh.numpy());drs.append(dr.numpy())
 sh=np.concatenate(shs)[z['d'][te]>0];dr=np.concatenate(drs)[z['d'][te]>0];vs=np.var(sh,axis=0).sum();vd=np.var(dr,axis=0).sum();row={'mode':'negative_binomial','seed':seed,'best_epoch':be,'epochs_run':len(hist),'test_mse':float(np.mean(vals)),'test_plugin_count_nll':float(np.mean(ls)),'shared_fraction':float(vs/(vs+vd))};(RUN/f'negative_binomial_{seed}_metrics.json').write_text(json.dumps(row,indent=2));pd.DataFrame(hist).to_csv(RUN/f'negative_binomial_{seed}_history.csv',index=False);torch.save({'state':state,'theta_unconstrained':disp,'config':cfg.__dict__,'seed':seed,'best_epoch':be},RUN/f'negative_binomial_{seed}.pt');print(row,flush=True)
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--seed',type=int,default=0);a=ap.parse_args();train(a.seed)
