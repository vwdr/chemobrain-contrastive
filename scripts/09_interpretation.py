"""Integrated-gradient sensitivity and complete mouse gene-set enrichment."""
from pathlib import Path
import sys,itertools,json
import numpy as np,pandas as pd,torch
from scipy.stats import spearmanr,hypergeom
from statsmodels.stats.multitest import multipletests
R=Path(__file__).resolve().parents[1];sys.path.insert(0,str(R));from src.models.mc_contrastive_vi import MCContrastiveVI,MCContrastiveVIConfig
RUN=R/'runs/corrected_20260920';T=R/'analysis/corrected_20260920';z=np.load(RUN/'inputs.npz');X=z['X'];genes=z['genes'];tr=z['train'];te=z['test'];rng=np.random.default_rng(314);torch.set_num_threads(2);rows=[];comp=[];arrays={};nodes,weights=np.polynomial.legendre.leggauss(256)
indices={'shared':np.concatenate([rng.choice(te[z['d'][te]==d],64,replace=False) for d in [1,2]]),'drug0':rng.choice(te[z['d'][te]==1],128,replace=False),'drug1':rng.choice(te[z['d'][te]==2],128,replace=False)}
np.savez(RUN/'attribution_indices.npz',**indices)
for seed in [0,1,2]:
 c=torch.load(RUN/f'full_{seed}.pt',map_location='cpu',weights_only=True);m=MCContrastiveVI(MCContrastiveVIConfig(**c['config']));m.load_state_dict(c['state']);m.eval()
 for axis,ix in indices.items():
  for base in ['zero','control_median']:
   controls=tr[z['d'][tr]==0]
   if axis!='shared':controls=controls[z['b'][controls]==(1 if axis=='drug0' else 0)]
   baseline=np.zeros(X.shape[1],dtype=np.float32) if base=='zero' else np.median(X[controls],axis=0)
   for target in ['coordinate_sum','norm_squared']:
    def fun(x):
     if axis=='shared':mu=m.shared_head(m.shared_trunk(x))[0]
     else:k=int(axis[-1]);mu=m.drug_heads[k](m.drug_trunks[k](x))[0]
     return mu.sum(1) if target=='coordinate_sum' else mu.square().sum(1)
    inp=torch.tensor(X[ix]);ref=torch.tensor(baseline)[None,:];delta=inp-ref;acc=torch.zeros_like(inp)
    with torch.no_grad():diff=(fun(inp)-fun(ref)).numpy()
    for order in [256,1024,4096]:
     nodes,weights=np.polynomial.legendre.leggauss(order);acc=torch.zeros_like(inp)
     for node,weight in zip(nodes,weights):
      xx=(ref+delta*float((node+1)/2)).detach().requires_grad_(True);grad=torch.autograd.grad(fun(xx).sum(),xx)[0];acc+=grad*float(weight/2)
     ig=(delta*acc).detach().numpy();err=np.abs(ig.sum(1)-diff)
     if np.median(err)/(np.mean(np.abs(diff))+1e-8)<.01:break
    comp.append({'seed':seed,'axis':axis,'baseline':base,'target':target,'median_absolute_error':np.median(err),'max_absolute_error':err.max(),'median_relative_error':np.median(err/(np.abs(diff)+1e-8)),'median_scaled_error':np.median(err)/(np.mean(np.abs(diff))+1e-8),'integration_rule':'adaptive_Gauss_Legendre','quadrature_points':order});arrays[(seed,axis,base,target)]=ig;np.savez_compressed(RUN/f'ig_{seed}_{axis}_{base}_{target}.npz',attributions=ig,indices=ix,baseline=baseline,output_difference=diff)
    for g,a,b in zip(genes,ig.mean(0),np.abs(ig).mean(0)):rows.append({'seed':seed,'axis':axis,'baseline':base,'target':target,'gene':g,'mean_ig':a,'mean_abs_ig':b})
 print('Attributions',seed,flush=True)
pd.DataFrame(rows).to_csv(T/'complete_corrected_attributions.csv',index=False);pd.DataFrame(comp).to_csv(T/'ig_completeness.csv',index=False)
def sim(a,b):
 A=set(np.argsort(a)[-50:]);B=set(np.argsort(b)[-50:]);return float(spearmanr(a,b).statistic),len(A&B)/len(A|B)
s=[]
for axis,base,target in itertools.product(indices,['zero','control_median'],['coordinate_sum','norm_squared']):
 for a,b in itertools.combinations([0,1,2],2):
  r,j=sim(np.abs(arrays[(a,axis,base,target)]).mean(0),np.abs(arrays[(b,axis,base,target)]).mean(0));s.append(dict(comparison='seed',axis=axis,baseline=base,target=target,a=a,b=b,spearman=r,top50_jaccard=j))
for seed,axis,target in itertools.product([0,1,2],indices,['coordinate_sum','norm_squared']):
 r,j=sim(np.abs(arrays[(seed,axis,'zero',target)]).mean(0),np.abs(arrays[(seed,axis,'control_median',target)]).mean(0));s.append(dict(comparison='baseline',axis=axis,baseline='paired',target=target,a=seed,b=seed,spearman=r,top50_jaccard=j))
for axis in indices:
 ar=np.abs(arrays[(0,axis,'zero','norm_squared')]);ref=ar.mean(0)
 for b in range(100):
  r,j=sim(ref,ar[rng.integers(len(ar),size=len(ar))].mean(0));s.append(dict(comparison='cell_bootstrap',axis=axis,baseline='zero',target='norm_squared',a=0,b=b,spearman=r,top50_jaccard=j))
pd.DataFrame(s).to_csv(T/'attribution_stability.csv',index=False)
u=set(genes);en=[];df=pd.DataFrame(rows);df=df[(df.baseline=='zero')&(df.target=='norm_squared')]
for axis in indices:
 top=set(df[df.axis==axis].groupby('gene').mean_abs_ig.mean().nlargest(50).index)
 for lib,file in [('Reactome','m2.cp.reactome.v2025.1.Mm.symbols.gmt'),('GO_BP','m5.go.bp.v2025.1.Mm.symbols.gmt')]:
  subset=[]
  for line in (R/'data/evidence'/file).read_text().splitlines():
   a=line.split('\t');g=set(a[2:])&u
   if not 10<=len(g)<=500:continue
   overlap=top&g;k=len(overlap);subset.append(dict(axis=axis,library=lib,term=a[0],overlap=k,set_size_in_universe=len(g),list_size=50,universe_size=len(u),fold_enrichment=k/(50*len(g)/len(u)),p_value=hypergeom.sf(k-1,len(u),len(g),50),genes='|'.join(sorted(overlap))))
  qq=multipletests([d['p_value'] for d in subset],method='fdr_bh')[1]
  for d,q in zip(subset,qq):d['q_value']=q
  en.extend(subset)
pd.DataFrame(en).to_csv(T/'pathway_enrichment.csv',index=False)
