"""Check scientific invariants introduced by the corrected workflow."""
from pathlib import Path
import sys,importlib.util
import numpy as np,torch,yaml
R=Path(__file__).resolve().parents[1];sys.path.insert(0,str(R));import src.models.mc_contrastive_vi as mm
reg=yaml.safe_load((R/'configs/dataset_registry.yaml').read_text());cis=next(d for d in reg['datasets'] if d['accession']=='GSE216146');assert cis['condition_map']=={'control':['PN','PS'],'cisplatin':['CN','CS']};assert set(cis['rescue_values'])=={'PS','CS'};bulk=next(d for d in reg['datasets'] if d['accession']=='GSE286221');assert bulk['status']=='excluded_bulk'
sp=importlib.util.spec_from_file_location('pre',R/'scripts/01_preprocess.py');pre=importlib.util.module_from_spec(sp);sp.loader.exec_module(pre);assert pre._condition_from_dirname('Supp_cisplatin')=='cisplatin';assert pre._condition_from_dirname('Supp_PBS')=='control';assert pre._condition_from_dirname('Supp_DOXACY')=='rescue'
torch.manual_seed(31);m=mm.MCContrastiveVI(mm.MCContrastiveVIConfig(n_genes=20,n_drugs=3,hidden_dim=16));e=m.encode(torch.rand(12,20),torch.tensor([0]*4+[1]*4+[2]*4));assert torch.count_nonzero(e['z_shared_gated'][:4])==0;assert torch.count_nonzero(e['z_drug_cat'][:4])==0
x=torch.randn(11,4,requires_grad=True);y=torch.randn(11,3,requires_grad=True);old=mm.hsic(x,y);gx,gy=torch.autograd.grad(old,(x,y),retain_graph=True);K=mm._rbf_kernel(x);L=mm._rbf_kernel(y);K=K-K.mean(0)-K.mean(1,keepdim=True)+K.mean();L=L-L.mean(0)-L.mean(1,keepdim=True)+L.mean();fast=(K*L).sum()/100;fx,fy=torch.autograd.grad(fast,(x,y));assert torch.allclose(old,fast,atol=1e-7);assert torch.allclose(gx,fx,atol=1e-7);assert torch.allclose(gy,fy,atol=1e-7)
z=np.load(R/'runs/corrected_20260920/inputs.npz');sets=[set(z[k]) for k in ['train','val','test']];assert not sets[0]&sets[1] and not sets[0]&sets[2] and not sets[1]&sets[2];assert len(sets[0])==6000;assert len(z['genes'])==1500
import pandas as pd
o=pd.read_csv(R/'runs/corrected_20260920/cells.csv',index_col=0);assert len(o)==46857;assert not o.rescue.iloc[np.concatenate([z[k] for k in ['train','val','test']])].any();assert (o.loc[o.source_code.isin(['PN','PS']),'drug']=='control').all();assert (o.loc[o.source_code.isin(['CN','CS']),'drug']=='cisplatin').all();q=np.load(R/'runs/corrected_20260920/count_inputs.npz');assert np.allclose(q['C'],q['C'].round())
print('PASS corrected labels, rescue exclusions, bulk exclusion, gating, HSIC values and gradients, integer counts and disjoint splits')
