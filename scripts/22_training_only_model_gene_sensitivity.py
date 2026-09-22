"""Training-only model-gene sensitivity for the corrected benchmark.

The canonical recovered matrices applied a minimum gene-detection filter before
the cell-level train/validation/test split. This sensitivity analysis preserves
the recovered cells and common-gene matrix but removes genes that are detected
in fewer than ten canonical training cells within either study before selecting
the 1,500 model genes. It then refits the full MC-ContrastiveVI for seeds 0-2
and PCA(32).

This does not reconstruct genes removed upstream from the deposited
source-filtered matrices. It tests whether the manuscript's model-level results
depend materially on using the pre-split gene-detection filter to define the
candidate gene universe.
"""
import sys,json,copy,time
from pathlib import Path
import numpy as np,pandas as pd,torch,anndata as ad,scanpy as sc
from sklearn.model_selection import train_test_split
from sklearn.decomposition import PCA
from sklearn.metrics import mean_squared_error

R=Path(__file__).resolve().parents[1];sys.path.insert(0,str(R))
import src.models.mc_contrastive_vi as mm
O=R/"analysis"/"corrected_20260920";RUN=R/"runs"/"corrected_20260920"
RUN.mkdir(parents=True,exist_ok=True);torch.set_num_threads(2)

def fast_hsic(x,y):
    if len(x)<4:return x.new_tensor(0.)
    K=mm._rbf_kernel(x);L=mm._rbf_kernel(y)
    K=K-K.mean(0,keepdim=True)-K.mean(1,keepdim=True)+K.mean()
    L=L-L.mean(0,keepdim=True)-L.mean(1,keepdim=True)+L.mean()
    return (K*L).sum()/(len(x)-1)**2
mm.hsic=fast_hsic

def split(obs):
    pure=np.where(~obs.rescue.astype(bool).to_numpy())[0]
    strata=(obs.study.astype(str)+"_"+obs.drug.astype(str)).to_numpy()
    tr,te=train_test_split(pure,test_size=.2,stratify=strata[pure],random_state=1729)
    va,te=train_test_split(te,test_size=.5,stratify=strata[te],random_state=1729)
    rng=np.random.default_rng(1729)
    tr=np.concatenate([rng.choice(tr[strata[tr]==s],min(1500,int(np.sum(strata[tr]==s))),replace=False) for s in np.unique(strata[tr])])
    return tr.astype(int),va.astype(int),te.astype(int)

def means(m,x,d):
    e=m.encode(x,d);bg=e["mu_bg"];sh=e["mu_shared"]*(d>0)[:,None]
    dr=torch.cat([mu*(d==k+1)[:,None] for k,mu in enumerate(e["mu_drug_list"])],1)
    return bg,sh,dr,e

a=ad.read_h5ad(R/"data"/"processed"/"recovered_counts.h5ad");obs=a.obs.copy()
tr,va,te=split(obs)
keep=np.ones(a.n_vars,dtype=bool);detail=[]
for study in sorted(obs.study.astype(str).unique()):
    ix=tr[obs.iloc[tr].study.astype(str).to_numpy()==study]
    det=np.asarray((a[ix].X>0).sum(0)).ravel()
    detail.append(pd.DataFrame({"study":study,"gene":a.var_names.astype(str),"training_cells_detected":det.astype(int),"passes_10":det>=10}))
    keep &= det>=10
pd.concat(detail,ignore_index=True).to_csv(O/"training_only_model_gene_filter_by_study.csv",index=False)
candidate=a.var_names[keep]
v=a[tr,candidate].copy();sc.pp.highly_variable_genes(v,flavor="seurat_v3",n_top_genes=1500,batch_key="study")
genes=v.var_names[v.var.highly_variable].astype(str).tolist()
assert len(genes)==1500
pd.DataFrame({"gene":genes}).to_csv(O/"training_only_model_gene_universe.csv",index=False)
summary={"recovered_common_genes":int(a.n_vars),"training_only_candidate_genes":int(keep.sum()),"removed_candidate_genes":int((~keep).sum()),"model_genes":len(genes)}
(O/"training_only_model_gene_design.json").write_text(json.dumps(summary,indent=2)+"\n")

norm=a.copy();sc.pp.normalize_total(norm,target_sum=1e4);sc.pp.log1p(norm)
X_np=norm[:,genes].X.toarray().astype("float32");X=torch.from_numpy(X_np)
d=torch.from_numpy(obs.drug.map({"control":0,"doxorubicin":1,"cisplatin":2}).to_numpy(dtype="int64"))
b=torch.from_numpy((obs.study=="GSE271055").to_numpy(dtype="int64"))
rows=[]
for seed in [0,1,2]:
    torch.manual_seed(seed);rng=np.random.default_rng(seed)
    cfg=mm.MCContrastiveVIConfig(n_genes=1500,n_drugs=3,n_batches=2,hidden_dim=128,dropout=.1)
    m=mm.MCContrastiveVI(cfg);opt=torch.optim.AdamW(m.parameters(),lr=.001,weight_decay=1e-5)
    best=float("inf");bad=0;state=None;be=-1
    for epoch in range(100):
        m.train()
        for ix in np.array_split(rng.permutation(tr),int(np.ceil(len(tr)/256))):
            e=m(X[ix],d[ix],b[ix]);loss=mm.mc_contrastive_loss(m,e,X[ix],d[ix],1.,10.,10.,5.)["loss"]
            opt.zero_grad();loss.backward();torch.nn.utils.clip_grad_norm_(m.parameters(),5);opt.step()
        m.eval();sq=0;n=0
        with torch.no_grad():
            for ix in np.array_split(va,int(np.ceil(len(va)/512))):
                bg,sh,dr,_=means(m,X[ix],d[ix]);p=m.decode(bg,sh,dr,b[ix]);sq+=float(((p-X[ix])**2).sum());n+=X[ix].numel()
        vm=sq/n
        if vm<best-1e-5:best=vm;state=copy.deepcopy(m.state_dict());be=epoch;bad=0
        else:bad+=1
        if bad>=15:break
    m.load_state_dict(state);m.eval()
    mse=[];shs=[];drs=[];rawsh=[];rawdr=[]
    with torch.no_grad():
        for ix in np.array_split(te,int(np.ceil(len(te)/512))):
            bg,sh,dr,e=means(m,X[ix],d[ix]);p=m.decode(bg,sh,dr,b[ix]);mse.append(((p-X[ix])**2).mean(1).numpy());shs.append(sh.numpy());drs.append(dr.numpy());rawsh.append(e["mu_shared"].numpy());rawdr.append(torch.cat(e["mu_drug_list"],1).numpy())
    mse=np.concatenate(mse);shs=np.concatenate(shs);drs=np.concatenate(drs);rawsh=np.concatenate(rawsh);rawdr=np.concatenate(rawdr)
    treated=np.where(obs.iloc[te].drug.astype(str).to_numpy()!="control")[0]
    vs=np.var(shs[treated],0).sum();vd=np.var(drs[treated],0).sum()
    active_shared=int((np.var(rawsh[treated],0)>0.01).sum())
    rows.append({"seed":seed,"best_epoch":be,"test_mse":float(mse.mean()),"shared_fraction":float(vs/(vs+vd)),"active_shared_dimensions":active_shared})
    print(rows[-1],flush=True)

pca=PCA(32,random_state=1729).fit(X_np[tr]);pred=pca.inverse_transform(pca.transform(X_np[te]))
pca_mse=float(mean_squared_error(X_np[te],pred))
df=pd.DataFrame(rows);df.to_csv(O/"training_only_model_gene_sensitivity.csv",index=False)
sumrow=pd.DataFrame([{"mean_test_mse":df.test_mse.mean(),"min_test_mse":df.test_mse.min(),"max_test_mse":df.test_mse.max(),"mean_shared_fraction":df.shared_fraction.mean(),"min_shared_fraction":df.shared_fraction.min(),"max_shared_fraction":df.shared_fraction.max(),"mean_active_shared_dimensions":df.active_shared_dimensions.mean(),"pca_test_mse":pca_mse}])
sumrow.to_csv(O/"training_only_model_gene_sensitivity_summary.csv",index=False)
print("\n",sumrow.to_string(index=False),flush=True)
