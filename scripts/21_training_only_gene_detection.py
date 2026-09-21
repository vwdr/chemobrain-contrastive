"""Training-only gene-detection sensitivity for the canonical model universe.

The historical preprocessing required genes to be detected in at least ten
observations before the train/validation/test split. This diagnostic asks
whether the final 1,500-gene model universe would also satisfy that detection
criterion using canonical training cells alone, both overall and within each
source study. It does not change the archived fits.
"""
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split

R=Path(__file__).resolve().parents[1]
T=R/"analysis"/"corrected_20260920"
SEED=1729

def split(obs):
    pure=np.where(~obs.rescue.astype(bool).to_numpy())[0]
    strata=(obs.study.astype(str)+"_"+obs.drug.astype(str)).to_numpy()
    tr,te=train_test_split(pure,test_size=.2,stratify=strata[pure],random_state=SEED)
    va,te=train_test_split(te,test_size=.5,stratify=strata[te],random_state=SEED)
    rng=np.random.default_rng(SEED)
    tr=np.concatenate([
        rng.choice(tr[strata[tr]==s],min(1500,int(np.sum(strata[tr]==s))),replace=False)
        for s in np.unique(strata[tr])
    ])
    return tr.astype(int)

a=ad.read_h5ad(R/"data"/"processed"/"recovered_counts.h5ad")
tr=split(a.obs)
genes=pd.read_csv(T/"benchmark_gene_universe.csv").gene.astype(str).tolist()
x=a[tr,genes].X
overall=np.asarray((x>0).sum(0)).ravel()
rows=[]
for study in sorted(a.obs.study.astype(str).unique()):
    ix=tr[a.obs.iloc[tr].study.astype(str).to_numpy()==study]
    z=a[ix,genes].X
    det=np.asarray((z>0).sum(0)).ravel()
    for g,n in zip(genes,det):
        rows.append({"study":study,"gene":g,"training_cells_detected":int(n)})
pd.DataFrame(rows).to_csv(T/"training_only_gene_detection_by_study.csv",index=False)
summary=pd.DataFrame([
    {"scope":"all training cells","n_genes":len(genes),"minimum_detected_cells":int(overall.min()),"n_below_10":int((overall<10).sum())},
    *[
        {"scope":study,"n_genes":len(genes),"minimum_detected_cells":int(g.training_cells_detected.min()),"n_below_10":int((g.training_cells_detected<10).sum())}
        for study,g in pd.DataFrame(rows).groupby("study")
    ]
])
summary.to_csv(T/"training_only_gene_detection_summary.csv",index=False)
print(summary.to_string(index=False))
