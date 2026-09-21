"""Conservative reference transfer from source-annotated hippocampal controls."""
from pathlib import Path
import anndata as ad,scanpy as sc,numpy as np,pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score
R=Path(__file__).resolve().parents[1];T=R/'analysis/corrected_20260920';a=ad.read_h5ad(R/'data/processed/recovered_counts.h5ad');sc.pp.normalize_total(a,target_sum=1e4);sc.pp.log1p(a);o=a.obs.copy();ref=np.where((o.study=='GSE216146')&(o.arm=='control'))[0];v=a[ref].copy();sc.pp.highly_variable_genes(v,n_top_genes=1500,flavor='seurat');genes=v.var_names[v.var.highly_variable];x=a[:,genes].X.toarray();y=o.source_cell_type.to_numpy();cv=[]
for sid in o.sample_id.iloc[ref].unique():
 tr=ref[o.sample_id.iloc[ref].to_numpy()!=sid];te=ref[o.sample_id.iloc[ref].to_numpy()==sid];m=LogisticRegression(C=.1,max_iter=600,class_weight='balanced').fit(x[tr],y[tr]);cv.append({'sample_id':sid,'n_cells':len(te),'balanced_accuracy':balanced_accuracy_score(y[te],m.predict(x[te]))})
m=LogisticRegression(C=.1,max_iter=600,class_weight='balanced').fit(x[ref],y[ref]);prob=m.predict_proba(x);pred=m.classes_[prob.argmax(1)];confidence=prob.max(1)
markers={'Neuron':['Snap25','Syt1','Rbfox3','Slc17a7','Gad1','Gad2'],'Astrocyte':['Aldh1l1','Aqp4','Slc1a3','Gfap'],'Oligodendrocyte':['Mbp','Plp1','Mog','Mag'],'OPC':['Pdgfra','Cspg4','Vcan'],'Microglia':['P2ry12','Cx3cr1','Tmem119','C1qa'],'Endothelial':['Cldn5','Pecam1','Kdr','Flt1'],'Pericyte':['Pdgfrb','Rgs5','Kcnj8'],'Ependymal':['Foxj1','Pifo','Dynlrb2'],'Choroid Plexus':['Ttr','Klotho','Aqp1'],'Immune':['Ptprc','Cd3d','Nkg7'],'Progenitor':['Sox2','Dcx','Mki67']}
score=[];mr=[]
for typ,gg in markers.items():
 gg=[g for g in gg if g in a.var_names];vals=a[:,gg].X.toarray();ss=np.zeros(len(a))
 for study in o.study.unique():
  ix=np.where(o.study==study)[0];ss[ix]=((vals[ix]-vals[ix].mean(0))/(vals[ix].std(0)+1e-6)).mean(1)
 score.append(ss)
 for g in gg:mr.append({'cell_type':typ,'gene':g})
marker_label=np.array(list(markers))[np.array(score).argmax(0)];cis=o.study=='GSE216146';o['annotation']=np.where(cis,o.source_cell_type,np.where((confidence>=.8)&(pred==marker_label),pred,'uncertain'));o['annotation_origin']=np.where(cis,'source_study','provisional_transfer');o['reference_prediction']=pred;o['reference_probability']=confidence;o['marker_label']=marker_label;o.to_csv(R/'runs/corrected_20260920/annotated_cells.csv');pd.DataFrame(cv).to_csv(T/'annotation_reference_cv.csv',index=False);pd.DataFrame(mr).to_csv(T/'canonical_markers.csv',index=False);o.groupby(['study','arm','annotation','annotation_origin'],observed=True).size().rename('n_cells').reset_index().to_csv(T/'annotation_counts.csv',index=False);o.groupby(['study','annotation'],observed=True).reference_probability.agg(['count','mean','median']).to_csv(T/'annotation_confidence.csv')
rows=[]
for typ,gg in markers.items():
 for gene in gg:
  if gene not in a.var_names:continue
  val=a[:,gene].X.toarray().ravel()
  for (study,label),idx in o.groupby(['study','annotation'],observed=True).indices.items():rows.append({'study':study,'annotation':label,'marker_class':typ,'gene':gene,'mean_log_expression':val[idx].mean(),'fraction_positive':(val[idx]>0).mean()})
pd.DataFrame(rows).to_csv(T/'marker_validation.csv',index=False);print(o.groupby(['study','annotation']).size(),flush=True)
