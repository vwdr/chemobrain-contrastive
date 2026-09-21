"""Recover integer counts and correct source labels using archived GEO metadata."""
from pathlib import Path
import h5py,numpy as np,pandas as pd,anndata as ad,scanpy as sc
from scipy import sparse,io
import gzip
ad.settings.allow_write_nullable_strings=True
R=Path(__file__).resolve().parents[1];D=R/'data/raw';O=R/'analysis/corrected_20260920';Q=[]
with h5py.File(D/'GSE216146_chemo_brain.h5ad','r') as f:a=ad.AnnData(ad.io.read_elem(f['X']),obs=ad.io.read_elem(f['obs']),var=ad.io.read_elem(f['var']))
print('Cis columns',a.obs.columns.tolist(),flush=True)
def qc(a):
 a.var_names_make_unique();a.var['mt']=a.var_names.str.startswith(('mt-','MT-','Mt-'));sc.pp.calculate_qc_metrics(a,qc_vars=['mt'],inplace=True,percent_top=None,log1p=False);pre=a.obs.copy();sc.pp.filter_cells(a,min_genes=500);sc.pp.filter_genes(a,min_cells=10);a=a[a.obs.pct_counts_mt<10].copy()
 for (study,sid,arm),s in pre.groupby(['study','sample_id','arm'],observed=True):Q.append({'study':study,'sample_id':sid,'arm':arm,'before_qc':len(s),'after_qc':int((a.obs.sample_id==sid).sum()),'median_counts_before':s.total_counts.median(),'median_pct_mt_before':s.pct_counts_mt.median()})
 return a
# Original sample code is source metadata, not the previous registry interpretation.
a.obs['study']='GSE216146';a.obs['source_code']=a.obs.pathology.astype(str);a.obs['sample_id']=a.obs['biosample'].astype(str)
a.obs['source_cell_type']=a.obs['MajorCellType'].astype(str);a.obs['rescue']=a.obs.source_code.isin(['PS','CS']);a.obs['drug']=a.obs.source_code.map({'PN':'control','PS':'control','CN':'cisplatin','CS':'cisplatin'});assert a.obs.drug.notna().all();a.obs['historical_drug']=a.obs.source_code.map({'PN':'cisplatin','PS':'cisplatin','CN':'control','CS':'control'});a.obs['arm']=a.obs.source_code.map({'PN':'control','PS':'control_GENUS','CN':'cisplatin','CS':'cisplatin_GENUS'});a=qc(a)
ls=[]
for gsm,code,drug,rescue in [('GSM8367861','CNT','control',False),('GSM8367862','DOX','doxorubicin',False),('GSM8367863','DOXACY','doxorubicin',True)]:
 m=next((D/'dox').glob(gsm+'*matrix.mtx.gz'));f=next((D/'dox').glob(gsm+'*features.tsv.gz'));b=next((D/'dox').glob(gsm+'*barcodes.tsv.gz'));x=io.mmread(gzip.open(m)).T.tocsr();genes=pd.read_csv(f,sep='\t',header=None);bar=pd.read_csv(b,sep='\t',header=None)[0].astype(str);v=ad.AnnData(x,obs=pd.DataFrame(index=(bar+'_'+gsm).to_numpy()),var=pd.DataFrame(index=genes[1].astype(str).to_numpy()));v.var_names_make_unique();v.obs['study']='GSE271055';v.obs['sample_id']=gsm;v.obs['source_cell_type']='unannotated';v.obs['source_code']=code;v.obs['drug']=drug;v.obs['historical_drug']=drug;v.obs['rescue']=rescue;v.obs['arm']='doxorubicin_ACY1083' if rescue else drug;ls.append(v)
b=qc(ad.concat(ls,join='outer',fill_value=0));cols=['study','sample_id','source_cell_type','source_code','drug','historical_drug','rescue','arm','total_counts','pct_counts_mt','n_genes_by_counts'];a.obs=a.obs[cols].copy();b.obs=b.obs[cols].copy();a.var=pd.DataFrame(index=a.var_names);b.var=pd.DataFrame(index=b.var_names);v=ad.concat([a,b],join='inner',index_unique='-');v.X=v.X.astype('float32');assert np.allclose(v.X.data,np.round(v.X.data));(R/'data/processed').mkdir(exist_ok=True);v.write_h5ad(R/'data/processed/recovered_counts.h5ad',compression='gzip');pd.DataFrame(Q).to_csv(O/'cohort_qc.csv',index=False);v.obs.groupby(['study','sample_id','arm','source_cell_type'],observed=True).size().rename('n_cells').reset_index().to_csv(O/'source_annotation_counts.csv',index=False);v.obs.groupby(['study','source_code','drug','historical_drug','rescue'],observed=True).size().rename('n_cells').reset_index().to_csv(O/'label_audit.csv',index=False);print(v,flush=True)
