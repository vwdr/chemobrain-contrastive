from pathlib import Path
import requests,gzip,shutil,tarfile,concurrent.futures
R=Path(__file__).resolve().parents[1];D=R/'data/raw';E=R/'data/evidence';D.mkdir(parents=True,exist_ok=True);E.mkdir(parents=True,exist_ok=True)
def get(pair):
 u,p=pair
 if not p.exists():
  r=requests.get(u,stream=True,timeout=240);r.raise_for_status()
  with open(str(p)+'.part','wb') as f:
   for c in r.iter_content(2**20):f.write(c)
  Path(str(p)+'.part').rename(p)
 print(p.name,p.stat().st_size,flush=True)
jobs=[]
for acc in ['GSE216146','GSE271055','GSE286221']:jobs.append((f'https://ftp.ncbi.nlm.nih.gov/geo/series/{acc[:-3]}nnn/{acc}/soft/{acc}_family.soft.gz',E/f'{acc}.soft.gz'))
for acc,name in [('GSE216146','GSE216146_chemo_brain.h5ad.gz'),('GSE271055','GSE271055_RAW.tar')]:jobs.append((f'https://ftp.ncbi.nlm.nih.gov/geo/series/{acc[:-3]}nnn/{acc}/suppl/{name}',D/name))
for n in ['m2.cp.reactome.v2025.1.Mm.symbols.gmt','m5.go.bp.v2025.1.Mm.symbols.gmt','m8.all.v2025.1.Mm.symbols.gmt']:jobs.append(('https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2025.1.Mm/'+n,E/n))
with concurrent.futures.ThreadPoolExecutor(5) as ex:list(ex.map(get,jobs))
for p in E.glob('*.gz'):
 with gzip.open(p,'rb') as f,open(p.with_suffix(''),'wb') as g:shutil.copyfileobj(f,g)
p=D/'GSE216146_chemo_brain.h5ad'
if not p.exists():
 with gzip.open(str(p)+'.gz','rb') as f,open(p,'wb') as g:shutil.copyfileobj(f,g)
with tarfile.open(D/'GSE271055_RAW.tar') as t:t.extractall(D/'dox',filter='data')