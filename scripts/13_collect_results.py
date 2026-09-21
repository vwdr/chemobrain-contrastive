from pathlib import Path
import json
import pandas as pd
R=Path(__file__).resolve().parents[1]
rows=[json.loads(p.read_text()) for p in sorted((R/'runs/corrected_20260920').glob('*_metrics.json')) if p.name!='pca_metrics.json']
assert len(rows)==15
pd.DataFrame(rows).to_csv(R/'analysis/corrected_20260920/model_benchmarks.csv',index=False)
