"""Task C6: Ttr distribution in GSE216146 by source cell type and arm.

For each source cell type x arm (control, control_GENUS, cisplatin, cisplatin_GENUS; nonrescue and rescue
arms reported separately) and for Choroid Plexus vs all other cell types combined: number of cells,
fraction of cells with nonzero Ttr counts, mean raw Ttr counts, mean log-normalized Ttr (canonical
normalization: 1e4 over all genes, log1p).
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import anndata as ad  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from src.peerreview import core  # noqa: E402


def summarize(df, keys):
    return (df.groupby(keys, observed=True)
            .agg(n_cells=('raw', 'size'), fraction_nonzero=('raw', lambda v: float((v > 0).mean())),
                 mean_raw_counts=('raw', 'mean'), mean_log_normalized=('lognorm', 'mean')).reset_index())


def main():
    a = ad.read_h5ad(core.RAW_H5AD)
    a = a[a.obs.study == 'GSE216146']
    raw = np.asarray(a[:, 'Ttr'].X.toarray()).ravel().astype(np.float64)
    lib = np.asarray(a.X.sum(1)).ravel().astype(np.float64)
    df = pd.DataFrame({'cell_type': a.obs.source_cell_type.astype(str).to_numpy(), 'arm': a.obs.arm.astype(str).to_numpy(),
                       'rescue': a.obs.rescue.astype(bool).to_numpy(), 'raw': raw, 'lognorm': np.log1p(raw / lib * 1e4)})
    by_type = summarize(df, ['cell_type', 'arm'])
    df['cp_group'] = np.where(df.cell_type == 'Choroid Plexus', 'Choroid Plexus', 'All other cell types')
    cp = summarize(df, ['cp_group', 'arm'])
    cp_all = summarize(df, ['cp_group']).assign(arm='all_arms')
    by_type.to_csv(core.PR_OUT / 'C6_ttr_by_celltype_arm.csv', index=False)
    pd.concat([cp, cp_all], ignore_index=True).to_csv(core.PR_OUT / 'C6_ttr_choroid_plexus_vs_other.csv', index=False)
    pd.set_option('display.width', 200); pd.set_option('display.max_rows', 200)
    print(by_type.to_string(index=False)); print(pd.concat([cp, cp_all]).to_string(index=False))


if __name__ == '__main__':
    main()
