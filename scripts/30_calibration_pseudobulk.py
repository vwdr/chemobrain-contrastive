"""Phase 4 step 11: effect-size calibration from the Task 4 paired pseudobulk analysis.

Runs scripts/16_pseudobulk_paired.py unmodified in logic, with its output folder (module global ``T``)
redirected to runs/peerreview_20261003/pseudobulk_task4. Then computes the median |log2FC| of
direction-consistent, BH-significant cisplatin-vs-control rows using the support rule of
scripts/19_gene_pattern_audit.derive_cisplatin_de_support (wald_padj < 0.05 and all three paired
log2CPM differences in the direction of the fitted log2FoldChange).
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.peerreview.threads import configure  # noqa: E402
configure(1)

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from src.peerreview import core  # noqa: E402

T = core.PR_RUN / 'pseudobulk_task4'
GRID = [0.25, 0.5, 1.0, 2.0]


def main():
    T.mkdir(parents=True, exist_ok=True)
    res = T / 'pseudobulk_task4_paired_results.csv.gz'
    if not res.exists():
        m16 = core.load_script('pseudobulk_paired_16', '16_pseudobulk_paired.py')
        m16.T = T
        m16.main()
    x = pd.read_csv(res)
    x = x[x.contrast == 'cisplatin_vs_control'].copy()
    pair_cols = ['replicate_1_log2cpm_difference', 'replicate_2_log2cpm_difference', 'replicate_3_log2cpm_difference']
    same = (np.sign(x[pair_cols].to_numpy()) == np.sign(x.log2FoldChange.to_numpy())[:, None]).all(axis=1)
    x['supported_row'] = (x.wald_padj < 0.05) & same
    sup = x[x.supported_row]
    universe = set(pd.read_csv(core.CORR_OUT / 'benchmark_gene_universe.csv')['gene'].astype(str))
    sup_u = sup[sup.gene.isin(universe)]
    gene_max = sup.groupby('gene').log2FoldChange.apply(lambda v: float(np.max(np.abs(v))))
    gene_max_u = gene_max[gene_max.index.isin(universe)]
    out = {
        'n_supported_rows_all_genes': int(len(sup)), 'n_supported_unique_genes_all': int(sup.gene.nunique()),
        'median_abs_log2FC_supported_rows_all_genes': float(sup.log2FoldChange.abs().median()),
        'n_supported_rows_in_1500_universe': int(len(sup_u)), 'n_supported_unique_genes_in_1500_universe': int(sup_u.gene.nunique()),
        'median_abs_log2FC_supported_rows_in_1500_universe': float(sup_u.log2FoldChange.abs().median()),
        'median_of_per_gene_max_abs_log2FC_all': float(gene_max.median()),
        'median_of_per_gene_max_abs_log2FC_in_1500_universe': float(gene_max_u.median()),
        'quartiles_abs_log2FC_supported_rows_all_genes': [float(q) for q in sup.log2FoldChange.abs().quantile([.25, .5, .75])],
        'delta_grid': GRID,
        'support_rule': 'wald_padj < 0.05 and all three paired log2CPM differences share the sign of log2FoldChange',
        'pydeseq2_outputs': str(T.relative_to(core.R)),
    }
    for k in ['median_abs_log2FC_supported_rows_all_genes', 'median_abs_log2FC_supported_rows_in_1500_universe']:
        v = out[k]; lower = [g for g in GRID if g <= v]; upper = [g for g in GRID if g >= v]
        out[k + '_grid_position'] = {'nearest_grid_value': min(GRID, key=lambda g: abs(g - v)),
                                     'between': [max(lower) if lower else None, min(upper) if upper else None]}
    sup.to_csv(core.PR_OUT / 'D1_calibration_supported_rows.csv', index=False)
    core.write_json(core.PR_OUT / 'D1_calibration.json', out)
    print(out)


if __name__ == '__main__':
    main()
