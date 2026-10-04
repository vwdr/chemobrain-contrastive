"""Render the Phase 4 (D1) result tables as Markdown fragments from the CSV outputs of 31_semisynthetic.py.

Writes runs/peerreview_20261003/scratch/d1_tables.md (assembled into D1_semisynthetic.md).
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

R = Path(__file__).resolve().parents[1]
O = R / 'analysis' / 'peerreview_20261003'


def f(x, nd=3):
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return 'NA'
    return f'{x:.{nd}f}' if isinstance(x, (float, np.floating)) else str(x)


def main():
    s = pd.read_csv(O / 'D1_config_summary.csv')
    L = ['## Per-configuration summary (`D1_config_summary.csv`; means over seeds 0–2)\n',
         '| config | true shared (energy) | true shared (trace) | F_sh mean | F_sh SD | F_sh min–max | F_sh dox | F_sh cis | test MSE | AU sh/dox/cis | KL sh/dox/cis (nats) | CKA bg/sh/drug |',
         '|---|---|---|---|---|---|---|---|---|---|---|---|']
    for x in s.itertuples():
        L.append(f'| {x.config} | {f(x.true_shared_energy_fraction)} | {f(x.true_shared_trace_fraction)} | {f(x.F_sh_mean)} | {f(x.F_sh_sd)} | '
                 f'{f(x.F_sh_min)}–{f(x.F_sh_max)} | {f(x.F_sh_dox_mean)} | {f(x.F_sh_cis_mean)} | {x.test_mse_mean:.5f} | '
                 f'{x.shared_AU_mean:.2f}/{x.dox_AU_mean:.2f}/{x.cis_AU_mean:.2f} | {x.shared_KL_mean:.4f}/{x.dox_KL_mean:.4f}/{x.cis_KL_mean:.4f} | '
                 f'{x.mean_cross_seed_cka_bg:.3f}/{x.mean_cross_seed_cka_shared:.3f}/{x.mean_cross_seed_cka_drug:.3f} |')
    fits = pd.read_csv(O / 'D1_fit_metrics.csv')
    L += ['\n## Per-fit F_sh (`D1_fit_metrics.csv`)\n', '| config | seed 0 | seed 1 | seed 2 |', '|---|---|---|---|']
    for cfg, g in fits.groupby('config', sort=False):
        v = g.set_index('seed').shared_fraction
        L.append(f'| {cfg} | {v.get(0, np.nan):.4f} | {v.get(1, np.nan):.4f} | {v.get(2, np.nan):.4f} |')
    ig = pd.read_csv(O / 'D1_ig_per_config_mean.csv')
    L += ['\n## Integrated-gradient recovery, per fit, mean over seeds (`D1_ig_per_fit.csv`, `D1_ig_per_config_mean.csv`)\n',
          'Precision = hits/50, recall = hits/40 for the top 50 genes of each block. Cross-routing counts: A∪B genes in the shared '
          'block top 50; S genes in the cisplatin / doxorubicin block top 50.\n',
          '| config | baseline | shared→S P/R | cis→A P/R | dox→B P/R | A∪B in shared top50 | S in cis top50 | S in dox top50 |',
          '|---|---|---|---|---|---|---|---|']
    for x in ig.itertuples():
        L.append(f'| {x.config} | {x.baseline} | {x.shared_vs_S_precision:.3f}/{x.shared_vs_S_recall:.3f} | {x.cis_vs_A_precision:.3f}/{x.cis_vs_A_recall:.3f} | '
                 f'{x.dox_vs_B_precision:.3f}/{x.dox_vs_B_recall:.3f} | {x.shared_top_n_AB_genes:.2f} | {x.cis_top_n_S_genes:.2f} | {x.dox_top_n_S_genes:.2f} |')
    c = pd.read_csv(O / 'D1_ig_consensus.csv')
    L += ['\n## Integrated-gradient recovery, 4-of-6 consensus per configuration (`D1_ig_consensus.csv`)\n',
          'Precision = hits / consensus-list size; recall = hits / 40.\n',
          '| config | n shared/cis/dox | shared→S hits (P/R) | cis→A hits (P/R) | dox→B hits (P/R) | A∪B in shared | S in cis | S in dox |',
          '|---|---|---|---|---|---|---|---|']
    for x in c.itertuples():
        L.append(f'| {x.config} | {x.list_size_shared}/{x.list_size_cis}/{x.list_size_dox} | {x.shared_vs_S_hits} ({f(x.shared_vs_S_precision, 2)}/{x.shared_vs_S_recall:.2f}) | '
                 f'{x.cis_vs_A_hits} ({f(x.cis_vs_A_precision, 2)}/{x.cis_vs_A_recall:.2f}) | {x.dox_vs_B_hits} ({f(x.dox_vs_B_precision, 2)}/{x.dox_vs_B_recall:.2f}) | '
                 f'{x.shared_top_n_AB_genes} | {x.cis_top_n_S_genes} | {x.dox_top_n_S_genes} |')
    p = O / 'D1_permutation_summary.csv'
    if p.exists() and p.stat().st_size > 5:
        pr = pd.read_csv(p)
        L += ['\n## Label-exchange reference (`D1_permutation_summary.csv`, `D1_permutation_values_*.csv`)\n',
              '| config | observed F_sh (seed 0) | n perm | n perm ≥ observed | P | perm mean (SD) | perm 2.5% / 50% / 97.5% | perm min–max |',
              '|---|---|---|---|---|---|---|---|']
        for x in pr.itertuples():
            L.append(f'| {x.config} | {x.observed_seed0_F_sh:.4f} | {x.n_permutations} | {x.n_perm_ge_observed} | {x.p_value:.4f} | '
                     f'{x.perm_F_sh_mean:.4f} ({x.perm_F_sh_sd:.4f}) | {x.perm_F_sh_q025:.4f} / {x.perm_F_sh_median:.4f} / {x.perm_F_sh_q975:.4f} | '
                     f'{x.perm_F_sh_min:.4f}–{x.perm_F_sh_max:.4f} |')
    cal = O / 'D1_calibration.json'
    if cal.exists():
        j = json.loads(cal.read_text())
        L += ['\n## Calibration (`D1_calibration.json`)\n', '```', json.dumps({k: v for k, v in j.items()}, indent=1), '```']
    gs = pd.read_csv(O / 'D1_gene_sets.csv')
    L += ['\n## Injected gene sets (`D1_gene_sets.csv`)\n', '| set | direction | genes (decile) |', '|---|---|---|']
    for (st, dr), g in gs.groupby(['set', 'direction']):
        L.append(f'| {st} | {dr} | ' + ', '.join(f'{a} ({b})' for a, b in zip(g.gene, g.decile)) + ' |')
    out = R / 'runs' / 'peerreview_20261003' / 'scratch' / 'd1_tables.md'
    out.write_text('\n'.join(L) + '\n')
    print(out)


if __name__ == '__main__':
    sys.exit(main())
