"""Phase 4 diagnostic: presence of the injected signal in the model inputs vs its IG ranks (seed-0 fits).

For each configuration: among training cells, |mean(pseudo-treated) - mean(pseudo-control)| of the
log-normalized input per gene, computed within the study where the set is injected (S: both studies
pooled; A: GSE216146; B: GSE271055); the number of set genes in the top 50 of that ranking and their
median rank; the median IG rank (mean |IG|, zero baseline) of the set genes in the matching block; and
the metrics-JSON condition probes of the seed-0 fit.
"""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

R = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(R))
spec = importlib.util.spec_from_file_location('ss', R / 'scripts' / '31_semisynthetic.py')
ss = importlib.util.module_from_spec(spec); sys.argv = sys.argv[:1]; spec.loader.exec_module(ss)


def main():
    base = ss.load_base(); gs = pd.read_csv(ss.OUT / 'D1_gene_sets.csv'); rows = []
    for cfg in ss.all_configs():
        inp, _ = ss.build_inputs(cfg, base); X, d, b, tr = inp['X'], inp['d'], inp['b'], inp['train']
        m = json.loads((ss.RUN / cfg / f'{cfg}_s0_metrics.json').read_text())
        ig = pd.read_csv(ss.RUN / cfg / f'{cfg}_s0_ig_summary.csv.gz'); ig = ig[ig.baseline == 'zero']
        for st, axis, study in (('S', 'shared', None), ('A', 'drug1', 0), ('B', 'drug0', 1)):
            cells = tr if study is None else tr[b[tr] == study]
            diff = np.abs(X[cells[d[cells] > 0]].mean(0) - X[cells[d[cells] == 0]].mean(0))
            rank = (-diff).argsort().argsort()
            g = ig[ig.axis == axis].reset_index(drop=True)
            igrank = (-g.mean_abs_ig.to_numpy()).argsort().argsort()
            idx = gs[gs.set == st].gene_index.to_numpy()
            rows.append(dict(config=cfg, set=st, block=ss.core.AXIS_NAME[axis], n_set_in_top50_input_mean_difference=int((rank[idx] < 50).sum()),
                             median_rank_input_mean_difference=float(np.median(rank[idx])), n_set_in_top50_IG=int((igrank[idx] < 50).sum()),
                             median_rank_IG=float(np.median(igrank[idx])),
                             bg_condition_probe=m['bg_condition_balanced_accuracy'], shared_ungated_condition_probe=m['shared_ungated_condition_balanced_accuracy'],
                             drug_ungated_condition_probe=m['drug_ungated_condition_balanced_accuracy']))
    df = pd.DataFrame(rows); df.to_csv(ss.OUT / 'D1_signal_check_seed0.csv', index=False)
    pd.set_option('display.width', 250); print(df.to_string(index=False))


if __name__ == '__main__':
    main()
