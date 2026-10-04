"""R4: multiGroupVI on the semi-synthetic data (run with .venv-mgvi).

Groups: 0 pseudo-control (both studies), 1 pseudo-cisplatin (GSE216146 pseudo-treated), 2 pseudo-doxorubicin
(GSE271055 pseudo-treated). Phase 4 training cells (6,000) and test cells (1,429); injected raw counts from
scripts/41_r4_prepare.py. Settings as Task C4 (package defaults; max_epochs = package default rule as a Python int,
400 for 6,000 cells). Helper functions (make_adata, posterior, probe, per_dim_kl) are imported from
scripts/34_c4_multigroupvi.py.

Metrics per fit: canonical-scale MSE; active units / total KL of the private block of each treated group on that
group's test cells; within-study probes (pseudo-treated vs pseudo-control) on the ungated private block of the study's
treated group; integrated gradients of ||E q(t_g | x)||^2 for the treated group g with respect to the group encoder's
input (the package feeds raw counts to the group-specific encoders), zero and control-median baselines (median raw
counts of the study's training pseudo-controls), 128 test pseudo-treated cells per group drawn with default_rng(314)
(core.attribution_indices); top-50 hits against S u A (pseudo-cisplatin) or S u B (pseudo-doxorubicin).

Usage: .venv-mgvi/Scripts/python scripts/42_r4_multigroupvi.py --config shared_d1_rf1.0 --seed 0 [--threads=3]
"""
import sys
from pathlib import Path
R = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(R))
import importlib.util  # noqa: E402

spec = importlib.util.spec_from_file_location('c4', R / 'scripts' / '34_c4_multigroupvi.py')
c4 = importlib.util.module_from_spec(spec); spec.loader.exec_module(c4)   # sets thread env vars from --threads=

import argparse  # noqa: E402
import json  # noqa: E402
import time  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import torch  # noqa: E402

from src.peerreview.ig_generic import integrated_gradients, recovery  # noqa: E402

RUN = R / 'runs' / 'peerreview_20261004' / 'r4' / 'multigroupvi'
GS = pd.read_csv(R / 'analysis' / 'peerreview_20261003' / 'D1_gene_sets.csv')
SETS = {s: set(GS[GS.set == s].gene) for s in ('S', 'A', 'B')}


def fit(cfg, seed, smoke=False):
    global RUN
    if smoke:
        RUN = RUN / 'smoke'
    import scvi
    from multigroup_vi.model import MultiGroupVI
    RUN.mkdir(parents=True, exist_ok=True)
    name = f'{cfg}_s{seed}'; done = RUN / f'{name}_done.json'
    if done.exists():
        print('done marker exists', name); return
    z = np.load(R / 'runs' / 'peerreview_20261004' / 'r4' / f'inputs_{cfg}.npz')
    C, L, X, d, b, tr, te, genes = z['C'], z['L'], z['X'], z['d'], z['b'], z['train'], z['test'], z['genes']
    t0 = time.time()
    scvi.settings.seed = seed; torch.manual_seed(seed); np.random.seed(seed)
    atr = c4.make_adata(C, tr, d, b, genes); ate = c4.make_adata(C, te, d, b, genes)
    MultiGroupVI.setup_anndata(atr, group_key='group')
    gcodes = atr.obs['group'].cat.codes.to_numpy()
    model = MultiGroupVI(atr, n_groups=3)
    max_epochs = int(min(round((20000 / atr.n_obs) * 400), 400)) if not smoke else 2
    model.train([np.where(gcodes == g)[0] for g in range(3)], max_epochs=max_epochs, use_gpu=False)
    model.save(str(RUN / f'{name}_model'), save_anndata=False, overwrite=True)
    fit_s = time.time() - t0
    ptr = c4.posterior(model, atr); pte = c4.posterior(model, ate)
    pred = model.get_normalized_expression(ate, library_size=1e4, return_numpy=True)
    l_model = np.asarray(ate.X, dtype=np.float64).sum(1, keepdims=True)
    pred_can = np.log1p(np.asarray(pred, dtype=np.float64) * l_model / L[te][:, None])
    row = {'model': 'multiGroupVI', 'config': cfg, 'seed': seed, 'epochs': max_epochs, 'fit_seconds': fit_s,
           'canonical_log_mse': float(np.mean((pred_can - X[te]) ** 2))}
    npriv = model.module.n_private_latent
    igrows = []; comp = []
    for g, study_b, label, targets, idx_key in ((1, 0, 'cisplatin', {'S': SETS['S'], 'A': SETS['A']}, 'idx_drug1'),
                                                (2, 1, 'doxorubicin', {'S': SETS['S'], 'B': SETS['B']}, 'idx_drug0')):
        sl = slice(g * npriv, (g + 1) * npriv)
        # utilization on the group's test cells
        gm = (d[te] > 0) & (b[te] == study_b)
        m_, v_ = pte['qt_m_all'][gm, sl].astype(np.float64), pte['qt_v_all'][gm, sl].astype(np.float64)
        kl = c4.per_dim_kl(m_, v_).mean(0); var = m_.var(0)
        row[f'{label}_private_active_units'] = int((var > 0.01).sum()); row[f'{label}_private_total_kl_nats'] = float(kl.sum())
        row[f'{label}_private_dims_kl_gt_0p01'] = int((kl > 0.01).sum())
        # within-study probe on the ungated private block of the treated group
        mtr = b[tr] == study_b; mte = b[te] == study_b
        row[f'probe_within_{label}_study_private_ungated'] = c4.probe(ptr['qt_m_all'][mtr, sl], pte['qt_m_all'][mte, sl],
                                                                      (d[tr] > 0)[mtr].astype(int), (d[te] > 0)[mte].astype(int), max_iter=800)
        row[f'probe_within_{label}_study_all_private_ungated'] = c4.probe(ptr['qt_m_all'][mtr], pte['qt_m_all'][mte],
                                                                          (d[tr] > 0)[mtr].astype(int), (d[te] > 0)[mte].astype(int), max_iter=800)
        # attribution
        enc = model.module.t_encoders[g]; enc.eval(); model.module.eval()

        def fun(x, enc=enc):
            return enc(x)[0].square().sum(1)
        cells = z[idx_key]
        controls = tr[(d[tr] == 0) & (b[tr] == study_b)]
        for base_name, baseline in (('zero', np.zeros(C.shape[1], np.float32)), ('control_median', np.median(C[controls], axis=0).astype(np.float32))):
            ig, info = integrated_gradients(fun, C[cells], baseline)
            mabs = np.abs(ig).mean(0)
            rec, top = recovery(mabs, genes, targets)
            igrows.append(dict(model='multiGroupVI', config=cfg, seed=seed, group=label, baseline=base_name, encoder_input='raw counts',
                               **rec, **info, top50_genes='|'.join(top)))
            pd.DataFrame({'gene': genes, 'mean_abs_ig': mabs}).to_csv(RUN / f'{name}_{label}_{base_name}_ig.csv.gz', index=False, compression='gzip')
    pd.DataFrame(igrows).to_csv(RUN / f'{name}_ig_recovery.csv', index=False)
    row['seconds_total'] = time.time() - t0
    done.write_text(json.dumps(row, indent=2) + '\n'); print(row, flush=True)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--config', required=True); ap.add_argument('--seed', type=int, default=0); ap.add_argument('--smoke', action='store_true')
    a = ap.parse_args([x for x in sys.argv[1:] if not x.startswith('--threads=')])
    fit(a.config, a.seed, smoke=a.smoke)
