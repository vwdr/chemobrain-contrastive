"""Task C4: multiGroupVI baseline (Weinberger, Lopez, Huetter & Regev, MLCB 2022 / PMLR 200).

Run with the separate environment .venv-mgvi (scvi-tools 0.18.0 + multigroup_vi from
github.com/Genentech/multiGroupVI @ 3d001dbe6905803195165f936f879a1a166a6fc8, installed unmodified).

Groups: 0 = control (both studies), 1 = cisplatin, 2 = doxorubicin. Same canonical 6,000 training cells,
1,500 genes and canonical test cells; seeds 0-2; package defaults (n_hidden 128, n_layers 1, shared 10,
private 10, wasserstein_penalty 1, batch size 128, train_size 0.9, max_epochs = package default rule,
no early stopping); no hyperparameter search. Inputs are read from the regenerated canonical npz files so
the old anndata in .venv-mgvi never reads the h5ad.

Usage (with .venv-mgvi python):
  python scripts/34_c4_multigroupvi.py --fit --seed 0 [--threads=3]
  python scripts/34_c4_multigroupvi.py --collect
"""
import os
import sys
from pathlib import Path

N = 3
for a in sys.argv[1:]:
    if a.startswith('--threads='):
        N = int(a.split('=')[1])
for k in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS'):
    os.environ[k] = str(N)

import argparse  # noqa: E402
import itertools  # noqa: E402
import json  # noqa: E402
import time  # noqa: E402

import anndata as ad  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import torch  # noqa: E402
from sklearn.linear_model import LogisticRegression  # noqa: E402
from sklearn.metrics import balanced_accuracy_score  # noqa: E402

torch.set_num_threads(N)
R = Path(__file__).resolve().parents[1]
IN = R / 'runs' / 'peerreview_20261003' / 'canonical_inputs'
RUN = R / 'runs' / 'peerreview_20261003' / 'multigroupvi'
OUT = R / 'analysis' / 'peerreview_20261003'
GROUP_OF_D = {0: 0, 2: 1, 1: 2}          # d codes: 0 control, 1 doxorubicin, 2 cisplatin -> groups
GROUP_NAMES = ['control', 'cisplatin', 'doxorubicin']


def modeled_gene_log(counts):
    counts = np.asarray(counts, dtype=np.float64); total = counts.sum(1, keepdims=True); total[total <= 0] = 1.0
    return np.log1p(counts / total * 1e4).astype(np.float32)


def probe(xtr, xte, ytr, yte, max_iter=1500):
    c = LogisticRegression(max_iter=max_iter, class_weight='balanced').fit(xtr, ytr)
    return float(balanced_accuracy_score(yte, c.predict(xte)))


def cka(x, y):
    x = x - x.mean(0); y = y - y.mean(0)
    den = np.linalg.norm(x.T @ x, 'fro') * np.linalg.norm(y.T @ y, 'fro')
    return float(np.linalg.norm(x.T @ y, 'fro') ** 2 / den) if den > 0 else float('nan')


def per_dim_kl(m, v):
    return 0.5 * (m ** 2 + v - 1.0 - np.log(v))


def load():
    z = np.load(IN / 'inputs.npz'); q = np.load(IN / 'count_inputs.npz')
    return z, q['C'], q['library']


def make_adata(C, idx, d, b, genes):
    obs = pd.DataFrame({'group': pd.Categorical([GROUP_NAMES[GROUP_OF_D[int(x)]] for x in d[idx]], categories=GROUP_NAMES),
                        'study': np.where(b[idx] == 1, 'GSE271055', 'GSE216146')}, index=[f'c{i}' for i in idx])
    return ad.AnnData(X=np.rint(C[idx]).astype(np.float32), obs=obs, var=pd.DataFrame(index=list(genes)))


def posterior(model, adata):
    from scvi import REGISTRY_KEYS
    from scvi.dataloaders import AnnDataLoader
    adata = model._validate_anndata(adata)
    dl = model._make_data_loader(adata=adata, batch_size=512, shuffle=False, data_loader_class=AnnDataLoader)
    out = {k: [] for k in ('qz_m', 'qz_v', 'qt_m_all', 'qt_v_all', 'mask')}
    model.module.eval()
    with torch.no_grad():
        for t in dl:
            o = model.module.inference(x=t[REGISTRY_KEYS.X_KEY], group_labels=t['group'])
            for k in out:
                out[k].append(o[k].cpu().numpy())
    return {k: np.concatenate(v) for k, v in out.items()}


def fit(seed, smoke=False):
    global RUN
    if smoke:
        RUN = RUN / 'smoke'
    import scvi
    from multigroup_vi.model import MultiGroupVI as MultiGroupVIModel
    RUN.mkdir(parents=True, exist_ok=True)
    done = RUN / f'mgvi_{seed}_done.json'
    if done.exists():
        return
    z, C, lib = load(); d, b, tr, te, genes = z['d'], z['b'], z['train'], z['test'], z['genes']
    t0 = time.time()
    scvi.settings.seed = seed; torch.manual_seed(seed); np.random.seed(seed)
    atr = make_adata(C, tr, d, b, genes); ate = make_adata(C, te, d, b, genes)
    MultiGroupVIModel.setup_anndata(atr, group_key='group')
    gcodes = atr.obs['group'].cat.codes.to_numpy()
    group_indices = [np.where(gcodes == g)[0] for g in range(3)]
    model = MultiGroupVIModel(atr, n_groups=3)
    # package default rule min(round(20000 / n_cells * 400), 400), passed as a Python int because the package
    # computes it with np.min (numpy integer), which pytorch-lightning 1.7.7 rejects (assert isinstance(int)).
    max_epochs = int(min(round((20000 / atr.n_obs) * 400), 400)) if not smoke else 2
    model.train(group_indices, max_epochs=max_epochs, use_gpu=False)
    model_dir = RUN / f'mgvi_{seed}_model'
    if not model_dir.exists():
        model.save(str(model_dir), save_anndata=False)
    epochs = int(model.trainer.current_epoch) + 1 if hasattr(model, 'trainer') else None
    # posteriors
    ptr = posterior(model, atr); pte = posterior(model, ate)
    pred = model.get_normalized_expression(ate, library_size=1e4, return_numpy=True)
    obs_mod = modeled_gene_log(ate.X)
    pred_mod = np.log1p(np.asarray(pred, dtype=np.float64)).astype(np.float32)
    l_model = np.asarray(ate.X, dtype=np.float64).sum(1, keepdims=True)
    pred_can = np.log1p(np.asarray(pred, dtype=np.float64) * l_model / lib[te][:, None]).astype(np.float32)
    np_ = model.module.n_private_latent
    reps_tr = {'shared': ptr['qz_m'], 'private_gated': ptr['qt_m_all'] * ptr['mask'], 'private_ungated': ptr['qt_m_all']}
    reps_te = {'shared': pte['qz_m'], 'private_gated': pte['qt_m_all'] * pte['mask'], 'private_ungated': pte['qt_m_all']}
    row = {'method': 'multiGroupVI', 'seed': seed, 'epochs_run': epochs, 'n_train': len(tr), 'n_test': len(te),
           'modeled_gene_log_mse': float(np.mean((pred_mod - obs_mod) ** 2)),
           'canonical_log_mse': float(np.mean((pred_can - z['X'][te]) ** 2))}
    for nm, (xr, xe) in {k: (reps_tr[k], reps_te[k]) for k in reps_tr}.items():
        row[f'{nm}_study_balanced_accuracy'] = probe(xr, xe, b[tr], b[te])
        row[f'{nm}_condition_balanced_accuracy'] = probe(xr, xe, d[tr], d[te])
        for study, bv in (('cisplatin_within_GSE216146', 0), ('doxorubicin_within_GSE271055', 1)):
            mtr = b[tr] == bv; mte = b[te] == bv
            row[f'{nm}_treatment_probe_{study}_balanced_accuracy'] = probe(xr[mtr], xe[mte], (d[tr] > 0)[mtr].astype(int),
                                                                           (d[te] > 0)[mte].astype(int), max_iter=800)
    # utilization on held-out cells
    dims = []
    gte = np.array([GROUP_OF_D[int(x)] for x in d[te]])
    blocks = [('shared', 'all_test', pte['qz_m'], pte['qz_v']),
              ('shared', 'treated_test', pte['qz_m'][d[te] > 0], pte['qz_v'][d[te] > 0])]
    for g in range(3):
        sl = slice(g * np_, (g + 1) * np_)
        blocks.append((f'private_{GROUP_NAMES[g]}', f'{GROUP_NAMES[g]}_test', pte['qt_m_all'][gte == g, sl], pte['qt_v_all'][gte == g, sl]))
    for blk, pop, m, v in blocks:
        kl = per_dim_kl(m.astype(np.float64), v.astype(np.float64)).mean(0); var = m.var(0)
        row[f'{blk}__{pop}__active_units_var_gt_0p01'] = int((var > 0.01).sum())
        row[f'{blk}__{pop}__dims_kl_gt_0p01'] = int((kl > 0.01).sum())
        row[f'{blk}__{pop}__total_kl_nats'] = float(kl.sum())
        row[f'{blk}__{pop}__n_dims'] = int(len(kl))
        for j in range(len(kl)):
            dims.append(dict(seed=seed, block=blk, population=pop, dimension=j, n_cells=len(m), mean_kl_nats=float(kl[j]),
                             posterior_mean_variance=float(var[j]), active_unit_var_gt_0p01=bool(var[j] > 0.01)))
    pd.DataFrame(dims).to_csv(RUN / f'mgvi_{seed}_usage_per_dimension.csv', index=False)
    np.savez_compressed(RUN / f'mgvi_{seed}_latents.npz', **{f'test_{k}': v for k, v in reps_te.items()})
    hist = getattr(model, 'history', None)
    if hist:
        pd.concat([v for v in hist.values()], axis=1).to_csv(RUN / f'mgvi_{seed}_history.csv')
    row['fit_seconds'] = time.time() - t0; row['threads'] = N
    done.write_text(json.dumps(row, indent=2) + '\n'); print(row, flush=True)


def collect():
    rows = [json.loads((RUN / f'mgvi_{s}_done.json').read_text()) for s in (0, 1, 2)]
    df = pd.DataFrame(rows); df.to_csv(OUT / 'C4_multigroupvi_metrics.csv', index=False)
    lat = {s: np.load(RUN / f'mgvi_{s}_latents.npz') for s in (0, 1, 2)}
    ck = []
    for a, c in itertools.combinations((0, 1, 2), 2):
        for k in ('test_shared', 'test_private_gated', 'test_private_ungated'):
            ck.append(dict(seed_a=a, seed_b=c, representation=k.replace('test_', ''), linear_cka=cka(lat[a][k], lat[c][k])))
    pd.DataFrame(ck).to_csv(OUT / 'C4_multigroupvi_cross_seed_cka.csv', index=False)
    pd.concat([pd.read_csv(RUN / f'mgvi_{s}_usage_per_dimension.csv') for s in (0, 1, 2)]).to_csv(
        OUT / 'C4_multigroupvi_usage_per_dimension.csv', index=False)
    pd.set_option('display.width', 250)
    print(df.T.to_string()); print(pd.DataFrame(ck).to_string(index=False))


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--fit', action='store_true'); ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--collect', action='store_true'); ap.add_argument('--smoke', action='store_true')
    a = ap.parse_args([x for x in sys.argv[1:] if not x.startswith('--threads=')])
    if a.fit:
        fit(a.seed, smoke=a.smoke)
    if a.collect:
        collect()
