"""R4: pairwise contrastiveVI on the semi-synthetic data (main venv, scvi-tools 1.5.1).

Per study: background = pseudo-control cells, target = pseudo-treated cells of that study; Phase 4 training cells
of the study (3,000) and test cells of the study; injected raw counts from scripts/41_r4_prepare.py. Settings as
Task C3 / scripts/20 (n_hidden 128, background 16, salient 8, 2 layers, dropout 0.1, observed library size,
wasserstein_penalty 0, 100 epochs, batch size 128, train_size 1.0, no early stopping).

Metrics: canonical-scale MSE (background expression for pseudo-controls, salient for pseudo-treated; conversion
log1p(rho*1e4 * L_model / L_full)); salient active units / total KL on test pseudo-treated cells; within-study probe
(pseudo-treated vs pseudo-control, LogisticRegression balanced, max_iter=800) on the ungated salient means;
integrated gradients of ||E q(s | x)||^2 with respect to the salient encoder input (log1p counts, as fed by the model),
zero and control-median baselines (median log1p counts of the study's training pseudo-controls), 128 test
pseudo-treated cells (core.attribution_indices, default_rng(314)); top-50 hits vs S u A (GSE216146) or S u B (GSE271055).

Usage: python scripts/43_r4_contrastivevi.py --config shared_d1_rf1.0 --study GSE216146 --seed 0 [--threads=3]
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.peerreview.threads import configure, strip_threads_arg  # noqa: E402
N_THREADS = configure(3)

import argparse  # noqa: E402
import json  # noqa: E402
import time  # noqa: E402

import anndata as ad  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import torch  # noqa: E402
from sklearn.linear_model import LogisticRegression  # noqa: E402
from sklearn.metrics import balanced_accuracy_score  # noqa: E402

from src.peerreview.ig_generic import integrated_gradients, recovery  # noqa: E402

torch.set_num_threads(N_THREADS)
R = Path(__file__).resolve().parents[1]
RUN = R / 'runs' / 'peerreview_20261004' / 'r4' / 'contrastivevi'
GS = pd.read_csv(R / 'analysis' / 'peerreview_20261003' / 'D1_gene_sets.csv')
SETS = {s: set(GS[GS.set == s].gene) for s in ('S', 'A', 'B')}
STUDY = {'GSE216146': (0, 'idx_drug1', 'cisplatin', ('S', 'A')), 'GSE271055': (1, 'idx_drug0', 'doxorubicin', ('S', 'B'))}


def per_dim_kl(m, v):
    return 0.5 * (m ** 2 + v - 1.0 - np.log(v))


def make(C, idx, treated, genes):
    obs = pd.DataFrame({'pseudo': np.where(treated[idx], 'treated', 'control')}, index=[f'c{i}' for i in idx])
    return ad.AnnData(X=np.rint(C[idx]).astype(np.int64), obs=obs, var=pd.DataFrame(index=list(genes)))


def fit(cfg, study, seed, smoke=False):
    import scvi
    scvi.settings.num_threads = N_THREADS
    run = RUN / 'smoke' if smoke else RUN; run.mkdir(parents=True, exist_ok=True)
    name = f'{cfg}_{study}_s{seed}'; done = run / f'{name}_done.json'
    if done.exists():
        print('done marker exists', name); return
    z = np.load(R / 'runs' / 'peerreview_20261004' / 'r4' / f'inputs_{cfg}.npz')
    C, L, X, d, b, tr, te, genes, treated = z['C'], z['L'], z['X'], z['d'], z['b'], z['train'], z['test'], z['genes'], z['treated']
    sb, idx_key, label, set_names = STUDY[study]
    trs = tr[b[tr] == sb]; tes = te[b[te] == sb]
    t0 = time.time()
    scvi.settings.seed = seed; torch.manual_seed(seed); np.random.seed(seed)
    atr = make(C, trs, treated, genes); ate = make(C, tes, treated, genes)
    scvi.external.ContrastiveVI.setup_anndata(atr)
    bidx = np.where(atr.obs.pseudo.to_numpy() == 'control')[0]; tidx = np.where(atr.obs.pseudo.to_numpy() == 'treated')[0]
    model = scvi.external.ContrastiveVI(atr, n_hidden=128, n_background_latent=16, n_salient_latent=8, n_layers=2,
                                        dropout_rate=0.1, use_observed_lib_size=True, wasserstein_penalty=0)
    model.train(background_indices=bidx, target_indices=tidx, max_epochs=2 if smoke else 100, train_size=1.0, validation_size=0.0,
                early_stopping=False, accelerator='cpu', devices=1, batch_size=128, enable_progress_bar=False)
    model.save(str(run / f'{name}_model'), overwrite=True, save_anndata=False)
    fit_s = time.time() - t0
    ytr = treated[trs].astype(int); yte = treated[tes].astype(int)
    parts = model.get_normalized_expression(ate, library_size=1e4, n_samples=1, return_mean=True, return_numpy=True)
    pred = np.asarray(parts['background'], dtype=np.float64).copy(); sal = np.asarray(parts['salient'], dtype=np.float64)
    pred[yte == 1] = sal[yte == 1]
    l_model = np.asarray(ate.X, dtype=np.float64).sum(1, keepdims=True)
    pred_can = np.log1p(pred * l_model / L[tes][:, None])
    module = model.module; module.eval()
    def post(adata):
        with torch.no_grad():
            o = module._generic_inference(x=torch.tensor(np.asarray(adata.X, dtype=np.float32)), batch_index=torch.zeros((adata.n_obs, 1), dtype=torch.long))
        return o['qs_m'].numpy(), o['qs_v'].numpy()
    mtr, _ = post(atr); mte, vte = post(ate)
    tm = yte == 1
    kl = per_dim_kl(mte[tm].astype(np.float64), vte[tm].astype(np.float64)).mean(0); var = mte[tm].var(0)
    clf = LogisticRegression(max_iter=800, class_weight='balanced').fit(mtr, ytr)
    row = {'model': 'contrastiveVI', 'config': cfg, 'study': study, 'seed': seed, 'fit_seconds': fit_s,
           'canonical_log_mse': float(np.mean((pred_can - X[tes]) ** 2)), 'salient_active_units': int((var > 0.01).sum()),
           'salient_total_kl_nats': float(kl.sum()), 'salient_dims_kl_gt_0p01': int((kl > 0.01).sum()),
           'probe_within_study_salient_ungated': float(balanced_accuracy_score(yte, clf.predict(mte)))}
    enc = module.s_encoder; enc.eval()

    def fun(x):
        return enc(x, torch.zeros((x.shape[0], 1), dtype=torch.long))[0].square().sum(1)
    cells = z[idx_key]
    logC = np.log1p(C).astype(np.float32)
    controls = trs[treated[trs] == 0]
    igrows = []
    targets = {k: SETS[k] for k in set_names}
    for base_name, baseline in (('zero', np.zeros(C.shape[1], np.float32)), ('control_median', np.median(logC[controls], axis=0).astype(np.float32))):
        ig, info = integrated_gradients(fun, logC[cells], baseline)
        mabs = np.abs(ig).mean(0); rec, top = recovery(mabs, genes, targets)
        igrows.append(dict(model='contrastiveVI', config=cfg, seed=seed, group=label, study=study, baseline=base_name,
                           encoder_input='log1p counts', **rec, **info, top50_genes='|'.join(top)))
        pd.DataFrame({'gene': genes, 'mean_abs_ig': mabs}).to_csv(run / f'{name}_{base_name}_ig.csv.gz', index=False, compression='gzip')
    pd.DataFrame(igrows).to_csv(run / f'{name}_ig_recovery.csv', index=False)
    row['seconds_total'] = time.time() - t0
    done.write_text(json.dumps(row, indent=2) + '\n'); print(row, flush=True)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--config', required=True); ap.add_argument('--study', required=True)
    ap.add_argument('--seed', type=int, default=0); ap.add_argument('--smoke', action='store_true')
    a = ap.parse_args(strip_threads_arg(sys.argv[1:]))
    fit(a.config, a.study, a.seed, smoke=a.smoke)
