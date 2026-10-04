"""Task C3: matched scVI / pairwise contrastiveVI baselines with outputs redirected, plus extensions.

Reuses scripts/20_peerreview_matched_baselines.py (imported; its helper functions canonical_split, probe,
linear_cka, modeled_gene_log, make_adata and its model settings) with outputs written to
runs/peerreview_20261003/baselines and analysis/peerreview_20261003. Additions:

* canonical-scale log-expression MSE (full-library normalization, as in the manuscript) for every baseline,
  next to the script's modeled-gene-scale MSE. For a scvi-tools prediction normalized to 1e4 within the 1,500
  modeled genes (rho * 1e4), the canonical-scale prediction is log1p(rho * L_model / L_full * 1e4), where
  L_model is the cell's observed count total over the modeled genes (the library size the model conditions on,
  use_observed_lib_size=True) and L_full is its total over all 18,271 genes;
* MC-ContrastiveVI (fresh Phase 2 fits) and PCA(32) on both scales. A canonical-scale prediction p is moved to
  the modeled-gene scale as log1p(e / sum(e) * 1e4) with e = max(expm1(p), 0). PCA(32) is additionally fitted
  directly on modeled-gene-scale training data;
* within-study treatment probes for scVI (cisplatin vs PBS within GSE216146; doxorubicin vs control within
  GSE271055; 11_diagnostics definition) in addition to the pooled three-class probe;
* latent-utilization diagnostics for contrastiveVI's salient space (Var[E q(s|x)] > 0.01 active units and
  per-dimension KL to N(0, I)) on held-out target (treated) and background (control) cells.

Usage:
  python scripts/29_baselines_extended.py --job scvi --seed 0 [--threads=3]
  python scripts/29_baselines_extended.py --job cvi --pair cisplatin --seed 0
  python scripts/29_baselines_extended.py --collect
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.peerreview.threads import configure, strip_threads_arg  # noqa: E402
N_THREADS = configure(3)

import argparse  # noqa: E402
import json  # noqa: E402
import time  # noqa: E402
from itertools import combinations  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import torch  # noqa: E402

from src.peerreview import core  # noqa: E402

torch.set_num_threads(N_THREADS)
B20 = core.load_script('matched_baselines_20', '20_peerreview_matched_baselines.py')
RUN = core.PR_RUN / 'baselines'; OUT = core.PR_OUT; IN = core.PR_RUN / 'canonical_inputs'
RUN.mkdir(parents=True, exist_ok=True)
PAIRS = {'cisplatin': ('GSE216146', 'cisplatin'), 'doxorubicin': ('GSE271055', 'doxorubicin')}


def per_dim_kl(m, v):
    return 0.5 * (m ** 2 + v - 1.0 - np.log(v))


def setup():
    import anndata as ad
    import scvi
    scvi.settings.num_threads = N_THREADS
    raw = ad.read_h5ad(core.RAW_H5AD)
    train, val, test = B20.canonical_split(raw.obs)
    z = np.load(IN / 'inputs.npz')
    assert (train == z['train']).all() and (test == z['test']).all() and (val == z['val']).all()
    genes = pd.read_csv(core.CORR_OUT / 'benchmark_gene_universe.csv')['gene'].astype(str).tolist()
    lib_full = np.load(IN / 'count_inputs.npz')['library']
    return raw, train, test, genes, z, lib_full


def canonical_from_rho(pred_1e4, counts, lib_full_cells):
    l_model = np.asarray(counts, dtype=np.float64).sum(1, keepdims=True)
    return np.log1p(np.asarray(pred_1e4, dtype=np.float64) * l_model / lib_full_cells[:, None]).astype(np.float32)


def job_scvi(seed):
    done = RUN / f'scvi_{seed}_done.json'
    if done.exists():
        return
    import scvi
    raw, train, test, genes, z, lib_full = setup()
    t0 = time.time()
    scvi.settings.seed = seed; torch.manual_seed(seed); np.random.seed(seed)
    atr = B20.make_adata(raw, train, genes); ate = B20.make_adata(raw, test, genes)
    scvi.model.SCVI.setup_anndata(atr, batch_key='study')
    model = scvi.model.SCVI(atr, n_hidden=128, n_latent=32, n_layers=2, dropout_rate=0.1, dispersion='gene',
                            gene_likelihood='nb', use_observed_lib_size=True)
    model.train(max_epochs=B20.EPOCHS, train_size=1.0, validation_size=0.0, early_stopping=False, accelerator='cpu',
                devices=1, batch_size=256, enable_progress_bar=False)
    ztr = model.get_latent_representation(atr); zte = model.get_latent_representation(ate)
    pred = model.get_normalized_expression(ate, library_size=1e4, n_samples=1, return_mean=True, return_numpy=True)
    obs_log = B20.modeled_gene_log(ate.X)
    pred_log = np.log1p(np.asarray(pred, dtype=np.float64)).astype(np.float32)
    canon_pred = canonical_from_rho(pred, ate.X, lib_full[test])
    row = {'method': 'scVI', 'study_or_pair': 'combined', 'seed': seed,
           'modeled_gene_log_mse': float(np.mean((pred_log - obs_log) ** 2)),
           'canonical_log_mse': float(np.mean((canon_pred - z['X'][test]) ** 2)),
           'condition_balanced_accuracy': B20.probe(ztr, zte, atr.obs['drug'].astype(str).to_numpy(), ate.obs['drug'].astype(str).to_numpy()),
           'study_balanced_accuracy': B20.probe(ztr, zte, atr.obs['study'].astype(str).to_numpy(), ate.obs['study'].astype(str).to_numpy())}
    # within-study treatment probes (11_diagnostics definition)
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import balanced_accuracy_score
    btr = z['b'][train]; bte = z['b'][test]; ytr = (z['d'][train] > 0).astype(int); yte = (z['d'][test] > 0).astype(int)
    for name, bv in [('cisplatin_within_GSE216146', 0), ('doxorubicin_within_GSE271055', 1)]:
        c = LogisticRegression(max_iter=800, class_weight='balanced').fit(ztr[btr == bv], ytr[btr == bv])
        row[f'treatment_probe_{name}_balanced_accuracy'] = float(balanced_accuracy_score(yte[bte == bv], c.predict(zte[bte == bv])))
    row['fit_seconds'] = time.time() - t0; row['threads'] = N_THREADS
    np.savez_compressed(RUN / f'scvi_{seed}_latents.npz', train=ztr, test=zte)
    done.write_text(json.dumps(row, indent=2) + '\n'); print(row, flush=True)


def job_cvi(pair, seed):
    done = RUN / f'cvi_{pair}_{seed}_done.json'
    if done.exists():
        return
    import scvi
    raw, train, test, genes, z, lib_full = setup()
    study, treated = PAIRS[pair]; t0 = time.time()
    tr_mask = (raw.obs.iloc[train].study.astype(str).to_numpy() == study) & np.isin(raw.obs.iloc[train].drug.astype(str).to_numpy(), ['control', treated])
    te_mask = (raw.obs.iloc[test].study.astype(str).to_numpy() == study) & np.isin(raw.obs.iloc[test].drug.astype(str).to_numpy(), ['control', treated])
    pair_train = train[tr_mask]; pair_test = test[te_mask]; assert len(pair_train) == 3000
    scvi.settings.seed = seed; torch.manual_seed(seed); np.random.seed(seed)
    atr = B20.make_adata(raw, pair_train, genes); ate = B20.make_adata(raw, pair_test, genes)
    scvi.external.ContrastiveVI.setup_anndata(atr)
    bidx = np.where(atr.obs['drug'].astype(str).to_numpy() == 'control')[0]
    tidx = np.where(atr.obs['drug'].astype(str).to_numpy() == treated)[0]
    model = scvi.external.ContrastiveVI(atr, n_hidden=128, n_background_latent=16, n_salient_latent=8, n_layers=2,
                                        dropout_rate=0.1, use_observed_lib_size=True, wasserstein_penalty=0)
    model.train(background_indices=bidx, target_indices=tidx, max_epochs=B20.EPOCHS, train_size=1.0, validation_size=0.0,
                early_stopping=False, accelerator='cpu', devices=1, batch_size=128, enable_progress_bar=False)
    zbg_tr = model.get_latent_representation(atr, representation_kind='background')
    zsal_tr = model.get_latent_representation(atr, representation_kind='salient')
    zbg_te = model.get_latent_representation(ate, representation_kind='background')
    zsal_te = model.get_latent_representation(ate, representation_kind='salient')
    parts = model.get_normalized_expression(ate, library_size=1e4, n_samples=1, return_mean=True, return_numpy=True)
    ytr = (atr.obs['drug'].astype(str).to_numpy() == treated).astype(int); yte = (ate.obs['drug'].astype(str).to_numpy() == treated).astype(int)
    pred = np.asarray(parts['background'], dtype=np.float64).copy(); sal = np.asarray(parts['salient'], dtype=np.float64)
    pred[yte == 1] = sal[yte == 1]
    obs_log = B20.modeled_gene_log(ate.X); pred_log = np.log1p(pred).astype(np.float32)
    canon_pred = canonical_from_rho(pred, ate.X, lib_full[pair_test])
    # salient posterior moments for utilization diagnostics
    module = model.module; module.eval()
    with torch.no_grad():
        out = module._generic_inference(x=torch.tensor(np.asarray(ate.X, dtype=np.float32)),
                                        batch_index=torch.zeros((len(ate), 1), dtype=torch.long))
    qs_m = out['qs_m'].cpu().numpy(); qs_v = out['qs_v'].cpu().numpy()
    assert np.allclose(qs_m, zsal_te, atol=1e-5)
    usage = {}
    dim_rows = []
    for grp, msk in [('target_test', yte == 1), ('background_test', yte == 0)]:
        kl = per_dim_kl(qs_m[msk].astype(np.float64), qs_v[msk].astype(np.float64)).mean(0); var = qs_m[msk].var(0)
        usage[f'salient_{grp}_active_units_var_gt_0p01'] = int((var > 0.01).sum())
        usage[f'salient_{grp}_dims_kl_gt_0p01'] = int((kl > 0.01).sum())
        usage[f'salient_{grp}_total_kl_nats'] = float(kl.sum())
        usage[f'salient_{grp}_posterior_mean_var_sum'] = float(var.sum())
        for j in range(len(kl)):
            dim_rows.append(dict(pair=pair, seed=seed, group=grp, dimension=j, n_cells=int(msk.sum()), mean_kl_nats=float(kl[j]),
                                 posterior_mean_variance=float(var[j]), active_unit_var_gt_0p01=bool(var[j] > 0.01)))
    pd.DataFrame(dim_rows).to_csv(RUN / f'cvi_{pair}_{seed}_salient_usage_per_dimension.csv', index=False)
    row = {'method': 'contrastiveVI', 'study_or_pair': pair, 'seed': seed, 'n_train': len(atr), 'n_test': len(ate),
           'modeled_gene_log_mse': float(np.mean((pred_log - obs_log) ** 2)),
           'canonical_log_mse': float(np.mean((canon_pred - z['X'][pair_test]) ** 2)),
           'background_condition_balanced_accuracy': B20.probe(zbg_tr, zbg_te, ytr, yte),
           'salient_condition_balanced_accuracy': B20.probe(zsal_tr, zsal_te, ytr, yte),
           'fit_seconds': time.time() - t0, 'threads': N_THREADS, **usage}
    np.savez_compressed(RUN / f'cvi_{pair}_{seed}_latents.npz', bg_train=zbg_tr, sal_train=zsal_tr, bg_test=zbg_te,
                        sal_test=zsal_te, qs_v_test=qs_v, test_indices=pair_test)
    done.write_text(json.dumps(row, indent=2) + '\n'); print(row, flush=True)


def canon_to_modeled(p):
    e = np.clip(np.expm1(np.asarray(p, dtype=np.float64)), 0, None)
    s = e.sum(1, keepdims=True); s[s <= 0] = 1.0
    return np.log1p(e / s * 1e4).astype(np.float32)


def mccvi_and_pca():
    from sklearn.decomposition import PCA
    z = np.load(IN / 'inputs.npz'); C = np.load(IN / 'count_inputs.npz')['C']
    X, d, b, tr, te = z['X'], z['d'], z['b'], z['train'], z['test']
    obs_mod = B20.modeled_gene_log(C[te])
    subsets = {'combined': np.ones(len(te), bool), 'cisplatin': (b[te] == 0), 'doxorubicin': (b[te] == 1)}
    rows = []
    for seed in (0, 1, 2):
        m, _ = core.load_model(core.PR_RUN / 'canonical' / f'full_{seed}.pt')
        with torch.no_grad():
            xt = torch.from_numpy(X[te]); dt = torch.from_numpy(d[te]); bt = torch.from_numpy(b[te])
            bg, sh, dr, _ = core.means(m, xt, dt); p = m.decode(bg, sh, dr, bt).numpy()
        pm = canon_to_modeled(p)
        for sub, msk in subsets.items():
            rows.append(dict(method='MC-ContrastiveVI (fresh)', study_or_pair=sub, seed=seed,
                             canonical_log_mse=float(np.mean((p[msk] - X[te][msk]) ** 2)),
                             modeled_gene_log_mse=float(np.mean((pm[msk] - obs_mod[msk]) ** 2))))
    pca = PCA(32, random_state=1729).fit(X[tr]); p = pca.inverse_transform(pca.transform(X[te])); pm = canon_to_modeled(p)
    Xm_tr = B20.modeled_gene_log(C[tr]); pca2 = PCA(32, random_state=1729).fit(Xm_tr)
    pm2 = pca2.inverse_transform(pca2.transform(obs_mod))
    for sub, msk in subsets.items():
        rows.append(dict(method='PCA(32) fit on canonical scale', study_or_pair=sub, seed=None,
                         canonical_log_mse=float(np.mean((p[msk] - X[te][msk]) ** 2)),
                         modeled_gene_log_mse=float(np.mean((pm[msk] - obs_mod[msk]) ** 2))))
        rows.append(dict(method='PCA(32) fit on modeled-gene scale', study_or_pair=sub, seed=None, canonical_log_mse=None,
                         modeled_gene_log_mse=float(np.mean((pm2[msk] - obs_mod[msk]) ** 2))))
    return pd.DataFrame(rows)


def collect():
    sc_rows = [json.loads((RUN / f'scvi_{s}_done.json').read_text()) for s in (0, 1, 2)]
    cv_rows = [json.loads((RUN / f'cvi_{p}_{s}_done.json').read_text()) for p in PAIRS for s in (0, 1, 2)]
    scvi_df = pd.DataFrame(sc_rows); cvi_df = pd.DataFrame(cv_rows)
    scvi_df.to_csv(OUT / 'C3_scvi_metrics.csv', index=False); cvi_df.to_csv(OUT / 'C3_contrastivevi_metrics.csv', index=False)
    z = np.load(IN / 'inputs.npz'); te = z['test']
    lat_s = {s: np.load(RUN / f'scvi_{s}_latents.npz')['test'] for s in (0, 1, 2)}
    lat_c = {(p, s): np.load(RUN / f'cvi_{p}_{s}_latents.npz') for p in PAIRS for s in (0, 1, 2)}
    ck = []
    mc = {s: np.load(core.PR_RUN / 'canonical' / f'full_{s}_latents.npz') for s in (0, 1, 2)}
    for a, c in combinations((0, 1, 2), 2):
        ck.append(dict(method='scVI', study_or_pair='combined', representation='latent', seed_a=a, seed_b=c,
                       linear_cka=B20.linear_cka(lat_s[a], lat_s[c])))
        for p in PAIRS:
            for rep, key in (('background', 'bg_test'), ('salient', 'sal_test')):
                ck.append(dict(method='contrastiveVI', study_or_pair=p, representation=rep, seed_a=a, seed_b=c,
                               linear_cka=B20.linear_cka(lat_c[(p, a)][key], lat_c[(p, c)][key])))
            # salient restricted to target cells, for comparability with gated MC-ContrastiveVI blocks
            ti = lat_c[(p, a)]['test_indices']; tmask = z['d'][ti] > 0
            ck.append(dict(method='contrastiveVI', study_or_pair=p, representation='salient_target_cells_only', seed_a=a, seed_b=c,
                           linear_cka=B20.linear_cka(lat_c[(p, a)]['sal_test'][tmask], lat_c[(p, c)]['sal_test'][tmask])))
        for key in ('bg', 'shared', 'drug'):
            ck.append(dict(method='MC-ContrastiveVI (fresh)', study_or_pair='combined', representation=key, seed_a=a, seed_b=c,
                           linear_cka=core.cka(mc[a][key][te], mc[c][key][te])))
    ck = pd.DataFrame(ck); ck.to_csv(OUT / 'C3_cross_seed_cka.csv', index=False)
    mp = mccvi_and_pca(); mp.to_csv(OUT / 'C3_mccvi_pca_both_scales.csv', index=False)
    # combined MSE table (both scales)
    allm = pd.concat([scvi_df[['method', 'study_or_pair', 'seed', 'canonical_log_mse', 'modeled_gene_log_mse']],
                      cvi_df[['method', 'study_or_pair', 'seed', 'canonical_log_mse', 'modeled_gene_log_mse']], mp], ignore_index=True)
    allm.to_csv(OUT / 'C3_mse_both_scales.csv', index=False)
    summ = (allm.groupby(['method', 'study_or_pair'], dropna=False)
            .agg(n=('modeled_gene_log_mse', 'size'), canonical_mean=('canonical_log_mse', 'mean'), canonical_sd=('canonical_log_mse', 'std'),
                 modeled_mean=('modeled_gene_log_mse', 'mean'), modeled_sd=('modeled_gene_log_mse', 'std')).reset_index())
    summ.to_csv(OUT / 'C3_mse_summary.csv', index=False)
    dims = pd.concat([pd.read_csv(RUN / f'cvi_{p}_{s}_salient_usage_per_dimension.csv') for p in PAIRS for s in (0, 1, 2)])
    dims.to_csv(OUT / 'C3_contrastivevi_salient_usage_per_dimension.csv', index=False)
    import scvi
    core.write_json(OUT / 'C3_design.json', {
        'scvi_tools_version': scvi.__version__, 'reused_script': 'scripts/20_peerreview_matched_baselines.py (imported, unmodified)',
        'seeds': [0, 1, 2], 'epochs': int(B20.EPOCHS), 'threads_per_fit': N_THREADS,
        'canonical_scale_conversion': 'log1p(rho*1e4 * L_model / L_full); L_model = observed modeled-gene count total, L_full = all-gene total',
        'modeled_scale_conversion_for_canonical_predictions': 'log1p(e/sum(e)*1e4), e = max(expm1(p), 0)',
        'within_study_probe': "LogisticRegression(max_iter=800, class_weight='balanced'), canonical train -> test cells of the study, y = treated",
        'salient_usage': 'qs_m/qs_v from module._generic_inference on held-out pair test cells; KL to N(0,I); active unit Var[qs_m] > 0.01'})
    pd.set_option('display.width', 250)
    print(scvi_df.to_string(index=False)); print(cvi_df.to_string(index=False)); print(summ.to_string(index=False))
    print(ck.groupby(['method', 'study_or_pair', 'representation']).linear_cka.agg(['mean', 'min', 'max']).to_string())


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--job'); ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--pair'); ap.add_argument('--collect', action='store_true')
    a = ap.parse_args(strip_threads_arg(sys.argv[1:]))
    if a.job == 'scvi':
        job_scvi(a.seed)
    elif a.job == 'cvi':
        job_cvi(a.pair, a.seed)
    if a.collect:
        collect()
