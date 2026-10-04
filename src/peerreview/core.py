"""Path-parameterized copies of the canonical MC-ContrastiveVI fitting and diagnostics.

The logic is copied (not modified) from the canonical scripts so that the new
analyses can write into ``runs/peerreview_20261003`` and
``analysis/peerreview_20261003`` without touching protected folders:

* ``scripts/06_validation.py``   -> canonical split, input preparation, ``full``-mode
  training loop, gated posterior means, metrics JSON probes;
* ``scripts/11_diagnostics.py``  -> latent probes, linear CKA, F_sh;
* ``scripts/14_latent_usage.py`` -> active units (Var[E q(z|x)] > 0.01) and per-dimension KL;
* ``scripts/09_interpretation.py`` -> attribution cell indices and adaptive
  Gauss-Legendre integrated gradients;
* ``scripts/19_gene_pattern_audit.py`` -> frozen consensus lists and enrichment
  (imported directly, see ``load_script``).
"""
from __future__ import annotations

import copy
import importlib.util
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score
from sklearn.model_selection import train_test_split

import src.models.mc_contrastive_vi as mm

R = Path(__file__).resolve().parents[2]
CORR_RUN = R / "runs" / "corrected_20260920"
CORR_OUT = R / "analysis" / "corrected_20260920"
PR_RUN = R / "runs" / "peerreview_20261003"
PR_OUT = R / "analysis" / "peerreview_20261003"
RAW_H5AD = R / "data" / "processed" / "recovered_counts.h5ad"
VED = R.parent / "ved_archive"
SPLIT_SEED = 1729
AU_VAR_THRESHOLD = 0.01
KL_SENSITIVITY_THRESHOLD = 0.01
DRUG_CODE = {"control": 0, "doxorubicin": 1, "cisplatin": 2}


def load_script(name: str, filename: str):
    """Import an existing script as a module without executing ``__main__``."""
    spec = importlib.util.spec_from_file_location(name, R / "scripts" / filename)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# --------------------------------------------------------------------- 06 logic
def fast_hsic(x, y):
    if len(x) < 4:
        return x.new_tensor(0.)
    K = mm._rbf_kernel(x)
    L = mm._rbf_kernel(y)
    K = K - K.mean(0, keepdim=True) - K.mean(1, keepdim=True) + K.mean()
    L = L - L.mean(0, keepdim=True) - L.mean(1, keepdim=True) + L.mean()
    return (K * L).sum() / (len(x) - 1) ** 2


# 06_validation.py replaces the module-level HSIC with the fast equivalent.
mm.hsic = fast_hsic


def canonical_split(strata: np.ndarray, eligible: np.ndarray, seed: int = SPLIT_SEED,
                    cap: int = 1500):
    """80/10/10 stratified split, then at most ``cap`` training cells per stratum."""
    tr, te = train_test_split(eligible, test_size=.2, stratify=strata[eligible], random_state=seed)
    va, te = train_test_split(te, test_size=.5, stratify=strata[te], random_state=seed)
    rng = np.random.default_rng(seed)
    tr = np.concatenate([rng.choice(tr[strata[tr] == s], min(cap, sum(strata[tr] == s)), replace=False)
                         for s in np.unique(strata[tr])])
    return tr, va, te


def prepare_canonical(a, genes=None, n_top_genes: int = 1500):
    """Copy of ``06_validation.prepare`` returning arrays instead of writing them.

    ``a`` is normalized in place (as in the original). If ``genes`` is given it
    is used as the modeled universe; the HVG selection is still computed and
    returned for comparison.
    """
    import scanpy as sc
    obs = a.obs.copy()
    pure = np.where(~obs.rescue.values.astype(bool))[0]
    strata = (obs.study.astype(str) + '_' + obs.drug.astype(str)).values
    tr, va, te = canonical_split(strata, pure)
    v = a[tr].copy()
    sc.pp.highly_variable_genes(v, flavor='seurat_v3', n_top_genes=n_top_genes, batch_key='study')
    hvg = v.var_names[v.var.highly_variable]
    if genes is None:
        genes = hvg
    genes = pd.Index(list(genes))
    counts = a[:, genes].X.toarray().astype('float32')
    library = np.asarray(a.X.sum(1)).ravel().astype('float32')
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    X = a[:, genes].X.toarray().astype('float32')
    d = obs.drug.map(DRUG_CODE).values.astype('int64')
    b = (obs.study == 'GSE271055').values.astype('int64')
    return dict(X=X, d=d, b=b, train=tr, val=va, test=te,
                genes=np.asarray(genes.tolist(), dtype='U'), C=counts, library=library,
                hvg=np.asarray(list(hvg), dtype='U'), obs=obs)


def means(m, x, d):
    """Gated posterior means (06_validation.means, mode 'full')."""
    e = m.encode(x, d)
    bg = e['mu_bg']
    sh = e['mu_shared'] * (d > 0)[:, None]
    dr = torch.cat([mu * (d == k + 1)[:, None] for k, mu in enumerate(e['mu_drug_list'])], 1)
    return bg, sh, dr, e


def make_model(n_genes: int):
    cfg = mm.MCContrastiveVIConfig(n_genes=n_genes, n_drugs=3, n_batches=2, hidden_dim=128, dropout=.1)
    return cfg, mm.MCContrastiveVI(cfg)


def load_model(path):
    c = torch.load(path, map_location='cpu', weights_only=True)
    m = mm.MCContrastiveVI(mm.MCContrastiveVIConfig(**c['config']))
    m.load_state_dict(c['state'])
    m.eval()
    return m, c


def train_full(X: np.ndarray, d: np.ndarray, b: np.ndarray, tr, va, te, seed: int, prefix,
               epochs: int = 100, log=print, save_checkpoint: bool = True):
    """Copy of ``06_validation.train`` for mode ``full`` with parameterized inputs/outputs.

    Returns (model, history DataFrame, best_epoch, seconds).
    """
    prefix = Path(prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    X = torch.from_numpy(np.ascontiguousarray(X, dtype=np.float32))
    d = torch.from_numpy(np.asarray(d, dtype=np.int64))
    b = torch.from_numpy(np.asarray(b, dtype=np.int64))
    torch.manual_seed(seed)
    rng = np.random.default_rng(seed)
    cfg, m = make_model(X.shape[1])
    opt = torch.optim.AdamW(m.parameters(), lr=.001, weight_decay=1e-5)
    best = float('inf'); bad = 0; history = []; start = time.time(); state = None; be = -1
    for epoch in range(epochs):
        m.train(); vals = []
        for ix in np.array_split(rng.permutation(tr), int(np.ceil(len(tr) / 256))):
            x, dd, bb = X[ix], d[ix], b[ix]
            e = m(x, dd, bb)
            loss = mm.mc_contrastive_loss(m, e, x, dd, 1., 10, 10, 5)['loss']
            opt.zero_grad(); loss.backward()
            torch.nn.utils.clip_grad_norm_(m.parameters(), 5); opt.step()
            vals.append(float(loss.detach()))
        m.eval(); sq = 0; n = 0
        with torch.no_grad():
            for ix in np.array_split(va, int(np.ceil(len(va) / 512))):
                bg, sh, dr, _ = means(m, X[ix], d[ix])
                p = m.decode(bg, sh, dr, b[ix])
                sq += float(((p - X[ix]) ** 2).sum()); n += X[ix].numel()
        vm = sq / n
        history.append({'epoch': epoch, 'train_objective': np.mean(vals), 'val_mse': vm})
        if vm < best - 1e-5:
            best = vm; state = copy.deepcopy(m.state_dict()); be = epoch; bad = 0
        else:
            bad += 1
        if epoch % 20 == 0:
            log(f'full seed={seed} epoch={epoch} val_mse={vm:.5f} seconds={time.time() - start:.0f}')
        if bad >= 15:
            break
    m.load_state_dict(state); m.eval()
    if save_checkpoint:
        torch.save({'state': state, 'config': cfg.__dict__, 'seed': seed, 'mode': 'full', 'best_epoch': be},
                   str(prefix) + '.pt')
    hist = pd.DataFrame(history)
    hist.to_csv(str(prefix) + '_history.csv', index=False)
    return m, hist, be, time.time() - start


def compute_latents(m, X: np.ndarray, d: np.ndarray, b: np.ndarray):
    """Latent arrays with the same keys as the canonical ``*_latents.npz``, plus log-variances."""
    Xt = torch.from_numpy(np.ascontiguousarray(X, dtype=np.float32))
    dt = torch.from_numpy(np.asarray(d, dtype=np.int64))
    bt = torch.from_numpy(np.asarray(b, dtype=np.int64))
    keys = ['bg', 'shared', 'drug', 'raw_shared', 'raw_drug', 'mse', 'nll', 'lv_shared', 'lv_drug']
    arrays = {k: [] for k in keys}
    m.eval()
    with torch.no_grad():
        for ix in np.array_split(np.arange(len(Xt)), int(np.ceil(len(Xt) / 512))):
            bg, sh, dr, e = means(m, Xt[ix], dt[ix])
            p = m.decode(bg, sh, dr, bt[ix])
            arrays['bg'].append(bg.numpy()); arrays['shared'].append(sh.numpy()); arrays['drug'].append(dr.numpy())
            arrays['raw_shared'].append(e['mu_shared'].numpy())
            arrays['raw_drug'].append(torch.cat(e['mu_drug_list'], 1).numpy())
            arrays['mse'].append(((p - Xt[ix]) ** 2).mean(1).numpy())
            arrays['nll'].append(mm.gaussian_nll(Xt[ix], p, m.log_sigma).numpy())
            arrays['lv_shared'].append(e['lv_shared'].numpy())
            arrays['lv_drug'].append(torch.cat(e['lv_drug_list'], 1).numpy())
    return {k: np.concatenate(v) for k, v in arrays.items()}


def variance_fraction(a, c):
    va = np.var(a, axis=0).sum(); vc = np.var(c, axis=0).sum()
    return float(va / (va + vc))


def shared_fractions(lat, d, te):
    """Pooled and per-drug F_sh on held-out treated test cells (gated posterior means)."""
    out = {}
    for name, k in [('pooled', None), ('doxorubicin', 1), ('cisplatin', 2)]:
        ii = te[d[te] > 0] if k is None else te[d[te] == k]
        out[name] = variance_fraction(lat['shared'][ii], lat['drug'][ii])
    return out


def per_dim_kl(mu, logvar):
    return 0.5 * (mu ** 2 + np.exp(logvar) - 1.0 - logvar)


def latent_usage(lat, d, te, latent_dim_drug: int = 4):
    """14_latent_usage.py logic on held-out treated test cells. Returns (summary, per-dim rows)."""
    treated = te[d[te] > 0]; dox = te[d[te] == 1]; cis = te[d[te] == 2]
    blocks = {
        ('shared', 'all_treated'): (lat['raw_shared'][treated], lat['lv_shared'][treated]),
        ('drug_specific', 'doxorubicin'): (lat['raw_drug'][dox, :latent_dim_drug], lat['lv_drug'][dox, :latent_dim_drug]),
        ('drug_specific', 'cisplatin'): (lat['raw_drug'][cis, latent_dim_drug:], lat['lv_drug'][cis, latent_dim_drug:]),
    }
    rows = []; summ = {}
    short = {('shared', 'all_treated'): 'shared', ('drug_specific', 'doxorubicin'): 'dox',
             ('drug_specific', 'cisplatin'): 'cis'}
    for key, (mu, lv) in blocks.items():
        kl = per_dim_kl(mu.astype(np.float64), lv.astype(np.float64)).mean(0)
        var = np.var(mu, axis=0)
        for j in range(mu.shape[1]):
            rows.append(dict(block=key[0], condition=key[1], dimension=j, n_cells=len(mu),
                             mean_kl_nats=float(kl[j]), posterior_mean_variance=float(var[j]),
                             active_unit_var_gt_0p01=bool(var[j] > AU_VAR_THRESHOLD),
                             kl_gt_0p01=bool(kl[j] > KL_SENSITIVITY_THRESHOLD)))
        s = short[key]
        summ[f'{s}_active_units_var_gt_0p01'] = int((var > AU_VAR_THRESHOLD).sum())
        summ[f'{s}_dims_kl_gt_0p01'] = int((kl > KL_SENSITIVITY_THRESHOLD).sum())
        summ[f'{s}_mean_total_kl_nats'] = float(kl.sum())
        summ[f'{s}_var_sum'] = float(var.sum())
        summ[f'{s}_max_dim_kl_nats'] = float(kl.max())
        summ[f'{s}_max_dim_posterior_mean_variance'] = float(var.max())
    return summ, rows


def json_probes(lat, b, d, tr, te):
    """Probes stored in the canonical metrics JSON (06_validation.train)."""
    out = {}
    for name, rep in [('bg', lat['bg']), ('shared_ungated', lat['raw_shared']), ('drug_gated', lat['drug']),
                      ('drug_ungated', lat['raw_drug'])]:
        for label, y in [('study', b), ('condition', d)]:
            clf = LogisticRegression(max_iter=500, class_weight='balanced').fit(rep[tr], y[tr])
            out[name + '_' + label + '_balanced_accuracy'] = float(balanced_accuracy_score(y[te], clf.predict(rep[te])))
    return out


def latent_probes_11(lat, b, d, tr, te, cell_type, seed):
    """Rows of latent_probes.csv (11_diagnostics.py logic) for one full-model seed."""
    rows = []
    aa = tr[b[tr] == 0]; bb = te[b[te] == 0]
    cell_type = np.asarray(cell_type)
    for label, rep in [('background', lat['bg']), ('shared_ungated', lat['raw_shared']), ('drug_ungated', lat['raw_drug'])]:
        c = LogisticRegression(max_iter=800, class_weight='balanced').fit(rep[aa], cell_type[aa])
        rows.append(dict(seed=seed, representation=label, task='cis_source_cell_type',
                         balanced_accuracy=balanced_accuracy_score(cell_type[bb], c.predict(rep[bb]))))
        for study, bv in [('doxorubicin', 1), ('cisplatin', 0)]:
            ti = tr[b[tr] == bv]; vi = te[b[te] == bv]; y = (d > 0).astype(int)
            c = LogisticRegression(max_iter=800, class_weight='balanced').fit(rep[ti], y[ti])
            rows.append(dict(seed=seed, representation=label, task=study + '_treatment',
                             balanced_accuracy=balanced_accuracy_score(y[vi], c.predict(rep[vi]))))
    return rows


def cka(x, y):
    x = x - x.mean(0); y = y - y.mean(0)
    return float(np.linalg.norm(x.T @ y, 'fro') ** 2 /
                 (np.linalg.norm(x.T @ x, 'fro') * np.linalg.norm(y.T @ y, 'fro') + 1e-20))


# --------------------------------------------------------------------- 09 logic
def attribution_indices(d, te, seed: int = 314):
    rng = np.random.default_rng(seed)
    return {'shared': np.concatenate([rng.choice(te[d[te] == dd], 64, replace=False) for dd in [1, 2]]),
            'drug0': rng.choice(te[d[te] == 1], 128, replace=False),
            'drug1': rng.choice(te[d[te] == 2], 128, replace=False)}


def integrated_gradients(m, X, d, b, tr, indices, baselines=('zero', 'control_median'),
                         targets=('norm_squared',)):
    """Adaptive Gauss-Legendre integrated gradients exactly as in 09_interpretation.py.

    Returns ({(axis, baseline, target): ig array}, completeness rows).
    """
    arrays = {}; comp = []
    m.eval()
    for axis, ix in indices.items():
        for base in baselines:
            controls = tr[d[tr] == 0]
            if axis != 'shared':
                controls = controls[b[controls] == (1 if axis == 'drug0' else 0)]
            baseline = np.zeros(X.shape[1], dtype=np.float32) if base == 'zero' else np.median(X[controls], axis=0)
            for target in targets:
                def fun(x):
                    if axis == 'shared':
                        mu = m.shared_head(m.shared_trunk(x))[0]
                    else:
                        k = int(axis[-1]); mu = m.drug_heads[k](m.drug_trunks[k](x))[0]
                    return mu.sum(1) if target == 'coordinate_sum' else mu.square().sum(1)
                inp = torch.tensor(X[ix]); ref = torch.tensor(baseline)[None, :]; delta = inp - ref
                with torch.no_grad():
                    diff = (fun(inp) - fun(ref)).numpy()
                for order in [256, 1024, 4096]:
                    nodes, weights = np.polynomial.legendre.leggauss(order); acc = torch.zeros_like(inp)
                    for node, weight in zip(nodes, weights):
                        xx = (ref + delta * float((node + 1) / 2)).detach().requires_grad_(True)
                        grad = torch.autograd.grad(fun(xx).sum(), xx)[0]; acc += grad * float(weight / 2)
                    ig = (delta * acc).detach().numpy(); err = np.abs(ig.sum(1) - diff)
                    if np.median(err) / (np.mean(np.abs(diff)) + 1e-8) < .01:
                        break
                comp.append({'axis': axis, 'baseline': base, 'target': target,
                             'median_absolute_error': float(np.median(err)), 'max_absolute_error': float(err.max()),
                             'median_scaled_error': float(np.median(err) / (np.mean(np.abs(diff)) + 1e-8)),
                             'quadrature_points': order, 'mean_abs_output_difference': float(np.mean(np.abs(diff)))})
                arrays[(axis, base, target)] = ig
    return arrays, comp


AXIS_NAME = {'shared': 'shared', 'drug0': 'doxorubicin', 'drug1': 'cisplatin'}


def consensus_lists(rankings, genes, top: int = 50, min_count: int = 4):
    """Paper rule: top ``top`` in at least ``min_count`` of the seed-by-baseline rankings.

    ``rankings`` maps axis -> list of mean-|IG| vectors (one per seed x baseline).
    Genes are ordered by the mean of the ranking vectors (descending).
    Returns {axis: [genes]} and {axis: counts vector}.
    """
    genes = np.asarray(genes)
    out = {}; cnts = {}
    for axis, vecs in rankings.items():
        V = np.vstack(vecs)
        counts = np.zeros(V.shape[1], dtype=int)
        for v in V:
            counts[np.argsort(v)[-top:]] += 1
        sel = np.where(counts >= min_count)[0]
        sel = sel[np.argsort(-V.mean(0)[sel], kind='stable')]
        out[axis] = genes[sel].tolist(); cnts[axis] = counts
    return out, cnts


def top_genes(vec, genes, top: int = 50):
    return set(np.asarray(genes)[np.argsort(vec)[-top:]].tolist())


def jaccard(a, b):
    a, b = set(a), set(b)
    return len(a & b) / len(a | b) if (a | b) else float('nan')


def write_json(path, obj):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(obj, indent=2, default=lambda o: o.item() if hasattr(o, 'item') else str(o)) + '\n')
