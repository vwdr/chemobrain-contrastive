"""One resumable MC-ContrastiveVI fit: train, latents, metrics, latent usage, IG, done marker."""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

from src.peerreview import core


def run_fit(inputs: dict, seed: int, outdir, name: str, *, do_ig: bool = True, save_ig_arrays: bool = False,
            save_checkpoint: bool = True, save_latents: bool = True, epochs: int = 100, log=print,
            extra: dict | None = None):
    """Fit the canonical full model on ``inputs`` (keys X, d, b, train, val, test, genes).

    Writes ``{outdir}/{name}_*`` files and ``{outdir}/{name}_done.json`` last. If the marker
    exists the fit is skipped and the marker contents are returned.
    """
    outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    done = outdir / f'{name}_done.json'
    if done.exists():
        log(f'{name}: done marker exists, skipping')
        return json.loads(done.read_text())
    X, d, b = inputs['X'], inputs['d'], inputs['b']
    tr, va, te = inputs['train'], inputs['val'], inputs['test']
    genes = np.asarray(inputs['genes'])
    prefix = outdir / name
    m, hist, be, secs = core.train_full(X, d, b, tr, va, te, seed, prefix, epochs=epochs, log=log,
                                        save_checkpoint=save_checkpoint)
    lat = core.compute_latents(m, X, d, b)
    if save_latents:
        np.savez_compressed(str(prefix) + '_latents.npz', **lat)
    fsh = core.shared_fractions(lat, d, te)
    usage, dim_rows = core.latent_usage(lat, d, te)
    pd.DataFrame(dim_rows).to_csv(str(prefix) + '_latent_usage_per_dimension.csv', index=False)
    res = {'name': name, 'seed': seed, 'best_epoch': be, 'epochs_run': len(hist), 'fit_seconds': secs,
           'test_mse': float(lat['mse'][te].mean()), 'validation_mse': float(lat['mse'][va].mean()),
           'test_plugin_nll': float(lat['nll'][te].mean()), 'shared_fraction': fsh['pooled'],
           'shared_fraction_doxorubicin': fsh['doxorubicin'], 'shared_fraction_cisplatin': fsh['cisplatin']}
    res.update(usage)
    res.update(core.json_probes(lat, b, d, tr, te))
    if extra:
        res.update(extra)
    if do_ig:
        t0 = time.time()
        idx = core.attribution_indices(d, te)
        arrays, comp = core.integrated_gradients(m, X, d, b, tr, idx)
        rows = []
        for (axis, base, target), ig in arrays.items():
            for g, a_, b_ in zip(genes, ig.mean(0), np.abs(ig).mean(0)):
                rows.append({'axis': axis, 'baseline': base, 'target': target, 'gene': g, 'mean_ig': float(a_),
                             'mean_abs_ig': float(b_)})
        pd.DataFrame(rows).to_csv(str(prefix) + '_ig_summary.csv.gz', index=False, compression='gzip')
        pd.DataFrame(comp).to_csv(str(prefix) + '_ig_completeness.csv', index=False)
        if save_ig_arrays:
            np.savez_compressed(str(prefix) + '_ig_arrays.npz',
                                **{f'{a}__{bb}__{t}': v for (a, bb, t), v in arrays.items()},
                                **{f'indices__{k}': v for k, v in idx.items()})
        res['ig_seconds'] = time.time() - t0
        res['ig_max_quadrature_points'] = int(max(c['quadrature_points'] for c in comp))
        res['ig_max_median_scaled_error'] = float(max(c['median_scaled_error'] for c in comp))
    (outdir / f'{name}_metrics.json').write_text(json.dumps(res, indent=2) + '\n')
    done.write_text(json.dumps(res, indent=2) + '\n')
    log(f'{name}: done test_mse={res["test_mse"]:.6f} F_sh={res["shared_fraction"]:.4f} '
        f'fit={secs:.0f}s ig={res.get("ig_seconds", 0):.0f}s')
    return res


def load_ig_rankings(outdir, names, target='norm_squared'):
    """axis -> list of mean|IG| vectors over (fit name x baseline), plus the gene order."""
    rankings = {}; genes = None
    for nm in names:
        df = pd.read_csv(Path(outdir) / f'{nm}_ig_summary.csv.gz')
        df = df[df.target == target]
        for (axis, base), g in df.groupby(['axis', 'baseline'], sort=True):
            if genes is None:
                genes = g.gene.to_numpy()
            assert (g.gene.to_numpy() == genes).all()
            rankings.setdefault(core.AXIS_NAME[axis], []).append(g.mean_abs_ig.to_numpy())
    return rankings, genes
