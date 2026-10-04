"""Phase 0/1: regenerate canonical benchmark inputs and compare Ved's archived files.

Regeneration reuses the logic of ``scripts/06_validation.py --prepare`` (via
``src.peerreview.core.prepare_canonical``) with outputs redirected to
``runs/peerreview_20261003/canonical_inputs``. The modeled genes are the committed
``benchmark_gene_universe.csv`` (the HVG reselection is recomputed and compared).

The archive in ``../ved_archive`` is opened read-only.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.peerreview.threads import configure  # noqa: E402  (must precede numpy)
N_THREADS = configure(3)

import json  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import torch  # noqa: E402
import anndata as ad  # noqa: E402

from src.peerreview import core  # noqa: E402

torch.set_num_threads(N_THREADS)
IN = core.PR_RUN / 'canonical_inputs'
OUT = core.PR_OUT
IN.mkdir(parents=True, exist_ok=True)


def regenerate():
    if (IN / 'inputs.npz').exists() and (IN / 'cells.csv').exists():
        print('canonical inputs already regenerated', flush=True)
        return
    a = ad.read_h5ad(core.RAW_H5AD)
    genes = pd.read_csv(core.CORR_OUT / 'benchmark_gene_universe.csv')['gene'].astype(str).tolist()
    p = core.prepare_canonical(a, genes=genes)
    np.savez_compressed(IN / 'inputs.npz', X=p['X'], d=p['d'], b=p['b'], train=p['train'], val=p['val'],
                        test=p['test'], genes=p['genes'])
    np.savez_compressed(IN / 'count_inputs.npz', C=p['C'], library=p['library'])
    p['obs'].to_csv(IN / 'cells.csv')
    hv = set(p['hvg'].tolist()); gu = set(genes)
    core.write_json(IN / 'hvg_reselection_check.json', {
        'n_hvg_reselected': len(hv), 'n_committed_universe': len(gu), 'n_intersection': len(hv & gu),
        'reselected_not_in_committed': sorted(hv - gu), 'committed_not_in_reselected': sorted(gu - hv),
        'order_identical': list(p['hvg']) == genes})
    print('regenerated', p['X'].shape, len(p['train']), len(p['val']), len(p['test']), flush=True)


def compare_inputs():
    v = np.load(core.VED / 'inputs.npz'); g = np.load(IN / 'inputs.npz')
    rows = []
    for k in sorted(set(v.files) | set(g.files)):
        if k not in v.files or k not in g.files:
            rows.append(dict(array=k, status='missing_in_' + ('ved' if k not in v.files else 'regenerated')))
            continue
        a, b = v[k], g[k]
        row = dict(array=k, ved_dtype=str(a.dtype), regen_dtype=str(b.dtype), ved_shape=str(a.shape),
                   regen_shape=str(b.shape))
        if a.shape != b.shape:
            row['status'] = 'shape_differs'
        elif a.dtype.kind in 'iuUSb':
            row['n_mismatch'] = int((a != b).sum()); row['status'] = 'identical' if row['n_mismatch'] == 0 else 'differs'
        else:
            diff = np.abs(a.astype(np.float64) - b.astype(np.float64))
            row['max_abs_diff'] = float(diff.max()); row['n_nonzero_diff'] = int((diff > 0).sum())
            row['status'] = 'identical' if row['max_abs_diff'] == 0 else ('within_1e-6' if row['max_abs_diff'] <= 1e-6 else 'differs')
        rows.append(row)
    return pd.DataFrame(rows)


def model_checks(inputs_path, label):
    z = np.load(inputs_path)
    X, d, b, tr, va, te = z['X'], z['d'], z['b'], z['train'], z['val'], z['test']
    m, c = core.load_model(core.VED / 'full_0.pt')
    lat = core.compute_latents(m, X, d, b)
    treated = te[d[te] > 0]
    vs = np.var(lat['shared'][treated], axis=0).sum(); vd = np.var(lat['drug'][treated], axis=0).sum()
    res = {'inputs': label, 'checkpoint_seed': c.get('seed'), 'checkpoint_mode': c.get('mode'),
           'best_epoch': c.get('best_epoch'), 'test_mse': float(lat['mse'][te].mean()),
           'test_plugin_nll': float(lat['nll'][te].mean()), 'shared_fraction': float(vs / (vs + vd)),
           'validation_mse': float(lat['mse'][va].mean())}
    res.update(core.json_probes(lat, b, d, tr, te))
    return res, lat, z


def main():
    regenerate()
    inp = compare_inputs(); inp.to_csv(OUT / 'A0_inputs_comparison.csv', index=False); print(inp.to_string(), flush=True)
    committed = json.loads((core.CORR_RUN / 'full_0_metrics.json').read_text())
    hist = pd.read_csv(core.CORR_RUN / 'full_0_history.csv')
    committed['validation_mse'] = float(hist.loc[hist.epoch == committed['best_epoch'], 'val_mse'].iloc[0])
    rows = []; lat_by = {}
    for label, path in [('ved_inputs', core.VED / 'inputs.npz'), ('regenerated_inputs', IN / 'inputs.npz')]:
        res, lat, z = model_checks(path, label); lat_by[label] = (lat, z)
        for k, cv in committed.items():
            if k in ('mode',):
                continue
            rv = res.get(k)
            rows.append(dict(inputs=label, metric=k, committed=cv, recomputed=rv,
                             difference=(rv - cv) if isinstance(cv, (int, float)) and rv is not None else None))
    met = pd.DataFrame(rows); met.to_csv(OUT / 'A0_full0_metrics_comparison.csv', index=False)
    pd.set_option('display.precision', 17); print(met.to_string(), flush=True)

    # Latents vs archived latents.
    vl = np.load(core.VED / 'full_0_latents.npz'); lrows = []
    for label, (lat, _) in lat_by.items():
        for k in vl.files:
            a = vl[k].astype(np.float64); b2 = lat[k].astype(np.float64)
            lrows.append(dict(inputs=label, array=k, shape_archived=str(vl[k].shape), shape_recomputed=str(lat[k].shape),
                              max_abs_diff=float(np.abs(a - b2).max()) if a.shape == b2.shape else None,
                              max_abs_value=float(np.abs(a).max())))
    lt = pd.DataFrame(lrows); lt.to_csv(OUT / 'A0_full0_latents_comparison.csv', index=False); print(lt.to_string(), flush=True)

    # Seed-0 latent_probes rows (11_diagnostics logic) vs committed.
    cells = pd.read_csv(IN / 'cells.csv', index_col=0)
    probe_rows = []
    committed_probes = pd.read_csv(core.CORR_OUT / 'latent_probes.csv'); committed_probes = committed_probes[committed_probes.seed == 0]
    for label, (lat, z) in list(lat_by.items()) + [('archived_latents_file', (dict(vl), lat_by['ved_inputs'][1]))]:
        rr = core.latent_probes_11(lat, z['b'], z['d'], z['train'], z['test'], cells.source_cell_type.to_numpy(), 0)
        for r in rr:
            cv = committed_probes[(committed_probes.representation == r['representation']) & (committed_probes.task == r['task'])].balanced_accuracy.iloc[0]
            probe_rows.append(dict(latents=label, representation=r['representation'], task=r['task'], committed=cv,
                                   recomputed=r['balanced_accuracy'], difference=r['balanced_accuracy'] - cv))
    pr = pd.DataFrame(probe_rows); pr.to_csv(OUT / 'A0_latent_probes_seed0_comparison.csv', index=False); print(pr.to_string(), flush=True)

    # cells.csv vs regenerated metadata.
    vc = pd.read_csv(core.VED / 'cells.csv', index_col=0); crow = []
    crow.append(dict(check='n_rows', ved=len(vc), regenerated=len(cells), equal=len(vc) == len(cells)))
    crow.append(dict(check='index_identical', equal=bool(vc.index.equals(cells.index))))
    crow.append(dict(check='columns', ved='|'.join(vc.columns), regenerated='|'.join(cells.columns),
                     equal=list(vc.columns) == list(cells.columns)))
    for col in vc.columns:
        if col not in cells.columns:
            continue
        a, b2 = vc[col], cells[col]
        if pd.api.types.is_numeric_dtype(a) and pd.api.types.is_numeric_dtype(b2):
            dmax = float(np.nanmax(np.abs(a.to_numpy(float) - b2.to_numpy(float))))
            crow.append(dict(check='column:' + col, max_abs_diff=dmax, equal=dmax <= 1e-9))
        else:
            n = int((a.astype(str).to_numpy() != b2.astype(str).to_numpy()).sum())
            crow.append(dict(check='column:' + col, n_mismatch=n, equal=n == 0))
    cc = pd.DataFrame(crow); cc.to_csv(OUT / 'A0_cells_comparison.csv', index=False); print(cc.to_string(), flush=True)


if __name__ == '__main__':
    main()
