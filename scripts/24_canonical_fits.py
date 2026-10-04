"""Phase 2: fresh canonical full MC-ContrastiveVI fits (seeds 0-2) and diagnostics.

Usage:
  python scripts/24_canonical_fits.py --fit --seed 0 [--threads=2]
  python scripts/24_canonical_fits.py --analyze

Training/diagnostic logic is copied from 06/09/11/14 (see src/peerreview/core.py); all
outputs go to runs/peerreview_20261003/canonical and analysis/peerreview_20261003.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.peerreview.threads import configure, strip_threads_arg  # noqa: E402
N_THREADS = configure(2)

import argparse  # noqa: E402
import itertools  # noqa: E402
import json  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import torch  # noqa: E402

from src.peerreview import core, fitjob  # noqa: E402

torch.set_num_threads(N_THREADS)
IN = core.PR_RUN / 'canonical_inputs'
RUN = core.PR_RUN / 'canonical'
OUT = core.PR_OUT


def load_inputs():
    z = np.load(IN / 'inputs.npz')
    return {k: z[k] for k in z.files}


def fit(seed):
    fitjob.run_fit(load_inputs(), seed, RUN, f'full_{seed}', save_ig_arrays=True,
                   extra={'threads': N_THREADS})


def analyze():
    z = load_inputs(); d, b, tr, te = z['d'], z['b'], z['train'], z['test']
    cells = pd.read_csv(IN / 'cells.csv', index_col=0)
    g19 = core.load_script('gene_pattern_audit', '19_gene_pattern_audit.py')
    frozen = g19.CONSENSUS
    lats = {s: dict(np.load(RUN / f'full_{s}_latents.npz')) for s in (0, 1, 2)}
    fresh = {s: json.loads((RUN / f'full_{s}_metrics.json').read_text()) for s in (0, 1, 2)}
    pd.DataFrame([fresh[s] for s in (0, 1, 2)]).to_csv(OUT / 'A1_fresh_metrics.csv', index=False)

    # --- metric comparison with committed values
    rows = []
    vr = pd.read_csv(core.CORR_OUT / 'variance_robustness.csv')
    lu = pd.read_csv(core.CORR_OUT / 'latent_usage_committed_summary.csv')
    lu_map = {'shared_active_units_var_gt_0p01': 'shared_active_units', 'dox_active_units_var_gt_0p01': 'dox_active_units',
              'cis_active_units_var_gt_0p01': 'cis_active_units', 'shared_dims_kl_gt_0p01': 'shared_dims_kl_gt_0p01',
              'dox_dims_kl_gt_0p01': 'dox_dims_kl_gt_0p01', 'cis_dims_kl_gt_0p01': 'cis_dims_kl_gt_0p01',
              'shared_mean_total_kl_nats': 'shared_mean_total_kl_nats', 'dox_mean_total_kl_nats': 'dox_mean_total_kl_nats',
              'cis_mean_total_kl_nats': 'cis_mean_total_kl_nats', 'shared_var_sum': 'shared_var_sum',
              'dox_var_sum': 'dox_var_sum', 'cis_var_sum': 'cis_var_sum'}
    for s in (0, 1, 2):
        cm = json.loads((core.CORR_RUN / f'full_{s}_metrics.json').read_text())
        hist = pd.read_csv(core.CORR_RUN / f'full_{s}_history.csv')
        cm['validation_mse'] = float(hist.loc[hist.epoch == cm['best_epoch'], 'val_mse'].iloc[0])
        cm['shared_fraction_doxorubicin'] = float(vr[(vr.seed == s) & (vr.drug == 'doxorubicin')].fraction.iloc[0])
        cm['shared_fraction_cisplatin'] = float(vr[(vr.seed == s) & (vr.drug == 'cisplatin')].fraction.iloc[0])
        lrow = lu[lu.seed == s].iloc[0]
        for k, v in lu_map.items():
            cm[k] = float(lrow[v])
        for k, cv in cm.items():
            if k in ('mode', 'seed') or k not in fresh[s]:
                continue
            rows.append(dict(seed=s, metric=k, committed=cv, fresh=fresh[s][k], difference=fresh[s][k] - cv))
    # --- latent probes (11 logic)
    probes = []
    for s in (0, 1, 2):
        probes += core.latent_probes_11(lats[s], b, d, tr, te, cells.source_cell_type.to_numpy(), s)
    probes = pd.DataFrame(probes); probes.to_csv(OUT / 'A1_latent_probes_fresh.csv', index=False)
    cp = pd.read_csv(core.CORR_OUT / 'latent_probes.csv')
    for r in probes.itertuples():
        cv = cp[(cp.seed == r.seed) & (cp.representation == r.representation) & (cp.task == r.task)].balanced_accuracy.iloc[0]
        rows.append(dict(seed=r.seed, metric=f'probe:{r.representation}:{r.task}', committed=cv, fresh=r.balanced_accuracy,
                         difference=r.balanced_accuracy - cv))
    comp = pd.DataFrame(rows); comp.to_csv(OUT / 'A1_metric_comparison.csv', index=False)

    # --- per-dimension usage
    pdim = []
    for s in (0, 1, 2):
        x = pd.read_csv(RUN / f'full_{s}_latent_usage_per_dimension.csv'); x.insert(0, 'seed', s); pdim.append(x)
    pd.concat(pdim).to_csv(OUT / 'A1_latent_usage_per_dimension.csv', index=False)

    # --- cross-seed CKA
    cs = pd.read_csv(core.CORR_OUT / 'cross_seed_latent_stability.csv'); ck = []
    for s, t in itertools.combinations((0, 1, 2), 2):
        for key in ('bg', 'shared', 'drug'):
            cv = cs[(cs.seed_a == s) & (cs.seed_b == t) & (cs.space == key)].test_cka.iloc[0]
            fv = core.cka(lats[s][key][te], lats[t][key][te])
            ck.append(dict(seed_a=s, seed_b=t, space=key, committed_test_cka=cv, fresh_test_cka=fv))
    ck = pd.DataFrame(ck); ck.to_csv(OUT / 'A1_cross_seed_cka.csv', index=False)

    # --- consensus lists
    rankings, genes = fitjob.load_ig_rankings(RUN, [f'full_{s}' for s in (0, 1, 2)])
    assert (genes == z['genes']).all()
    cons, counts = core.consensus_lists(rankings, genes)
    crow = []
    for axis, gl in cons.items():
        for rank, g in enumerate(gl, 1):
            crow.append(dict(axis=axis, rank_by_mean_abs_ig=rank, gene=g,
                             n_top50_of_6=int(counts[axis][list(genes).index(g)]), in_frozen_list=g in frozen[axis]))
    pd.DataFrame(crow).to_csv(OUT / 'A1_fresh_consensus_genes.csv', index=False)
    jac = []
    for axis in ('shared', 'doxorubicin', 'cisplatin'):
        A, B = set(cons[axis]), set(frozen[axis])
        jac.append(dict(axis=axis, n_fresh=len(A), n_frozen=len(B), n_overlap=len(A & B), jaccard=core.jaccard(A, B),
                        overlap_genes='|'.join(sorted(A & B)), fresh_only='|'.join(sorted(A - B)),
                        frozen_only='|'.join(sorted(B - A))))
    jac = pd.DataFrame(jac); jac.to_csv(OUT / 'A1_consensus_jaccard.csv', index=False)

    # --- per-seed top-50 Jaccard between fresh and committed IG is not available (committed IG table is not
    # in the repository); per-ranking agreement among fresh rankings:
    rj = []
    for axis, vecs in rankings.items():
        labels = [f'seed{s}_{bse}' for s in (0, 1, 2) for bse in ('control_median', 'zero')]
        for (i, a), (j, c) in itertools.combinations(enumerate(vecs), 2):
            rj.append(dict(axis=axis, ranking_a=labels[i], ranking_b=labels[j],
                           top50_jaccard=core.jaccard(core.top_genes(a, genes), core.top_genes(c, genes))))
    pd.DataFrame(rj).to_csv(OUT / 'A1_fresh_ranking_top50_jaccard.csv', index=False)

    # --- fresh seed 0 vs ved full_0
    vm, _ = core.load_model(core.VED / 'full_0.pt')
    vlat = core.compute_latents(vm, z['X'], d, b)
    vrow = []
    vf = core.shared_fractions(vlat, d, te); vu, _ = core.latent_usage(vlat, d, te)
    ved_metrics = {'test_mse': float(vlat['mse'][te].mean()), 'validation_mse': float(vlat['mse'][z['val']].mean()),
                   'test_plugin_nll': float(vlat['nll'][te].mean()), 'shared_fraction': vf['pooled'],
                   'shared_fraction_doxorubicin': vf['doxorubicin'], 'shared_fraction_cisplatin': vf['cisplatin']}
    ved_metrics.update(vu); ved_metrics.update(core.json_probes(vlat, b, d, tr, te))
    for k, v in ved_metrics.items():
        if k in fresh[0]:
            vrow.append(dict(quantity=k, ved_full_0=v, fresh_seed0=fresh[0][k], difference=fresh[0][k] - v))
    for key in ('bg', 'shared', 'drug', 'raw_shared', 'raw_drug'):
        vrow.append(dict(quantity=f'test_cka:{key}', ved_full_0=None, fresh_seed0=None,
                         difference=None, cka=core.cka(lats[0][key][te], vlat[key][te])))
    vv = pd.DataFrame(vrow); vv.to_csv(OUT / 'A1_seed0_vs_ved.csv', index=False)
    core.write_json(OUT / 'A1_fresh_consensus.json', cons)
    pd.set_option('display.width', 250); pd.set_option('display.max_rows', 500)
    for nm, df in [('comparison', comp), ('cka', ck), ('jaccard', jac[['axis', 'n_fresh', 'n_frozen', 'n_overlap', 'jaccard']]), ('ved', vv)]:
        print('=====', nm); print(df.to_string(index=False), flush=True)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--fit', action='store_true'); ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--analyze', action='store_true'); a = ap.parse_args(strip_threads_arg(sys.argv[1:]))
    if a.fit:
        fit(a.seed)
    if a.analyze:
        analyze()
