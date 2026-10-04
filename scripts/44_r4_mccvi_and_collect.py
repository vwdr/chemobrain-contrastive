"""R4 steps 4 and 6: matched MC-ContrastiveVI attribution recovery (existing Phase 4 grid checkpoints, no refits) and
collection of all R4 outputs.

MC-ContrastiveVI treatment-related latent for group g: concatenation of the ungated shared posterior mean and the
matching drug-head posterior mean (cisplatin: drug head 1; doxorubicin: drug head 0); target = squared norm of the
concatenation; input = log-normalized expression (the encoders' input); baselines zero and control median (median of
the study's training pseudo-controls, as 09_interpretation.py uses study-matched controls for drug axes); 128 test
pseudo-treated cells per group (core.attribution_indices). Probes / active units / KL for MC-ContrastiveVI are taken
from existing files (R1 probes; Phase 4 per-fit metrics).

Usage:
  python scripts/44_r4_mccvi_and_collect.py --mccvi
  python scripts/44_r4_mccvi_and_collect.py --collect
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.peerreview.threads import configure  # noqa: E402
configure(3)

import argparse  # noqa: E402
import json  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import torch  # noqa: E402

from src.peerreview import core  # noqa: E402
from src.peerreview.ig_generic import integrated_gradients, recovery  # noqa: E402

torch.set_num_threads(3)
R = core.R
R4 = R / 'runs' / 'peerreview_20261004' / 'r4'; OUT = R / 'analysis' / 'peerreview_20261004'
CONFIGS = ['null', 'shared_d1_rf1.0', 'specific_d1_rf1.0', 'mixed_d1_rf1.0', 'shared_d2_rf1.0', 'specific_d2_rf1.0', 'mixed_d2_rf1.0']
GS = pd.read_csv(core.PR_OUT / 'D1_gene_sets.csv')
SETS = {s: set(GS[GS.set == s].gene) for s in ('S', 'A', 'B')}
GROUPS = (('cisplatin', 0, 1, 'idx_drug1', ('S', 'A')), ('doxorubicin', 1, 0, 'idx_drug0', ('S', 'B')))


def mccvi():
    out = R4 / 'mccvi'; out.mkdir(parents=True, exist_ok=True)
    rows = []
    for cfg in CONFIGS:
        z = np.load(R4 / f'inputs_{cfg}.npz'); X, d, b, tr, genes = z['X'], z['d'], z['b'], z['train'], z['genes']
        for seed in (0, 1, 2):
            m, _ = core.load_model(core.PR_RUN / 'semisynthetic' / cfg / f'{cfg}_s{seed}.pt')
            for label, study_b, k, idx_key, set_names in GROUPS:
                def fun(x, k=k):
                    sh = m.shared_head(m.shared_trunk(x))[0]; dr = m.drug_heads[k](m.drug_trunks[k](x))[0]
                    return torch.cat([sh, dr], 1).square().sum(1)
                controls = tr[(d[tr] == 0) & (b[tr] == study_b)]
                for base_name, baseline in (('zero', np.zeros(X.shape[1], np.float32)), ('control_median', np.median(X[controls], axis=0))):
                    ig, info = integrated_gradients(fun, X[z[idx_key]], baseline)
                    mabs = np.abs(ig).mean(0); rec, top = recovery(mabs, genes, {s: SETS[s] for s in set_names})
                    rows.append(dict(model='MC-ContrastiveVI', config=cfg, seed=seed, group=label, baseline=base_name,
                                     encoder_input='log-normalized expression', **rec, **info, top50_genes='|'.join(top)))
        print('done', cfg, flush=True)
    pd.DataFrame(rows).to_csv(out / 'mccvi_ig_recovery.csv', index=False)


def collect():
    rd = lambda p: pd.read_csv(p, keep_default_na=False, na_values=[''])   # keep the configuration name 'null' as a string
    ig = [rd(R4 / 'mccvi' / 'mccvi_ig_recovery.csv')]
    ig += [rd(p) for p in sorted((R4 / 'multigroupvi').glob('*_ig_recovery.csv'))]
    ig += [rd(p) for p in sorted((R4 / 'contrastivevi').glob('*_ig_recovery.csv'))]
    ig = pd.concat(ig, ignore_index=True); ig.to_csv(OUT / 'R4_ig_recovery_per_fit.csv', index=False)
    mg = pd.DataFrame([json.loads(p.read_text()) for p in sorted((R4 / 'multigroupvi').glob('*_done.json'))])
    cv = pd.DataFrame([json.loads(p.read_text()) for p in sorted((R4 / 'contrastivevi').glob('*_done.json'))])
    mg.to_csv(OUT / 'R4_multigroupvi_fits.csv', index=False); cv.to_csv(OUT / 'R4_contrastivevi_fits.csv', index=False)
    # long-format metrics per model x config x seed x group
    rows = []
    for r in mg.itertuples():
        for label in ('cisplatin', 'doxorubicin'):
            rows.append(dict(model='multiGroupVI', config=r.config, seed=r.seed, group=label, canonical_log_mse=r.canonical_log_mse,
                             probe_within_study=getattr(r, f'probe_within_{label}_study_private_ungated'),
                             active_units=getattr(r, f'{label}_private_active_units'), n_dims=10,
                             total_kl_nats=getattr(r, f'{label}_private_total_kl_nats')))
    for r in cv.itertuples():
        rows.append(dict(model='contrastiveVI', config=r.config, seed=r.seed, group={'GSE216146': 'cisplatin', 'GSE271055': 'doxorubicin'}[r.study],
                         canonical_log_mse=r.canonical_log_mse, probe_within_study=r.probe_within_study_salient_ungated,
                         active_units=r.salient_active_units, n_dims=8, total_kl_nats=r.salient_total_kl_nats))
    p1 = pd.read_csv(OUT / 'R1_probes_per_fit.csv', keep_default_na=False, na_values=[''])
    p1 = p1[p1.label == 'pseudo_treated_vs_pseudo_control']
    for cfg in CONFIGS:
        for seed in (0, 1, 2):
            met = json.loads((core.PR_RUN / 'semisynthetic' / cfg / f'{cfg}_s{seed}_metrics.json').read_text())
            for label, study in (('cisplatin', 'GSE216146'), ('doxorubicin', 'GSE271055')):
                pr = p1[(p1.config == cfg) & (p1.seed == seed) & (p1.study == study)].set_index('representation').balanced_accuracy
                short = 'cis' if label == 'cisplatin' else 'dox'
                rows.append(dict(model='MC-ContrastiveVI', config=cfg, seed=seed, group=label, canonical_log_mse=met['test_mse'],
                                 probe_within_study=pr['raw_drug'], probe_within_study_raw_shared=pr['raw_shared'],
                                 active_units=met['shared_active_units_var_gt_0p01'] + met[f'{short}_active_units_var_gt_0p01'], n_dims=12,
                                 total_kl_nats=met['shared_mean_total_kl_nats'] + met[f'{short}_mean_total_kl_nats']))
    met = pd.DataFrame(rows); met.to_csv(OUT / 'R4_metrics_per_fit.csv', index=False)
    igm = ig.groupby(['model', 'config', 'group', 'baseline'], sort=False).agg(
        top50_hits_union_mean=('top50_hits_union', 'mean'), top50_hits_union_min=('top50_hits_union', 'min'),
        top50_hits_union_max=('top50_hits_union', 'max'), top50_hits_S_mean=('top50_hits_S', 'mean'),
        top50_hits_A_mean=('top50_hits_A', 'mean'), top50_hits_B_mean=('top50_hits_B', 'mean'),
        precision_union_mean=('precision_union', 'mean'), recall_union_mean=('recall_union', 'mean'),
        random_expectation_hits=('random_expectation_hits', 'first'), max_quadrature=('quadrature_points', 'max'),
        max_median_scaled_error=('median_scaled_error', 'max'), n_fits=('seed', 'nunique')).reset_index()
    igm.to_csv(OUT / 'R4_ig_recovery_summary.csv', index=False)
    mm = met.groupby(['model', 'config', 'group'], sort=False).agg(
        n_fits=('seed', 'nunique'), canonical_log_mse_mean=('canonical_log_mse', 'mean'), probe_within_study_mean=('probe_within_study', 'mean'),
        probe_within_study_min=('probe_within_study', 'min'), probe_within_study_max=('probe_within_study', 'max'),
        active_units_mean=('active_units', 'mean'), n_dims=('n_dims', 'first'), total_kl_nats_mean=('total_kl_nats', 'mean')).reset_index()
    mm.to_csv(OUT / 'R4_metrics_summary.csv', index=False)
    comb = mm.merge(igm[igm.baseline == 'zero'].drop(columns=['baseline', 'n_fits']), on=['model', 'config', 'group'], how='left') \
             .merge(igm[igm.baseline == 'control_median'][['model', 'config', 'group', 'top50_hits_union_mean']]
                    .rename(columns={'top50_hits_union_mean': 'top50_hits_union_mean_control_median'}), on=['model', 'config', 'group'], how='left')
    comb.to_csv(OUT / 'R4_comparison_table.csv', index=False)
    pd.set_option('display.width', 250); pd.set_option('display.max_rows', 200)
    print(comb[['model', 'config', 'group', 'n_fits', 'probe_within_study_mean', 'active_units_mean', 'total_kl_nats_mean',
                'top50_hits_union_mean', 'top50_hits_union_mean_control_median', 'random_expectation_hits']].round(3).to_string(index=False))


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--mccvi', action='store_true'); ap.add_argument('--collect', action='store_true')
    a = ap.parse_args([x for x in sys.argv[1:] if not x.startswith('--threads=')])
    if a.mccvi:
        mccvi()
    if a.collect:
        collect()
