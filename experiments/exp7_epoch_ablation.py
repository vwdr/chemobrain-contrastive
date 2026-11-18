"""Exp 7: epoch ablation. Does the collapse of the treatment blocks persist once the model has converged?

All 75 archived semi-synthetic fits stopped at the 100-epoch cap with the best validation MSE at epochs 96-99, i.e.
still improving. exp7 fits the canonical model (free_bits = 0, config shared_d1_rf1.0) with epoch caps
{100, 200, 300, 500} x 3 seeds = 12 fits, early-stopping patience 30 as in exp1, and records at each budget: active
units (per block and total), shared KL per dimension, shared fraction, recovery (top-50 enrichment, AUROC, AP),
epochs_run, best_epoch and converged. Fits at different caps with the same seed follow the same training trajectory
until early stopping or their cap (same seed, same data order), so differences between budgets are the effect of
training longer, but they are separate fits and not bitwise nested.

<out>/jobs/exp7/*.json: one per fit. <out>/epoch_ablation.csv: per-fit table; <out>/epoch_ablation_summary.csv: mean and
sd over seeds per budget.

Usage:
  python experiments/exp7_epoch_ablation.py --out runs/experiments/exp7 [--seeds 3] [--resume] [--threads=3]
  python experiments/exp7_epoch_ablation.py --smoke
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from experiments import common  # noqa: E402

import pandas as pd  # noqa: E402

CONFIG = 'shared_d1_rf1.0'
BUDGETS = [100, 200, 300, 500]
PATIENCE = 30


def flat_row(rec):
    ig = [r for r in rec.get('ig_recovery', []) if r['baseline'] == common.BASELINE_MEAN and r['target_definition'] == 'with_shared'
          and r['axis'] in ('cisplatin_group', 'doxorubicin_group')]
    mean = lambda k: float(pd.Series([r[k] for r in ig], dtype=float).mean()) if ig else float('nan')
    enr = pd.Series([r['top_hits'] / r['random_expectation_hits'] for r in ig], dtype=float).mean() if ig else float('nan')
    return dict(max_epochs=rec['max_epochs'], seed=rec['seed'], epochs_run=rec['epochs_run'], best_epoch=rec['best_epoch'],
                converged=rec['converged'], shared_active_units=rec['shared_active_units_var_gt_0p01'],
                dox_active_units=rec['dox_active_units_var_gt_0p01'], cis_active_units=rec['cis_active_units_var_gt_0p01'],
                shared_kl_per_dim=rec['shared_kl_per_dim'], shared_fraction=rec['shared_fraction'], test_mse=rec['test_mse'],
                recovery_top50_enrichment=float(enr), recovery_auroc=mean('auroc'), recovery_ap=mean('average_precision'))


def main(argv=None):
    args = common.parse(common.make_parser(__doc__, default_seeds=3, default_out='runs/experiments/exp7'), argv)
    H = common.Harness(); base = H.prepare(args)
    jobs = [(f'{CONFIG}_e{e}_s{seed}', dict(budget=e, seed=seed)) for e in BUDGETS for seed in range(args.seeds)]
    jobs = common.select_jobs(jobs, args)
    runner = common.Runner('exp7', args)

    def fn(job_id, budget, seed):
        return common.mccvi_fit_record(H, base, args, cfg=CONFIG, seed=seed, free_bits=0.0, fit_dir=args.out / 'fits' / 'exp7',
                                       experiment='exp7', epochs=budget, patience=PATIENCE, tag_extra=f'_e{budget}')

    runner.run(jobs, fn, resume=args.resume)
    results = runner.load_results([j for j, _ in jobs])
    if results:
        df = pd.DataFrame([flat_row(r) for r in results]).sort_values(['max_epochs', 'seed'])
        df.to_csv(args.out / 'epoch_ablation.csv', index=False)
        num = [c for c in df.columns if c not in ('max_epochs', 'seed', 'converged')]
        summ = df.groupby('max_epochs').agg(n=('seed', 'count'), n_converged=('converged', 'sum'),
                                            **{f'{c}_mean': (c, 'mean') for c in num}, **{f'{c}_sd': (c, 'std') for c in num}).reset_index()
        summ.to_csv(args.out / 'epoch_ablation_summary.csv', index=False)
        print(summ[['max_epochs', 'n', 'n_converged', 'epochs_run_mean', 'shared_active_units_mean', 'shared_kl_per_dim_mean',
                    'shared_fraction_mean', 'recovery_auroc_mean', 'recovery_top50_enrichment_mean']].round(3).to_string(index=False), flush=True)
    return runner.stats


if __name__ == '__main__':
    sys.exit(1 if main()['failed'] else 0)
