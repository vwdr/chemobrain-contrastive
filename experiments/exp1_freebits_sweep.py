"""Exp 1: free-bits sweep. free_bits in {0, .05, .1, .25, .5, 1} x {null, shared_d1_rf1.0, shared_d2_rf1.0} x 3 seeds = 54 fits.

One JSON per fit under <out>/jobs/exp1/: free_bits, config, seed, shared_fraction, per-block active-unit counts,
per-dimension KL and posterior-mean variance, probe balanced accuracies, test MSE and, for the injected
configurations, top-50 IG recovery / AUROC / AP against the known injected set (null: no truth, empty list).
Training budget: up to 300 epochs with early-stopping patience 30 (the archived fits all stopped at the 100-epoch cap
while still improving). Every record has epochs_run, best_epoch and converged (early stopping fired before the cap);
collect.py reports any fit with converged = False loudly.
free_bits=0 is the canonical objective. Free bits floor the per-dimension KL (batch-mean over the cells of the
block) of the shared and drug blocks only.

Usage:
  python experiments/exp1_freebits_sweep.py --out runs/experiments/exp1 [--seeds 3] [--resume] [--threads=3]
  python experiments/exp1_freebits_sweep.py --smoke
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from experiments import common  # noqa: E402

FREE_BITS = [0.0, 0.05, 0.1, 0.25, 0.5, 1.0]
MAX_EPOCHS = 300      # cap raised from the archived 100, which every archived fit hit (best epoch 96-99)
PATIENCE = 30         # early stopping on validation MSE; converged = stopped before the cap
CONFIGS = ['null', 'shared_d1_rf1.0', 'shared_d2_rf1.0']


def main(argv=None):
    args = common.parse(common.make_parser(__doc__, default_seeds=3, default_out='runs/experiments/exp1'), argv)
    H = common.Harness(); base = H.prepare(args)
    jobs = [(f'fb{fb:g}_{cfg}_s{seed}', dict(fb=fb, cfg=cfg, seed=seed))
            for fb in FREE_BITS for cfg in CONFIGS for seed in range(args.seeds)]
    jobs = common.select_jobs(jobs, args)
    runner = common.Runner('exp1', args)

    def fn(job_id, fb, cfg, seed):
        return common.mccvi_fit_record(H, base, args, cfg=cfg, seed=seed, free_bits=fb,
                                       fit_dir=args.out / 'fits' / 'exp1', experiment='exp1',
                                       epochs=MAX_EPOCHS, patience=PATIENCE)

    failed = runner.run(jobs, fn, resume=args.resume)
    common.write_summary(runner, [j for j, _ in jobs])
    return runner.stats


if __name__ == '__main__':
    sys.exit(1 if main()['failed'] else 0)
