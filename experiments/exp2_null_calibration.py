"""Exp 2: null calibration. 30 uninjected null fits (seeds 0-29), default settings (free_bits = 0), same per-fit
records as exp1. One JSON per fit under <out>/jobs/exp2/.

Usage:
  python experiments/exp2_null_calibration.py --out runs/experiments/exp2 [--seeds 30] [--resume] [--threads=3]
  python experiments/exp2_null_calibration.py --smoke
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from experiments import common  # noqa: E402


def main(argv=None):
    args = common.parse(common.make_parser(__doc__, default_seeds=30, default_out='runs/experiments/exp2'), argv)
    H = common.Harness(); base = H.prepare(args)
    jobs = common.select_jobs([(f'null_s{seed}', dict(seed=seed)) for seed in range(args.seeds)], args)
    runner = common.Runner('exp2', args)

    def fn(job_id, seed):
        return common.mccvi_fit_record(H, base, args, cfg='null', seed=seed, free_bits=0.0,
                                       fit_dir=args.out / 'fits' / 'exp2', experiment='exp2')

    failed = runner.run(jobs, fn, resume=args.resume)
    common.write_summary(runner, [j for j, _ in jobs])
    return runner.stats


if __name__ == '__main__':
    sys.exit(1 if main()['failed'] else 0)
