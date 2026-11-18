"""Tiered, resumable orchestrator for the experiments/ scripts.

Order (most valuable first, so a run that is killed partway still leaves the best results done):
  tier 0  exp6_fsh_noise_null, exp5_auroc_rescore      (no training)
  tier 1  exp1_freebits_sweep (the headline), exp7_epoch_ablation (does the collapse persist at convergence?)
  tier 2  exp2_null_calibration
  tier 3  exp3_effect_grid
  tier 4  exp4_more_methods

Flags:
  --tiers 0,1,2    which tiers to run (default all); always executed in the order above
  --max-hours N    no new job starts after N hours; the current job is finished and written first
  --smoke          pass --smoke to every experiment (synthetic data, 2 epochs, a few jobs, outputs under */smoke)
  --resume         pass --resume to every experiment (skip jobs whose JSON exists). Without it an experiment
                   that already has results aborts instead of mixing runs, and is reported as 'aborted'.
  --threads=N      thread count (read before numpy/torch load, as in scripts/)
SIGTERM: the current job is finished and written, nothing further starts, a summary is printed and the exit code is 0.
Failed jobs never stop the run; they are listed in each experiment's errors.log and retried by --resume.

Experiments run in this process (one after another), each through its own main(argv). exp3 contains jobs for a
second environment (.venv-mgvi); in this one its multiGroupVI jobs are logged as FAIL and retried later with
`.venv-mgvi/bin/python experiments/exp3_effect_grid.py --method multigroupvi --resume`.

Usage:
  python experiments/run_all.py --tiers 0,1 --max-hours 12 --resume
  python experiments/run_all.py --smoke
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from experiments import common  # noqa: E402  (configures threads before numpy/torch)

import argparse  # noqa: E402
import importlib  # noqa: E402
import signal  # noqa: E402
import time  # noqa: E402
import traceback  # noqa: E402

from src.peerreview.threads import strip_threads_arg  # noqa: E402

TIERS = {
    0: ['exp6_fsh_noise_null', 'exp5_auroc_rescore'],
    1: ['exp1_freebits_sweep', 'exp7_epoch_ablation'],
    2: ['exp2_null_calibration'],
    3: ['exp3_effect_grid'],
    4: ['exp4_more_methods'],
}


def _on_sigterm(signum, frame):
    common.STOP = True
    print('SIGTERM received: finishing the current job, then stopping', flush=True)


def parse_tiers(text: str):
    tiers = sorted({int(t) for t in text.split(',') if t.strip() != ''})
    bad = [t for t in tiers if t not in TIERS]
    if bad:
        raise SystemExit(f'unknown tier(s) {bad}; valid tiers are {sorted(TIERS)}')
    return tiers


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--tiers', default=','.join(map(str, TIERS)), help='comma-separated tiers (default 0,1,2,3,4)')
    ap.add_argument('--max-hours', type=float, default=None)
    ap.add_argument('--smoke', action='store_true')
    ap.add_argument('--resume', action='store_true')
    args = ap.parse_args(strip_threads_arg(sys.argv[1:] if argv is None else argv))
    tiers = parse_tiers(args.tiers)

    signal.signal(signal.SIGTERM, _on_sigterm)
    t_start = time.time()
    if args.max_hours is not None:
        common.DEADLINE = t_start + args.max_hours * 3600.0
    child_argv = (['--smoke'] if args.smoke else []) + (['--resume'] if args.resume else [])

    rows = []
    for tier in tiers:
        for name in TIERS[tier]:
            row = dict(tier=tier, experiment=name, done=0, failed=0, skipped=0, pending=0, note='')
            rows.append(row)
            if common.should_stop():
                row['note'] = 'not started (stop requested)'
                continue
            print(f'\n===== tier {tier}: {name} =====', flush=True)
            try:
                stats = importlib.import_module(f'experiments.{name}').main(child_argv)
                row.update({k: stats[k] for k in ('done', 'failed', 'skipped', 'pending')})
                if stats['stopped']:
                    row['note'] = 'stopped early (signal or --max-hours)'
            except SystemExit as e:              # argument errors, missing prepared data, results exist without --resume
                row['note'] = f'aborted: {e.code}'
                print(f'{name} aborted: {e.code}', flush=True)
            except Exception as e:               # setup failure outside any job (e.g. import error); later tiers still run
                row['note'] = f'error: {type(e).__name__}: {e}'
                print(traceback.format_exc(), flush=True)

    elapsed = time.time() - t_start
    print('\n' + '=' * 78 + f'\nSummary ({elapsed / 3600:.2f} h elapsed' + (', smoke' if args.smoke else '') + ')')
    print(f'{"tier":<5}{"experiment":<26}{"done":>6}{"failed":>8}{"skipped":>9}{"not run":>9}  note')
    for r in rows:
        print(f'{r["tier"]:<5}{r["experiment"]:<26}{r["done"]:>6}{r["failed"]:>8}{r["skipped"]:>9}{r["pending"]:>9}  {r["note"]}')
    print('=' * 78, flush=True)
    return 0      # failures are reported, not fatal; a SIGTERM stop is a clean exit


if __name__ == '__main__':
    sys.exit(main())
