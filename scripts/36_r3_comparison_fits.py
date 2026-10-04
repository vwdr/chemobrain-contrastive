"""R3.1: comparison fits for the reconstruction benchmark on the canonical inputs.

Runs the unmodified functions of scripts/06_validation.py (modes no_hsic, no_gating, gaussian_vae; PCA(32))
and scripts/12_nb_sensitivity.py (negative-binomial sensitivity) with their module-level ``RUN`` (and ``O``)
globals redirected to runs/peerreview_20261004/comparison. Inputs are copies of the regenerated canonical
inputs from runs/peerreview_20261003/canonical_inputs (content-identical to the committed inputs, see A0).
Thread counts are those set by the original scripts (06: torch 2 threads; 12: torch 1 thread).

Usage:
  python scripts/36_r3_comparison_fits.py --job mode --mode no_hsic --seed 0
  python scripts/36_r3_comparison_fits.py --job nb --seed 0
  python scripts/36_r3_comparison_fits.py --job pca
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.peerreview.threads import configure, strip_threads_arg  # noqa: E402
configure(2)

import argparse  # noqa: E402
import json  # noqa: E402
import shutil  # noqa: E402
import time  # noqa: E402

from src.peerreview import core  # noqa: E402

RUN = core.R / 'runs' / 'peerreview_20261004' / 'comparison'
SRC = core.PR_RUN / 'canonical_inputs'


def stage_inputs():
    RUN.mkdir(parents=True, exist_ok=True)
    for f in ('inputs.npz', 'count_inputs.npz'):
        if not (RUN / f).exists():
            shutil.copy2(SRC / f, RUN / f)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--job', required=True); ap.add_argument('--mode', default='full')
    ap.add_argument('--seed', type=int, default=0); a = ap.parse_args(strip_threads_arg(sys.argv[1:]))
    stage_inputs()
    name = {'mode': f'{a.mode}_{a.seed}', 'nb': f'negative_binomial_{a.seed}', 'pca': 'pca'}[a.job]
    done = RUN / f'{name}_done.json'
    if done.exists():
        print('done marker exists', name); return
    t0 = time.time()
    if a.job in ('mode', 'pca'):
        m06 = core.load_script('validation_06', '06_validation.py')
        m06.RUN = RUN; m06.O = RUN
        if a.job == 'mode':
            m06.train(a.mode, a.seed, 100)
        else:
            m06.pca()
    else:
        m12 = core.load_script('nb_sensitivity_12', '12_nb_sensitivity.py')
        m12.RUN = RUN; m12.v.RUN = RUN; m12.v.O = RUN
        m12.train(a.seed)
    done.write_text(json.dumps({'name': name, 'seconds': time.time() - t0}, indent=2) + '\n')
    print('finished', name, time.time() - t0)


if __name__ == '__main__':
    main()
