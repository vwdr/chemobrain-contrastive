"""R4 step 0: injected raw counts for the baseline positive control (main venv).

For each fixed configuration, the injected counts are rebuilt with the Phase 4 code
(scripts/31_semisynthetic.build_inputs gives the programs; src.peerreview.semisynth.inject with
default_rng(20261003) reproduces the same thinning draws), and saved with the reduced library sizes,
the log-normalized inputs, labels, split and responder masks to runs/peerreview_20261004/r4/inputs_<config>.npz.
A check asserts that log-normalizing the saved counts reproduces build_inputs' X exactly.
"""
import importlib.util
import sys
from pathlib import Path

import numpy as np

R = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(R))
spec = importlib.util.spec_from_file_location('ss', R / 'scripts' / '31_semisynthetic.py')
ss = importlib.util.module_from_spec(spec); sys.argv = sys.argv[:1]; spec.loader.exec_module(ss)
from src.peerreview import semisynth, core  # noqa: E402

CONFIGS = ['null', 'shared_d1_rf1.0', 'specific_d1_rf1.0', 'mixed_d1_rf1.0', 'shared_d2_rf1.0', 'specific_d2_rf1.0', 'mixed_d2_rf1.0']
OUT = R / 'runs' / 'peerreview_20261004' / 'r4'


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    base = ss.load_base()
    idx = core.attribution_indices(base['d'], base['test'])
    for cfg in CONFIGS:
        inp, progs = ss.build_inputs(cfg, base)
        C = base['C']
        Cn = semisynth.inject(C, progs, np.random.default_rng(ss.SEED)) if progs else C
        Lnew = base['L'] - (C.astype(np.float64) - Cn.astype(np.float64)).sum(1)
        X = semisynth.log_normalize(Cn, Lnew)
        assert np.array_equal(X, inp['X']), cfg
        np.savez_compressed(OUT / f'inputs_{cfg}.npz', C=Cn.astype(np.float32), L=Lnew, X=X, d=base['d'], b=base['b'],
                            treated=base['treated'], train=base['train'], val=base['val'], test=base['test'], genes=base['genes'],
                            resp_1p0=base['resp_1p0'], idx_shared=idx['shared'], idx_drug0=idx['drug0'], idx_drug1=idx['drug1'])
        print('saved', cfg, flush=True)


if __name__ == '__main__':
    main()
