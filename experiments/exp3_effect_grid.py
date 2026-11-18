"""Exp 3: effect-size grid for contrastiveVI and multiGroupVI.

delta in {0.25, 0.5} x scenario {shared, specific, mixed} (responder fraction 1.0, the setting scripts/41-44 use for
their delta-1/2 grid) x seeds {0,1,2}. The fitting and scoring code is scripts/43_r4_contrastivevi.py (pairwise
contrastiveVI, one job per study) and scripts/42_r4_multigroupvi.py (one job per configuration), called as modules;
only their input directory / output directory are redirected. Injected inputs for the new delta values are built
with scripts/41_r4_prepare.py's own main().

The two methods live in different environments: contrastiveVI needs the main venv (scvi-tools 1.5.1), multiGroupVI
the .venv-mgvi one (scvi-tools 0.18.0). Run with --method contrastivevi in the first and --method multigroupvi in
the second; --method both (default) logs a clean FAIL for whichever package is missing and carries on, and
--resume later picks those jobs up.

Usage:
  python experiments/exp3_effect_grid.py --out runs/experiments/exp3 --method contrastivevi [--seeds 3] [--resume]
  .venv-mgvi/bin/python experiments/exp3_effect_grid.py --out runs/experiments/exp3 --method multigroupvi --resume
  python experiments/exp3_effect_grid.py --smoke
"""
import importlib.util
import json
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from experiments import common  # noqa: E402

import pandas as pd  # noqa: E402

DELTAS = [0.25, 0.5]
SCENARIOS = ['shared', 'specific', 'mixed']
STUDIES = ['GSE216146', 'GSE271055']        # contrastiveVI runs one pairwise model per study
R = common.R


def load_module(name, filename):
    spec = importlib.util.spec_from_file_location(name, R / 'scripts' / filename)
    mod = importlib.util.module_from_spec(spec)
    argv = list(sys.argv)
    try:
        spec.loader.exec_module(mod)
    finally:
        sys.argv = argv
    return mod


def main(argv=None):
    ap = common.make_parser(__doc__, default_seeds=3, default_out='runs/experiments/exp3')
    ap.add_argument('--method', choices=['contrastivevi', 'multigroupvi', 'both'], default='both')
    args = common.parse(ap, argv)
    H = common.Harness(); H.prepare(args)
    ss = H.ss
    methods = ['contrastivevi', 'multigroupvi'] if args.method == 'both' else [args.method]
    jobs = []
    for method in methods:
        for scen in SCENARIOS:
            for delta in DELTAS:
                cfg = ss.config_name(scen, delta, 1.0)
                for seed in range(args.seeds):
                    for study in (STUDIES if method == 'contrastivevi' else [None]):
                        jid = f'{method}_{cfg}' + (f'_{study}' if study else '') + f'_s{seed}'
                        jobs.append((jid, dict(method=method, cfg=cfg, study=study, seed=seed)))
    jobs = common.select_jobs(jobs, args)
    runner = common.Runner('exp3', args)
    inputs_dir = args.out / 'inputs'
    mods = {}

    def get_module(method):
        if method not in mods:
            if method == 'contrastivevi':
                mods[method] = load_module('r4_cvi', '43_r4_contrastivevi.py')
            else:
                try:
                    mods[method] = load_module('r4_mgvi', '42_r4_multigroupvi.py')
                except ImportError as e:
                    raise RuntimeError(f'multiGroupVI needs the .venv-mgvi environment (scvi-tools 0.18.0): {e}') from e
            m = mods[method]
            if args.smoke or not m.SETS:
                m.set_gene_sets(H.gene_sets())
        return mods[method]

    def fn(job_id, method, cfg, study, seed):
        mod = get_module(method)
        mod.INPUTS = H.ensure_inputs(cfg, args, inputs_dir)
        run = args.out / method
        mod.RUN = run
        scen, delta, rf = ss.parse_config(cfg)
        if method == 'contrastivevi':
            name = f'{cfg}_{study}_s{seed}'
            mod.fit(cfg, study, seed, smoke=args.smoke)
        else:
            name = f'{cfg}_s{seed}'
            mod.fit(cfg, seed, smoke=args.smoke)
        outdir = run / 'smoke' if args.smoke else run               # both scripts write to RUN/'smoke' in smoke mode
        done = outdir / f'{name}_done.json'
        if not done.exists():
            raise RuntimeError(f'{done} was not written')
        rec = dict(method=method, config=cfg, scenario=scen, delta=delta, responder_fraction=rf, study=study, seed=seed,
                   inputs_dir=str(mod.INPUTS), smoke=bool(args.smoke))
        rec['metrics'] = json.loads(done.read_text())
        rec['ig_recovery'] = pd.read_csv(outdir / f'{name}_ig_recovery.csv', keep_default_na=False, na_values=['']
                                         ).drop(columns=['top50_genes'], errors='ignore').to_dict('records')
        return rec

    failed = runner.run(jobs, fn, resume=args.resume)
    common.write_summary(runner, [j for j, _ in jobs])
    return runner.stats


if __name__ == '__main__':
    sys.exit(1 if main()['failed'] else 0)
