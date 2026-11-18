"""Exp 5: AUROC / average-precision re-score of existing IG rankings. NO training, nothing is refitted.

For every ``*_ig_summary.csv.gz`` under runs/ (or --runs DIR ...) whose name parses as a semi-synthetic
configuration of scripts/31 (``<config>_s<seed>_ig_summary.csv.gz``; a ``_fb<x>`` tag from exp1 is accepted) the full
1,500-gene mean|IG| ranking of every axis x baseline is scored against the genes actually injected in that
configuration. Fits whose name carries no injected truth (real-data fits, anything unparsable) are logged as skipped.

Truth sets, from scripts/31 (injected sets S, A, B come from D1_gene_sets.csv; which sets exist is
SETS_BY_SCENARIO[scenario], so a 'specific' fit has no shared-axis truth and a 'shared' fit has no own-set truth for
the drug axes):
  own_set      shared axis vs S, cisplatin (drug1) vs A, doxorubicin (drug0) vs B   (scripts/31 ig_scores)
  with_shared  shared vs S, drug axes vs S u own set                                (scripts/44 recovery)

Writes <out>/auroc_per_fit.csv (default analysis/experiments/): one row per fit x axis x baseline x truth definition.

Usage:
  python experiments/exp5_auroc_rescore.py [--out analysis/experiments] [--runs runs ...] [--resume]
  python experiments/exp5_auroc_rescore.py --smoke
"""
import re
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from experiments import common  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

NAME = re.compile(r'^(?P<cfg>null|(?:shared|specific|mixed)_d[0-9.]+_rf[0-9.]+)(?P<tag>(?:_fb[0-9.]+)?)_s(?P<seed>\d+)_ig_summary\.csv\.gz$')


def make_smoke_runs(H, out: Path) -> Path:
    """A synthetic fit directory so the smoke run has something to score (random IG ranking, injected genes boosted)."""
    gs = H.gene_sets(); genes = [f'g{i}' for i in range(common.SMOKE_GENES)]
    rng = np.random.default_rng(0); d = out / 'smoke_runs' / 'semisynthetic' / 'mixed_d1_rf1.0'; d.mkdir(parents=True, exist_ok=True)
    inj = set(gs.gene)
    rows = []
    for axis in ('shared', 'drug0', 'drug1'):
        for base in ('zero', 'control_median'):
            for g in genes:
                v = abs(rng.normal()) + (1.0 if g in inj else 0.0)
                rows.append(dict(axis=axis, baseline=base, target='norm_squared', gene=g, mean_ig=v, mean_abs_ig=v))
    pd.DataFrame(rows).to_csv(d / 'mixed_d1_rf1.0_s0_ig_summary.csv.gz', index=False, compression='gzip')
    pd.DataFrame(rows).to_csv(d.parent / 'real_data_full_0_ig_summary.csv.gz', index=False, compression='gzip')   # no truth: skipped
    return out / 'smoke_runs'


def main(argv=None):
    ap = common.make_parser(__doc__, default_seeds=0, default_out='analysis/experiments')
    ap.add_argument('--runs', type=Path, nargs='+', default=[common.R / 'runs'], help='directories searched recursively')
    args = common.parse(ap, argv)
    H = common.Harness(); H.prepare(args); ss, core = H.ss, H.core
    gs = H.gene_sets()
    roots = [make_smoke_runs(H, args.out)] if args.smoke else args.runs
    files = sorted({p for r in roots for p in Path(r).rglob('*_ig_summary.csv.gz')})
    jobs = []
    for p in files:
        rel = p.relative_to(common.R) if p.is_relative_to(common.R) else p
        jobs.append((str(rel).replace('/', '__').replace('\\', '__').removesuffix('.csv.gz'), dict(path=p)))
    jobs = common.select_jobs(jobs, args) if args.max_jobs else jobs
    runner = common.Runner('exp5', args)
    print(f'{len(files)} ig_summary files found', flush=True)

    def fn(job_id, path):
        m = NAME.match(path.name)
        if not m:
            return dict(path=str(path), rows=[], _status='SKIP-no-truth', note='file name is not a semi-synthetic configuration')
        cfg, seed = m.group('cfg'), int(m.group('seed'))
        scen, delta, rf = ss.parse_config(cfg)
        rows = common.ig_summary_rows(path, cfg, gs, ss, core) if scen != 'null' else []
        for r in rows:
            r.update(path=str(path), fit=path.name.removesuffix('_ig_summary.csv.gz'), config=cfg, scenario=scen, delta=delta,
                     responder_fraction=rf, seed=seed, tag=m.group('tag').lstrip('_'))
        return dict(path=str(path), rows=rows, _status='OK' if rows else 'SKIP-no-truth',
                    note='' if rows else 'null configuration: no injected genes')

    failed = runner.run(jobs, fn, resume=args.resume)
    rows = [r for res in runner.load_results([j for j, _ in jobs]) for r in res.get('rows', [])]
    if rows:
        cols = ['fit', 'config', 'scenario', 'delta', 'responder_fraction', 'seed', 'tag', 'axis', 'baseline', 'target_definition',
                'n_genes', 'n_positive', 'auroc', 'average_precision', 'average_precision_baseline', 'top_n', 'top_hits',
                'random_expectation_hits', 'path']
        df = pd.DataFrame(rows)[cols]
        args.out.mkdir(parents=True, exist_ok=True)
        df.to_csv(args.out / 'auroc_per_fit.csv', index=False)
        print(f'wrote {args.out / "auroc_per_fit.csv"} ({len(df)} rows from {df.fit.nunique()} fits)', flush=True)
    else:
        print('no scorable fits found; nothing written', flush=True)
    return runner.stats


if __name__ == '__main__':
    sys.exit(1 if main()['failed'] else 0)
