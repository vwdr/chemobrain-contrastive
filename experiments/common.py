"""Shared plumbing for the experiments/ scripts: CLI, resumable job runner, data access, scoring.

Every experiment script does ``from experiments import common`` first (this configures thread counts
before numpy/torch load, honouring ``--threads=N`` like the scripts/ do) and then:

    args = common.parse(common.make_parser(...))
    jobs = common.select_jobs([(job_id, kwargs), ...], args)
    common.Runner('exp1', args).run(jobs, fn, resume=args.resume)

Shared CLI: --out DIR, --seeds N, --smoke, --resume (plus --max-jobs N and --threads=N).
  * a job is skipped when ``<out>/jobs/<exp>/<job_id>.json`` exists and --resume is given; without --resume an
    existing result aborts the run so that two runs are never silently mixed;
  * every job runs inside try/except: the traceback goes to ``<out>/errors.log``, the run continues, and a
    failed job leaves no JSON (so --resume retries it);
  * ``<out>/progress.log`` gets one tab-separated line per job: timestamp, experiment, job id, status, seconds.
  * --smoke: synthetic data (200 cells, 50 genes, 2 epochs), outputs under ``<out>/smoke``, and only a few
    evenly spaced jobs unless --max-jobs is given. Nothing in smoke mode needs the recovered h5ad.
"""
from __future__ import annotations

import sys
from pathlib import Path

R = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(R))
from src.peerreview.threads import configure, strip_threads_arg  # noqa: E402

N_THREADS = configure(3)

import argparse  # noqa: E402
import datetime as dt  # noqa: E402
import importlib.util  # noqa: E402
import json  # noqa: E402
import math  # noqa: E402
import time  # noqa: E402
import traceback  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

# Cooperative stop for experiments/run_all.py: checked before every job, never mid-job.
STOP = False            # set by a SIGTERM handler
DEADLINE = None         # epoch seconds (time.time()) after which no new job is started


def should_stop() -> bool:
    return STOP or (DEADLINE is not None and time.time() >= DEADLINE)


SMOKE_CELLS, SMOKE_GENES, SMOKE_EPOCHS = 200, 50, 2
FULL_EPOCHS = 100            # run_fit / 06_validation default; every scripts/ fit uses it
SMOKE_JOBS = 3


# ------------------------------------------------------------------------------------------ CLI
def make_parser(description: str, default_seeds: int, default_out: str) -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=description)
    ap.add_argument('--out', type=Path, default=R / default_out, help=f'output directory (default {default_out})')
    ap.add_argument('--seeds', type=int, default=default_seeds, help=f'number of seeds, 0..N-1 (default {default_seeds})')
    ap.add_argument('--smoke', action='store_true', help='synthetic 200 cells x 50 genes, 2 epochs, a few jobs')
    ap.add_argument('--resume', action='store_true', help='skip jobs whose output JSON exists')
    ap.add_argument('--max-jobs', type=int, default=None, help='run only the first N jobs (smoke: N evenly spaced jobs)')
    return ap


def parse(ap: argparse.ArgumentParser, argv=None):
    args = ap.parse_args(strip_threads_arg(sys.argv[1:] if argv is None else argv))
    args.out = Path(args.out)
    if args.smoke:
        args.out = args.out / 'smoke'
    args.epochs = SMOKE_EPOCHS if args.smoke else FULL_EPOCHS
    return args


def select_jobs(jobs: list, args) -> list:
    n = args.max_jobs if args.max_jobs else (SMOKE_JOBS if args.smoke else None)
    if n is None or n >= len(jobs):
        return jobs
    if args.smoke and not args.max_jobs:       # spread over the grid so every code branch is touched
        idx = np.unique(np.linspace(0, len(jobs) - 1, n).round().astype(int))
        return [jobs[i] for i in idx]
    return jobs[:n]


# ------------------------------------------------------------------------------------------ runner
def _clean(o):
    """Make an object JSON-safe: numpy scalars/arrays to python, NaN/inf to null."""
    if isinstance(o, dict):
        return {str(k): _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple, set)):
        return [_clean(v) for v in o]
    if isinstance(o, np.ndarray):
        return _clean(o.tolist())
    if isinstance(o, np.generic):
        return _clean(o.item())
    if isinstance(o, float) and not math.isfinite(o):
        return None
    if isinstance(o, Path):
        return str(o)
    return o


def write_json(path, obj):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(_clean(obj), indent=2) + '\n')
    tmp.replace(path)


class Runner:
    def __init__(self, exp: str, args):
        self.exp = exp
        self.out = Path(args.out)
        self.jobs_dir = self.out / 'jobs' / exp
        self.jobs_dir.mkdir(parents=True, exist_ok=True)
        self.progress = self.out / 'progress.log'
        self.errors = self.out / 'errors.log'
        self.stats = dict(done=0, failed=0, skipped=0, pending=0, stopped=False)

    def result_path(self, job_id: str) -> Path:
        return self.jobs_dir / f'{job_id}.json'

    def _log(self, job_id, status, secs):
        with open(self.progress, 'a') as f:
            f.write(f'{dt.datetime.now().isoformat(timespec="seconds")}\t{self.exp}\t{job_id}\t{status}\t{secs:.1f}\n')

    def run(self, jobs, fn, resume: bool = False) -> int:
        """Run ``fn(job_id, **kwargs) -> dict`` for each (job_id, kwargs); returns the number of failures.

        ``fn`` may put '_status' in its record (popped before writing) to log something other than OK.
        """
        existing = [j for j, _ in jobs if self.result_path(j).exists()]
        if existing and not resume:
            raise SystemExit(f'{len(existing)} of {len(jobs)} job results already exist in {self.jobs_dir}; '
                             f'pass --resume to skip them (or choose another --out).')
        failed = 0
        for k, (job_id, kw) in enumerate(jobs, 1):
            p = self.result_path(job_id)
            if should_stop():
                self.stats['stopped'] = True
                self.stats['pending'] = len(jobs) - k + 1
                self._log(job_id, 'STOPPED', 0.0)
                print(f'stop requested (signal or --max-hours): {self.stats["pending"]} job(s) not started', flush=True)
                break
            if p.exists():
                self.stats['skipped'] += 1
                self._log(job_id, 'SKIP', 0.0); print(f'[{k}/{len(jobs)}] {job_id}: exists, skipped', flush=True); continue
            t0 = time.time()
            try:
                rec = fn(job_id, **kw)
                status = rec.pop('_status', 'OK') if isinstance(rec, dict) else 'OK'
                write_json(p, {'experiment': self.exp, 'job_id': job_id, **(rec or {})})
            except Exception:
                status = 'FAIL'; failed += 1
                with open(self.errors, 'a') as f:
                    f.write(f'--- {dt.datetime.now().isoformat(timespec="seconds")} {self.exp} {job_id}\n{traceback.format_exc()}\n')
            secs = time.time() - t0
            self.stats['failed' if status == 'FAIL' else 'skipped' if status.startswith('SKIP') else 'done'] += 1
            self._log(job_id, status, secs)
            print(f'[{k}/{len(jobs)}] {job_id}: {status} {secs:.1f}s', flush=True)
        if failed:
            print(f'{failed} job(s) failed; see {self.errors}', flush=True)
        return failed

    def load_results(self, job_ids=None) -> list:
        paths = [self.result_path(j) for j in job_ids] if job_ids is not None else sorted(self.jobs_dir.glob('*.json'))
        return [json.loads(p.read_text()) for p in paths if p.exists()]


def write_summary(runner: Runner, job_ids, name='summary.csv'):
    """One row per finished job with its scalar fields (nested lists/dicts are left in the JSON)."""
    rows = [{k: v for k, v in r.items() if not isinstance(v, (list, dict))} for r in runner.load_results(job_ids)]
    if rows:
        pd.DataFrame(rows).to_csv(runner.out / name, index=False)
    return len(rows)


# ------------------------------------------------------------------------------------------ harness
class Harness:
    """Handles on the existing code: scripts/41 (input builder), scripts/31 (semi-synthetic design), core, fitjob."""

    def __init__(self):
        argv = list(sys.argv)
        spec = importlib.util.spec_from_file_location('r4_prepare', R / 'scripts' / '41_r4_prepare.py')
        self.prep = importlib.util.module_from_spec(spec)
        try:
            spec.loader.exec_module(self.prep)      # 41 loads scripts/31 as .ss and trims sys.argv
        finally:
            sys.argv = argv
        self.ss = self.prep.ss
        from src.peerreview import core, fitjob      # noqa: E402
        import torch
        torch.set_num_threads(N_THREADS)
        self.core, self.fitjob = core, fitjob
        self.smoke_dir = None

    # -------------------------------------------------- data
    def prepare(self, args):
        """Returns the base dict of scripts/31. Smoke: a synthetic base and gene sets patched into the harness."""
        if args.smoke:
            return self._patch_smoke(Path(args.out))
        if not (self.ss.RUN / 'base.npz').exists():
            raise SystemExit(f'{self.ss.RUN / "base.npz"} missing: run `python scripts/31_semisynthetic.py --prepare` '
                             f'(needs data/processed/recovered_counts.h5ad) first.')
        return self.ss.load_base()

    def gene_sets(self) -> pd.DataFrame:
        p = Path(self.ss.OUT) / 'D1_gene_sets.csv'
        if not p.exists():
            raise SystemExit(f'{p} missing: run `python scripts/31_semisynthetic.py --prepare` first.')
        return pd.read_csv(p)

    def _patch_smoke(self, out_dir: Path):
        base, gs = synthetic_base(SMOKE_CELLS, SMOKE_GENES)
        self.smoke_dir = out_dir / 'smoke_data'; self.smoke_dir.mkdir(parents=True, exist_ok=True)
        gs.to_csv(self.smoke_dir / 'D1_gene_sets.csv', index=False)
        ss = self.ss
        ss.OUT = self.smoke_dir                      # build_inputs / true_shared_fraction read gene sets from ss.OUT
        ss.load_base = lambda: base                  # scripts/41 main() reads the base through ss.load_base
        # attribution_indices draws 64/128/128 cells; the smoke test set has only a handful per drug
        def small_indices(d, te, seed: int = 314, n: int = 4):
            rng = np.random.default_rng(seed)
            pick = lambda dd, k: rng.choice(te[d[te] == dd], min(k, int((d[te] == dd).sum())), replace=False)
            return {'shared': np.concatenate([pick(1, n), pick(2, n)]), 'drug0': pick(1, n), 'drug1': pick(2, n)}
        self.core.attribution_indices = small_indices
        # integrated gradients: 16 quadrature nodes instead of 256-4096 (real code path, coarse integral)
        real = np.polynomial.legendre.leggauss
        np.polynomial.legendre.leggauss = lambda n: real(min(int(n), 16))
        return base

    def ensure_inputs(self, cfg: str, args, run_inputs_dir: Path) -> Path:
        """Directory containing inputs_<cfg>.npz (injected counts + normalized X + split + IG cells).

        Reuses runs/peerreview_20261004/r4/ if it already holds the configuration (never in smoke), else builds
        the file with scripts/41's own main() into ``run_inputs_dir``.
        """
        if not args.smoke and (self.prep.OUT / f'inputs_{cfg}.npz').exists():
            return self.prep.OUT
        if not (run_inputs_dir / f'inputs_{cfg}.npz').exists():
            orig_out, orig_cfgs = self.prep.OUT, self.prep.CONFIGS
            self.prep.CONFIGS, self.prep.OUT = [cfg], run_inputs_dir
            try:
                self.prep.main()
            finally:
                self.prep.OUT, self.prep.CONFIGS = orig_out, orig_cfgs
        return run_inputs_dir


def synthetic_base(n_cells: int, n_genes: int, seed: int = 0):
    """Synthetic stand-in for runs/.../semisynthetic/base.npz (same keys) plus gene sets S, A, B of 8 genes each."""
    from src.peerreview import core
    rng = np.random.default_rng(seed)
    lam = np.exp(rng.normal(0.5, 0.8, n_genes))
    C = rng.poisson(lam, size=(n_cells, n_genes)).astype(np.float32) * 3     # enough counts to thin
    L = C.sum(1).astype(np.float64) * 4.0                                      # library over all genes
    b = (np.arange(n_cells) % 2).astype(np.int64)
    treated = np.zeros(n_cells, bool)
    for sb in (0, 1):
        idx = np.where(b == sb)[0]; treated[idx[: len(idx) // 2]] = True
    d = np.where(treated, np.where(b == 1, 1, 2), 0).astype(np.int64)
    study = np.where(b == 1, 'GSE271055', 'GSE216146')
    strata = study + '_' + np.where(treated, 'treated', 'control')
    tr, va, te = core.canonical_split(strata, np.arange(n_cells))
    resp_03 = np.zeros(n_cells, bool)
    for sb in (0, 1):
        pt = np.where(treated & (b == sb))[0]; resp_03[rng.choice(pt, int(round(0.3 * len(pt))), replace=False)] = True
    genes = np.asarray([f'g{i}' for i in range(n_genes)], dtype='U')
    base = dict(C=C, L=L, base_index=np.arange(n_cells), d=d, b=b, treated=treated, train=tr, val=va, test=te,
                genes=genes, resp_1p0=treated.copy(), resp_0p3=resp_03)
    rows = []; perm = rng.permutation(n_genes)
    for si, s in enumerate(('S', 'A', 'B')):
        for j, g in enumerate(perm[si * 8:(si + 1) * 8]):
            rows.append(dict(set=s, gene=genes[g], gene_index=int(g), decile=int(j // 2), direction='up' if j % 2 == 0 else 'down',
                             mean_log_norm_expression_base=0.0))
    return base, pd.DataFrame(rows)


# ------------------------------------------------------------------------------------------ scoring
def rank_metrics(score, genes, positives, top: int = 50) -> dict:
    """Top-``top`` hits, AUROC and average precision of ``score`` (higher = more important) for ``positives``."""
    from sklearn.metrics import average_precision_score, roc_auc_score
    genes = np.asarray(genes); score = np.asarray(score, dtype=float)
    y = np.isin(genes, list(positives)).astype(int)
    n_pos = int(y.sum())
    out = dict(n_genes=int(len(genes)), n_positive=n_pos, top_n=top,
               top_hits=int(y[np.argsort(score)[-top:]].sum()), random_expectation_hits=top * n_pos / max(len(genes), 1))
    if 0 < n_pos < len(y):
        out.update(auroc=float(roc_auc_score(y, score)), average_precision=float(average_precision_score(y, score)),
                   average_precision_baseline=n_pos / len(y))
    else:
        out.update(auroc=float('nan'), average_precision=float('nan'), average_precision_baseline=float('nan'))
    return out


def injected_targets(cfg: str, gs: pd.DataFrame, ss) -> dict:
    """Per attribution axis, the true injected gene set(s) of configuration ``cfg``.

    'own_set' follows scripts/31 ig_scores (shared axis vs S, cisplatin vs A, doxorubicin vs B); 'with_shared'
    follows scripts/44 (drug axis vs S u own set). Only sets actually injected in the scenario count, so a
    scenario that injects nothing on an axis gives no entry for it (not an empty-set score).
    """
    scen, _, _ = ss.parse_config(cfg)
    present = set(ss.SETS_BY_SCENARIO[scen])
    sets = {s: set(gs[gs.set == s].gene) for s in ('S', 'A', 'B')}
    own = {'shared': 'S', 'cisplatin': 'A', 'doxorubicin': 'B'}
    out = {'own_set': {}, 'with_shared': {}}
    for axis, s in own.items():
        if s in present:
            out['own_set'][axis] = sets[s]
        union = set().union(*[sets[x] for x in ({s} | {'S'}) if x in present])
        if union:
            out['with_shared'][axis] = union
    return out


BASELINE_MEAN = 'mean_of_baselines'


def ig_summary_rows(path, cfg: str, gs: pd.DataFrame, ss, core, top: int = 50) -> list:
    """AUROC / AP / top-``top`` hits of the full IG ranking of one fit, per axis x truth definition.

    Rows with baseline == 'mean_of_baselines' are the headline: the mean|IG| vectors of the two IG baselines are
    averaged gene by gene *before* scoring, and the average is scored once. Per-baseline rows are kept alongside.
    Extra rows axis == '<group>_group' (cisplatin_group, doxorubicin_group) score the study-group readout, the mean of the
    shared axis and that group's own drug axis, against S u own set (the readout compared across methods).
    """
    df = pd.read_csv(path)                                           # pandas infers the gzip from the suffix
    if 'target' in df.columns:
        df = df[df.target == 'norm_squared']
    targets = injected_targets(cfg, gs, ss)
    rows = []; axis_mean = {}
    for axis, g in df.groupby('axis', sort=True):
        name = core.AXIS_NAME[axis]
        wide = g.pivot_table(index='gene', columns='baseline', values='mean_abs_ig', aggfunc='first')
        vectors = {b: wide[b] for b in wide.columns}
        vectors[BASELINE_MEAN] = wide.mean(axis=1)
        axis_mean[name] = vectors[BASELINE_MEAN]
        for baseline, vec in vectors.items():
            for tdef, per_axis in targets.items():
                if name in per_axis:
                    rows.append(dict(axis=name, baseline=baseline, target_definition=tdef,
                                     **rank_metrics(vec.to_numpy(), vec.index.to_numpy(), per_axis[name], top)))
    # study-group readout: mean of the shared axis and the group's own drug axis (baselines already averaged), scored once
    # against S u own set. This is the readout compared across methods.
    for group in ('cisplatin', 'doxorubicin'):
        truth = targets['with_shared'].get(group)
        if truth and 'shared' in axis_mean and group in axis_mean:
            vec = (axis_mean['shared'] + axis_mean[group].reindex(axis_mean['shared'].index)) / 2
            rows.append(dict(axis=f'{group}_group', baseline=BASELINE_MEAN, target_definition='with_shared',
                             **rank_metrics(vec.to_numpy(), vec.index.to_numpy(), truth, top)))
    return rows


# ------------------------------------------------------------------------------------------ MC-ContrastiveVI fit job
def mccvi_fit_record(H: Harness, base: dict, args, *, cfg: str, seed: int, free_bits: float, fit_dir: Path,
                     experiment: str, epochs: int | None = None, patience: int = 15, tag_extra: str = '') -> dict:
    """Fit one canonical MC-ContrastiveVI on ``cfg`` and return the per-fit record (used by exp1 and exp2).

    free_bits == 0 leaves the model configuration untouched (canonical objective). ``epochs`` is the epoch cap
    (default args.epochs: 100, or 2 in smoke; smoke always wins) and ``patience`` the early-stopping patience
    on validation MSE; the record carries epochs_run, best_epoch and converged (early stopping fired before the cap).
    """
    ss, core, fitjob = H.ss, H.core, H.fitjob
    inp, _ = ss.build_inputs(cfg, base)
    epochs = args.epochs if (epochs is None or args.smoke) else epochs
    tag = f'{cfg}_fb{free_bits:g}{tag_extra}_s{seed}'
    kw = {'free_bits': float(free_bits)} if free_bits else None
    res = fitjob.run_fit(inp, seed, fit_dir, tag, save_checkpoint=False, save_latents=False, epochs=epochs, patience=patience,
                         extra={'config': cfg, 'free_bits': float(free_bits), 'threads': N_THREADS, 'experiment': experiment},
                         model_kwargs=kw, log=lambda s: print(f'    {s}', flush=True))
    scen, delta, rf = ss.parse_config(cfg)
    rec = dict(res)
    rec.update(free_bits=float(free_bits), config=cfg, seed=seed, scenario=scen, delta=delta, responder_fraction=rf,
               smoke=bool(args.smoke), **ss.true_shared_fraction(cfg, base))
    rec['max_epochs'] = epochs
    rec['shared_kl_per_dim'] = res['shared_mean_total_kl_nats'] / 8          # latent_dim_shared = 8
    rec['probe_balanced_accuracy'] = {k: v for k, v in res.items() if k.endswith('_balanced_accuracy')}
    rec['active_units'] = {k: v for k, v in res.items() if k.endswith('_active_units_var_gt_0p01')}
    per_dim = pd.read_csv(fit_dir / f'{tag}_latent_usage_per_dimension.csv')
    rec['per_dimension'] = per_dim[['block', 'condition', 'dimension', 'n_cells', 'mean_kl_nats',
                                    'posterior_mean_variance', 'active_unit_var_gt_0p01']].to_dict('records')
    igp = fit_dir / f'{tag}_ig_summary.csv.gz'
    rec['ig_recovery'] = ig_summary_rows(igp, cfg, H.gene_sets(), ss, core) if scen != 'null' and igp.exists() else []
    return rec
