"""Exp 6: what F_sh does on pure noise. NO training, pure simulation.

Collapsed latent blocks are modelled as posterior means that are iid Gaussian noise with a tiny, matched variance:
an 8-dimensional shared block (non-zero for every treated cell) and two 4-dimensional drug blocks that are gated
exactly as the model gates them (a cell's non-matching drug coordinates are zero: doxorubicin cells keep drug block 0,
cisplatin cells drug block 1). F_sh is computed by the source code itself (src.peerreview.core.shared_fractions ->
variance_fraction: sum of per-dimension variances of the gated shared block over that of the gated concatenated drug
block, on treated test cells, pooled and per drug) and the draw is repeated 10,000 times.

Sections of <out>/fsh_noise_null.csv:
  main            8 shared + 2 x 4 drug dims, equal per-dimension variance, compared with the observed null fits
  sweep_total16   the same 16 salient dims split differently: shared 2..14, the rest shared out over the 2 drug blocks
  sweep_shared8   shared fixed at 8, drug dims per block 1..16
Also written: fsh_noise_null_draws.csv (the 10,000 main draws) and fsh_noise_null.md (short summary).
With equal variances the expected value is n_shared / (n_shared + n_drug_per_block): the pooled drug variance of a
gated column is (fraction of cells in that drug) x sigma^2 and the two fractions add to 1.

Usage:
  python experiments/exp6_fsh_noise_null.py [--out analysis/experiments] [--reps 10000] [--resume]
  python experiments/exp6_fsh_noise_null.py --smoke
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from experiments import common  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

OBSERVED_NULL = [0.705, 0.802, 0.869]
SWEEP_TOTAL16 = [(2, 7), (4, 6), (6, 5), (8, 4), (10, 3), (12, 2), (14, 1)]
SWEEP_SHARED8 = [(8, 1), (8, 2), (8, 4), (8, 8), (8, 16)]
SEED = 20261003


def simulate(core, ds, dd, n_dox, n_cis, reps, rng, var_shared, var_drug):
    """reps draws of (pooled, doxorubicin, cisplatin) F_sh from gated iid Gaussian posterior means."""
    n = n_dox + n_cis
    d = np.concatenate([np.ones(n_dox, dtype=np.int64), 2 * np.ones(n_cis, dtype=np.int64)])
    te = np.arange(n)
    out = np.empty((reps, 3))
    for r in range(reps):
        shared = rng.normal(0.0, np.sqrt(var_shared), size=(n, ds))           # all cells treated: gate is 1
        drug = np.zeros((n, 2 * dd))
        drug[:n_dox, :dd] = rng.normal(0.0, np.sqrt(var_drug), size=(n_dox, dd))   # drug block 0 on doxorubicin cells
        drug[n_dox:, dd:] = rng.normal(0.0, np.sqrt(var_drug), size=(n_cis, dd))   # drug block 1 on cisplatin cells
        f = core.shared_fractions({'shared': shared, 'drug': drug}, d, te)
        out[r] = (f['pooled'], f['doxorubicin'], f['cisplatin'])
    return out


def summarize(draws, section, ds, dd, n_dox, n_cis, var_shared, var_drug):
    rows = []
    expected = ds * var_shared / (ds * var_shared + dd * var_drug)             # pooled; per-drug is the same in expectation
    for j, kind in enumerate(('pooled', 'doxorubicin', 'cisplatin')):
        v = draws[:, j]
        rows.append(dict(section=section, f_sh_kind=kind, n_shared_dims=ds, n_drug_dims_per_block=dd,
                         shared_to_drug_dim_ratio=ds / dd, n_dox_cells=n_dox, n_cis_cells=n_cis,
                         var_shared=var_shared, var_drug=var_drug, reps=len(v), mean=float(v.mean()), sd=float(v.std(ddof=1)),
                         min=float(v.min()), q025=float(np.quantile(v, .025)), median=float(np.median(v)),
                         q975=float(np.quantile(v, .975)), max=float(v.max()), analytic_expected=expected,
                         **{f'p_ge_{o}': float((v >= o).mean()) for o in OBSERVED_NULL}))
    return rows


def implied_variance_ratio(f, ds, dd):
    """Shared:drug per-dimension variance ratio that makes the expected F_sh equal f (analytic)."""
    return f * dd / (ds * (1.0 - f))


def write_markdown(path, df, args):
    main = df[df.section == 'main'].set_index('f_sh_kind'); p = main.loc['pooled']
    sweep = df[(df.section == 'sweep_total16') & (df.f_sh_kind == 'pooled')].sort_values('n_shared_dims')
    above = all(o > p['max'] for o in OBSERVED_NULL)
    lines = ['# F_sh on pure noise', '',
             f'Simulation: {int(p.reps)} draws, {int(p.n_dox_cells)} doxorubicin + {int(p.n_cis_cells)} cisplatin treated test cells, '
             f'collapsed blocks = iid Gaussian posterior means (variance {p.var_shared:g} in every block), drug blocks gated as in the model '
             f'(a cell\'s non-matching drug coordinates are zero), F_sh from `src.peerreview.core.shared_fractions`. '
             f'The cell counts are an assumption; the expected value does not depend on them, the spread does.', '',
             '## Observed null fits versus matched-variance noise', '',
             (f'**All three observed null values ({", ".join(map(str, OBSERVED_NULL))}) exceed the matched-variance simulation.** '
              f'With 8 shared + two 4-dim drug blocks at equal variance the simulated pooled F_sh has mean {p["mean"]:.4f} '
              f'(analytic 8/(8+4) = {p.analytic_expected:.4f}), 97.5th percentile {p.q975:.4f} and maximum {p["max"]:.4f} over {int(p.reps):,} draws, '
              f'so none of the draws reaches even the smallest observed value.' if above else
              f'Not all observed null values exceed the simulated maximum ({p["max"]:.4f}); see the P column.'), '',
             'Shared:drug per-dimension variance ratio r that would give each observed value (analytic, E[F] = 8r / (8r + 4) so '
             'r = 4F / (8(1 - F))); r = 1 is the matched-variance case simulated here:', '',
             '| observed null F_sh | P(noise >= observed) | required r |', '|---|---|---|']
    for o in OBSERVED_NULL:
        lines.append(f'| {o} | {p[f"p_ge_{o}"]:.4f} | {implied_variance_ratio(o, 8, 4):.2f} |')
    lines += ['', 'The ratio says how much more variance per dimension the shared block would need than a drug dimension for noise alone to '
              'produce the observed value; it does not say the shared block has it.', '',
              '## Dimension split alone: total latent budget 16, shared dimensions 2 to 14', '',
              'Shared dimensions s and two drug blocks of (16 - s) / 2 dimensions each, equal per-dimension variance, pooled F_sh. '
              'Simulated mean with its analytic expectation s / (s + (16 - s) / 2) alongside:', '',
              '| shared dims | drug dims per block | shared:drug dims | simulated mean | sd | 2.5% | 97.5% | analytic |', '|---|---|---|---|---|---|---|---|']
    for r in sweep.itertuples():
        lines.append(f'| {r.n_shared_dims} | {r.n_drug_dims_per_block} | {r.shared_to_drug_dim_ratio:.2f} | {r.mean:.4f} | {r.sd:.4f} | {r.q025:.4f} | {r.q975:.4f} | {r.analytic_expected:.4f} |')
    equiv = ', '.join(f'{16 * o / (2 - o):.1f}' for o in OBSERVED_NULL)      # F = 2s / (s + 16) at equal variances, 16 dims in total
    lines += ['', f'The model\'s own split (8 shared, 4 + 4 drug) sits at 0.667. At equal variances the observed null values would need '
              f'{equiv} shared dimensions out of 16 (F = 2s / (s + 16)), against the 8 the model has.', '',
              '## Main case, all three F_sh definitions', '',
              '| F_sh | mean | sd | 2.5% | median | 97.5% | analytic E |', '|---|---|---|---|---|---|---|']
    for kind in ('pooled', 'doxorubicin', 'cisplatin'):
        r = main.loc[kind]
        lines.append(f'| {kind} | {r["mean"]:.4f} | {r["sd"]:.4f} | {r.q025:.4f} | {r["median"]:.4f} | {r.q975:.4f} | {r.analytic_expected:.4f} |')
    s8 = df[(df.section == 'sweep_shared8') & (df.f_sh_kind == 'pooled')]
    lines += ['', '## Drug block size with shared fixed at 8 (pooled F_sh, equal variances)', '',
              '| shared dims | drug dims per block | shared:drug | mean | 2.5% | 97.5% | analytic |', '|---|---|---|---|---|---|---|']
    for r in s8.itertuples():
        lines.append(f'| {r.n_shared_dims} | {r.n_drug_dims_per_block} | {r.shared_to_drug_dim_ratio:.2f} | {r.mean:.4f} | {r.q025:.4f} | {r.q975:.4f} | {r.analytic_expected:.4f} |')
    Path(path).write_text('\n'.join(lines) + '\n')


def main(argv=None):
    ap = common.make_parser(__doc__, default_seeds=1, default_out='analysis/experiments')
    ap.add_argument('--reps', type=int, default=10000)
    ap.add_argument('--n-dox', type=int, default=358, help='treated doxorubicin test cells (assumed; real counts: scripts/31 base.npz)')
    ap.add_argument('--n-cis', type=int, default=357, help='treated cisplatin test cells (assumed)')
    ap.add_argument('--sigma2', type=float, default=1e-4, help='per-dimension noise variance (F_sh is scale free when matched)')
    args = common.parse(ap, argv)                                    # --seeds is accepted for CLI uniformity and unused here
    if args.smoke:
        args.reps, args.n_dox, args.n_cis = 100, 20, 20
    from src.peerreview import core
    runner = common.Runner('exp6', args)
    cases = [('main', 8, 4)] + [('sweep_total16', ds, dd) for ds, dd in SWEEP_TOTAL16] + [('sweep_shared8', ds, dd) for ds, dd in SWEEP_SHARED8]
    jobs = [(f'{sec}_ds{ds}_dd{dd}', dict(sec=sec, ds=ds, dd=dd, idx=i)) for i, (sec, ds, dd) in enumerate(cases)]
    jobs = common.select_jobs(jobs, args) if args.max_jobs else jobs

    def fn(job_id, sec, ds, dd, idx):
        rng = np.random.default_rng([SEED, idx])
        draws = simulate(core, ds, dd, args.n_dox, args.n_cis, args.reps, rng, args.sigma2, args.sigma2)
        rec = dict(section=sec, rows=summarize(draws, sec, ds, dd, args.n_dox, args.n_cis, args.sigma2, args.sigma2))
        if sec == 'main':
            rec['draws'] = draws.tolist()
        return rec

    failed = runner.run(jobs, fn, resume=args.resume)
    results = runner.load_results([j for j, _ in jobs])
    rows = [r for res in results for r in res['rows']]
    if rows:
        df = pd.DataFrame(rows); args.out.mkdir(parents=True, exist_ok=True)
        df.to_csv(args.out / 'fsh_noise_null.csv', index=False)
        for res in results:
            if 'draws' in res:
                pd.DataFrame(res['draws'], columns=['pooled', 'doxorubicin', 'cisplatin']).to_csv(args.out / 'fsh_noise_null_draws.csv', index=False)
        if {'main', 'sweep_total16', 'sweep_shared8'} <= set(df.section):
            write_markdown(args.out / 'fsh_noise_null.md', df, args)
        print(f'wrote {args.out / "fsh_noise_null.csv"}', flush=True)
    return runner.stats


if __name__ == '__main__':
    sys.exit(1 if main()['failed'] else 0)
