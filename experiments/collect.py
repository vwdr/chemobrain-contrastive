"""Collect every fit (new experiments and the archived ones) into analysis/experiments/master.csv and write STATS.md.

Nothing is trained or refitted. Sources (anything missing is skipped and listed under Caveats):
  new       runs/experiments/exp{1,2,3,4}/jobs/exp*/*.json            (--exp-root)
  archived  runs/peerreview_20261003/semisynthetic/<cfg>/<cfg>_s<seed>_{done.json,ig_summary.csv.gz}   MC-ContrastiveVI
            analysis/peerreview_20261004/R4_{metrics,ig_recovery}_per_fit.csv + runs/.../r4/*/ *_ig.csv.gz  contrastiveVI, multiGroupVI
            (--root is the repo root holding runs/ and analysis/; the R4 rows for MC-ContrastiveVI are not used because the
            same fits are re-read from their IG summaries, so each fit enters once)

master.csv: one row per fit. A "fit" is a MC-ContrastiveVI fit, a contrastiveVI fit (one per study), a multiGroupVI fit
(both groups averaged) or a scVI / PCA+Wilcoxon job. Required columns: experiment, method, config, scenario, delta,
responder_fraction, free_bits, seed, shared_fraction, active_units, n_dims, kl_per_dim, total_kl,
probe_balanced_accuracy, test_mse, recovery_top50_enrichment, recovery_auroc, recovery_ap. Appended: fit_id, group,
source, recovery_auroc_all_axes, recovery_auroc_shared_axis, epochs_run, best_epoch, max_epochs, converged.

Definitions (choices of this script, not facts from the original analysis):
  * baselines: the two IG baselines of a fit are averaged gene by gene (mean|IG|) and scored once; truth for a study group =
    S (if the scenario injects it) u that group's own set (cisplatin A, doxorubicin B, if injected). Top-50 enrichment and
    AUROC are always reported together. contrastiveVI / multiGroupVI per-gene files are matched with their own file-name
    patterns (<cfg>_<study>_s<seed>_<baseline>_ig.csv.gz vs <cfg>_s<seed>_<group>_<baseline>_ig.csv.gz).
  * regression gate: group-level medians of the archived fits (delta 1 and 2, rf 1.0) must reproduce the expected values
    (MC-ContrastiveVI 0.461 / 0.75, contrastiveVI 0.479 / 1.50, multiGroupVI 0.528 / 3.75) or STATS.md is not written.
  * recovery = the group readout shared by all methods: the cisplatin and doxorubicin axes/groups, truth = genes
    injected in that configuration (S, plus A for cisplatin / B for doxorubicin, only those present in the scenario).
    Averaged over groups/axes and over the zero / control_median baselines. null configurations have no recovery.
    For MC-ContrastiveVI the group ranking is the mean of the shared axis and the group's own drug axis of the
    ig_summary files (baselines averaged first); the single-axis AUROCs are reported as recovery_auroc_shared_axis and
    recovery_auroc_all_axes, and the robustness of the readout choice in STATS (six readouts).
  * recovery_top50_enrichment = injected genes in the top-50 divided by their random expectation 50 * n_injected / n_genes.
  * active_units / total_kl / n_dims: MC-ContrastiveVI = shared + doxorubicin + cisplatin blocks (16 dims, variance of
    posterior mean > 0.01; summed per-block mean KL); contrastiveVI = salient block (8); multiGroupVI = the two private
    blocks (20); scVI / PCA: not defined (NaN). kl_per_dim = total_kl / n_dims.
  * probe_balanced_accuracy: MC-ContrastiveVI = 3-class condition probe on the ungated drug block
    (drug_ungated_condition_balanced_accuracy); contrastiveVI / multiGroupVI = within-study treated-vs-control probe on
    their private block. These are different tasks and are not strictly comparable across methods.
  * test_mse: MC-ContrastiveVI test_mse; contrastiveVI / multiGroupVI canonical_log_mse (same canonical scale).

Usage:
  python experiments/collect.py [--out analysis/experiments] [--root .] [--exp-root runs/experiments]
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from experiments import common  # noqa: E402

import argparse  # noqa: E402
import json  # noqa: E402
import re  # noqa: E402
from types import SimpleNamespace  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

R = common.R
NAN = float('nan')
OWN = {'cisplatin': 'A', 'doxorubicin': 'B'}
GROUPS = ('cisplatin', 'doxorubicin')
STUDY_GROUP = {'GSE216146': 'cisplatin', 'GSE271055': 'doxorubicin'}
BLOCKS = ('shared', 'dox', 'cis')
N_GENES = 1500
OBSERVED_NULL = [0.705, 0.802, 0.869]
COLUMNS = ['experiment', 'method', 'config', 'scenario', 'delta', 'responder_fraction', 'free_bits', 'seed', 'shared_fraction',
           'active_units', 'n_dims', 'kl_per_dim', 'total_kl', 'probe_balanced_accuracy', 'test_mse',
           'recovery_top50_enrichment', 'recovery_auroc', 'recovery_ap',
           'fit_id', 'group', 'source', 'recovery_auroc_all_axes', 'recovery_auroc_shared_axis',
           'epochs_run', 'best_epoch', 'max_epochs', 'converged']
GROUP_ROWS = []      # group-level (fit x study group) scores of the archived fits, for the regression gate
EXPECTED = {         # medians on the archived fits (delta 1 and 2, responder fraction 1.0, group level); supplied by the user
    'MC-ContrastiveVI': (0.461, 0.75), 'contrastiveVI': (0.479, 1.50), 'multiGroupVI': (0.528, 3.75)}
EXPECTED_P = {   # paired one-sided Wilcoxon vs MC-ContrastiveVI, group level; supplied by the user as approximate
    ('contrastiveVI', 'auroc'): 1.8e-4, ('contrastiveVI', 'enrichment'): 5.7e-5,
    ('multiGroupVI', 'auroc'): 1.0e-10, ('multiGroupVI', 'enrichment'): 1.2e-7}
MC_VARIANTS = []     # MC-ContrastiveVI group scores under six readouts (archived fits), for the robustness table

# Mirrors of scripts/31 (importing it would pull in torch); keep in sync with SETS_BY_SCENARIO / parse_config there.
SETS_BY_SCENARIO = {'null': [], 'shared': ['S'], 'specific': ['A', 'B'], 'mixed': ['S', 'A', 'B']}


def parse_config(name):
    if name == 'null':
        return 'null', 0.0, None
    scen, d, r = name.split('_')
    return scen, float(d[1:]), float(r[2:])


SS = SimpleNamespace(SETS_BY_SCENARIO=SETS_BY_SCENARIO, parse_config=parse_config)
CORE = SimpleNamespace(AXIS_NAME={'shared': 'shared', 'drug0': 'doxorubicin', 'drug1': 'cisplatin'})


# ------------------------------------------------------------------------------------------ small helpers
def mean_or_nan(v):
    v = [x for x in v if x is not None and not (isinstance(x, float) and np.isnan(x))]
    return float(np.mean(v)) if v else NAN


def fmt_p(p):
    return 'nan' if p is None or (isinstance(p, float) and np.isnan(p)) else (f'{p:.2e}' if p < 1e-3 else f'{p:.4f}')


def md_table(df, floatfmt='{:.4g}'):
    if df is None or len(df) == 0:
        return '_(none)_\n'
    cols = list(df.columns)
    out = ['| ' + ' | '.join(cols) + ' |', '|' + '---|' * len(cols)]
    for _, r in df.iterrows():
        out.append('| ' + ' | '.join(floatfmt.format(v) if isinstance(v, (float, np.floating)) and not pd.isna(v) else
                                     ('' if pd.isna(v) else str(v)) for v in r.tolist()) + ' |')
    return '\n'.join(out) + '\n'


SMOKE = False   # --smoke: read the smoke outputs of run_all.py --smoke and skip the archived sources


def jobs_of(root: Path, exp: str):
    d = root / exp / ('smoke' if SMOKE else '') / 'jobs' / exp
    return [json.loads(p.read_text()) for p in sorted(d.glob('*.json'))] if d.exists() else []


def target_sets(cfg, group, sizes):
    scen = parse_config(cfg)[0]
    present = SETS_BY_SCENARIO[scen]
    return [s for s in ('S', OWN[group]) if s in present]


def recovery_from_hits(row, cfg, group, sizes):
    """Top-50 hits in the injected sets (from per-set hit columns) as enrichment over random expectation."""
    sets = target_sets(cfg, group, sizes)
    if not sets:
        return NAN
    hits = sum(row.get(f'top50_hits_{s}', 0) or 0 for s in sets)
    n_pos = sum(sizes[s] for s in sets)
    return hits / (50.0 * n_pos / N_GENES) if n_pos else NAN


# ------------------------------------------------------------------------------------------ MC-ContrastiveVI fits
def aggregate_ig_rows(rows):
    """(group-readout enrichment, AUROC, AP, all-axes AUROC, shared-axis AUROC, per-group list) from ig_summary_rows.

    Only baseline == 'mean_of_baselines' rows and truth = S u own set ('with_shared'). The group readout is the mean of
    the shared axis and the group's own drug axis (rows named <group>_group); all-axes averages the three single axes.
    """
    rows = [r for r in rows if r.get('target_definition') == 'with_shared' and r.get('baseline') == common.BASELINE_MEAN]
    g = [r for r in rows if r['axis'].endswith('_group')]
    axes = [r for r in rows if not r['axis'].endswith('_group')]
    enr = lambda r: r['top_hits'] / r['random_expectation_hits'] if r['random_expectation_hits'] > 0 else NAN
    return (mean_or_nan([enr(r) for r in g]), mean_or_nan([r['auroc'] for r in g]), mean_or_nan([r['average_precision'] for r in g]),
            mean_or_nan([r['auroc'] for r in axes]), mean_or_nan([r['auroc'] for r in axes if r['axis'] == 'shared']),
            [(r['axis'][:-len('_group')], enr(r), r['auroc'], r['average_precision']) for r in g])


def mc_readout_variants(igp, cfg, gs):
    """MC-ContrastiveVI study-group scores under six readouts of the same IG summary (baselines averaged first)."""
    scen, delta, rf = parse_config(cfg)
    df = pd.read_csv(igp)
    if 'target' in df.columns:
        df = df[df.target == 'norm_squared']
    avg = {a: df[df.axis == a].groupby('gene').mean_abs_ig.mean() for a in ('shared', 'drug0', 'drug1')}
    sh = avg['shared']; out = []
    truths = common.injected_targets(cfg, gs, SS)['with_shared']
    for axis, group in (('drug0', 'doxorubicin'), ('drug1', 'cisplatin')):
        truth = truths.get(group)
        if not truth:
            continue
        dr = avg[axis].reindex(sh.index)
        readouts = {'shared only': sh, 'drug only': dr, 'mean': (sh + dr) / 2, 'max': pd.concat([sh, dr], axis=1).max(axis=1),
                    'rank-average': (sh.rank() + dr.rank()) / 2, 'norm-concat': np.sqrt(sh ** 2 + dr ** 2)}
        for name, v in readouts.items():
            m = common.rank_metrics(v.to_numpy(), v.index.to_numpy(), truth)
            out.append(dict(variant=name, config=cfg, delta=delta, responder_fraction=rf, group=group, auroc=m['auroc'],
                            enrichment=m['top_hits'] / m['random_expectation_hits'] if m['random_expectation_hits'] > 0 else NAN))
    return out


def mccvi_row(res, experiment, source, ig_rows):
    cfg = res['config']; scen, delta, rf = parse_config(cfg)
    tot = sum(res.get(f'{b}_mean_total_kl_nats', NAN) for b in BLOCKS)
    n_dims = 16
    enr, auc, ap, auc_all, auc_sh, per_group = aggregate_ig_rows(ig_rows)
    seed = int(res['seed'])
    fb = float(res.get('free_bits', 0.0))
    if source == 'archive':
        for group, e, a, p_ in per_group:
            GROUP_ROWS.append(dict(method='MC-ContrastiveVI', config=cfg, delta=delta, responder_fraction=rf, seed=seed, group=group,
                                   enrichment=e, auroc=a, ap=p_))
    run = res.get('epochs_run')
    return dict(experiment=experiment, method='MC-ContrastiveVI', config=cfg, scenario=scen, delta=delta, responder_fraction=rf,
                free_bits=fb, seed=seed, shared_fraction=res.get('shared_fraction', NAN),
                active_units=sum(res.get(f'{b}_active_units_var_gt_0p01', 0) for b in BLOCKS), n_dims=n_dims,
                kl_per_dim=tot / n_dims, total_kl=tot,
                probe_balanced_accuracy=res.get('drug_ungated_condition_balanced_accuracy', NAN), test_mse=res.get('test_mse', NAN),
                recovery_top50_enrichment=enr, recovery_auroc=auc, recovery_ap=ap, fit_id=f'{experiment}:{cfg}:fb{fb:g}:s{seed}',
                group='both', source=source, recovery_auroc_all_axes=auc_all, recovery_auroc_shared_axis=auc_sh,
                epochs_run=run if run is not None else NAN, best_epoch=res.get('best_epoch', NAN),
                max_epochs=res.get('epochs_requested', res.get('max_epochs', 100)), converged=NAN if not_converged(res) is None else (not not_converged(res)))


def not_converged(res):
    """True when the fit did not plateau. New fits carry `converged`; archived ones ran the 100-epoch cap with patience 15,
    so reaching the cap means early stopping never fired."""
    if 'converged' in res:
        return not res['converged']
    run, cap = res.get('epochs_run'), res.get('epochs_requested', 100)
    return None if run is None else bool(run >= cap)


def convergence_flags(res):
    """(did not converge, IG quadrature not converged) for one MC-ContrastiveVI fit."""
    nc = not_converged(res)
    ig_bad = bool((res.get('ig_max_median_scaled_error') or 0) >= 0.01)
    return bool(nc), ig_bad


def load_archive_mccvi(root, gs, notes):
    base = root / 'runs' / 'peerreview_20261003' / 'semisynthetic'
    rows, flags = [], []
    if not base.exists():
        notes.append(f'archived MC-ContrastiveVI fits not found ({base}); the 72-fit injected set is missing.')
        return rows, flags
    for done in sorted(base.glob('*/*_s[0-9]_done.json')):
        res = json.loads(done.read_text())
        cfg = res.get('config') or done.parent.name
        res['config'] = cfg
        igp = done.with_name(done.name.replace('_done.json', '_ig_summary.csv.gz'))
        ig = common.ig_summary_rows(igp, cfg, gs, SS, CORE) if (gs is not None and igp.exists() and parse_config(cfg)[0] != 'null') else []
        rows.append(mccvi_row(res, 'archive_semisynthetic', 'archive', ig)); flags.append(('archive_semisynthetic', *convergence_flags(res)))
        if gs is not None and igp.exists() and parse_config(cfg)[0] != 'null':
            for v in mc_readout_variants(igp, cfg, gs):
                MC_VARIANTS.append(dict(v, seed=int(res['seed'])))
    return rows, flags


def load_new_mccvi(exp_root, notes):
    rows, flags = [], []
    for exp in ('exp1', 'exp2', 'exp7'):
        recs = jobs_of(exp_root, exp)
        if not recs:
            notes.append(f'no {exp} results under {exp_root / exp}.')
        for r in recs:
            rows.append(mccvi_row(r, exp, 'new', r.get('ig_recovery', []))); flags.append((exp, *convergence_flags(r)))
    return rows, flags


# ------------------------------------------------------------------------------------------ external methods
IG_CVI = re.compile(r'^(?P<cfg>.+)_(?P<study>GSE\d+)_s(?P<seed>\d+)_(?P<base>zero|control_median)_ig\.csv\.gz$')
IG_MGVI = re.compile(r'^(?P<cfg>.+)_s(?P<seed>\d+)_(?P<group>cisplatin|doxorubicin)_(?P<base>zero|control_median)_ig\.csv\.gz$')


def vector_scores(root_dir, gs):
    """{(method, cfg, seed, group): (AUROC, AP, top-50 enrichment)} from the per-gene *_ig.csv.gz files of both methods.

    The two baselines of a fit are averaged gene by gene (mean|IG|) and the average is scored once, against S u the
    group's own injected set. File names differ: contrastiveVI <cfg>_<study>_s<seed>_<baseline>_ig.csv.gz,
    multiGroupVI <cfg>_s<seed>_<group>_<baseline>_ig.csv.gz; each method is matched with its own pattern.
    """
    acc = {}
    for method, sub, pat in (('contrastiveVI', 'contrastivevi', IG_CVI), ('multiGroupVI', 'multigroupvi', IG_MGVI)):
        d = root_dir / sub / ('smoke' if SMOKE else '')
        for p in sorted(d.glob('*_ig.csv.gz')) if d.exists() else []:
            m = pat.match(p.name)
            if not m:
                continue
            cfg, seed = m.group('cfg'), int(m.group('seed'))
            group = STUDY_GROUP[m.group('study')] if method == 'contrastiveVI' else m.group('group')
            if parse_config(cfg)[0] == 'null':
                continue
            acc.setdefault((method, cfg, seed, group), []).append(pd.read_csv(p).set_index('gene').mean_abs_ig)
    out = {}
    for (method, cfg, seed, group), vs in acc.items():
        truth = common.injected_targets(cfg, gs, SS)['with_shared'].get(group)
        if not truth:
            continue
        v = pd.concat(vs, axis=1).mean(axis=1)
        rm = common.rank_metrics(v.to_numpy(), v.index.to_numpy(), truth)
        enr = rm['top_hits'] / rm['random_expectation_hits'] if rm['random_expectation_hits'] > 0 else NAN
        out[(method, cfg, seed, group)] = (rm['auroc'], rm['average_precision'], enr)
    return out


def external_rows(met_long, vec, experiment, source):
    """Fit-level rows for contrastiveVI (one per study group) and multiGroupVI (the two groups averaged)."""
    rows = []
    if not vec:
        return rows
    met = met_long.set_index(['method', 'config', 'seed', 'group']) if not met_long.empty else None
    per_group = []
    for (method, cfg, seed, group), (auc, ap, enr) in vec.items():
        m = met.loc[(method, cfg, seed, group)] if met is not None and (method, cfg, seed, group) in met.index else None
        per_group.append(dict(method=method, config=cfg, seed=int(seed), group=group, enrich=enr, auroc=auc, ap=ap,
                              active_units=m.active_units if m is not None else NAN, n_dims=m.n_dims if m is not None else NAN,
                              total_kl=m.total_kl_nats if m is not None else NAN,
                              probe=m.probe_within_study if m is not None else NAN, mse=m.canonical_log_mse if m is not None else NAN))
        if source == 'archive':
            scen, delta, rf = parse_config(cfg)
            GROUP_ROWS.append(dict(method=method, config=cfg, delta=delta, responder_fraction=rf, seed=int(seed), group=group,
                                   enrichment=enr, auroc=auc, ap=ap))
    pg = pd.DataFrame(per_group)
    for (method, cfg, seed), sub in pg.groupby(['method', 'config', 'seed']):
        scen, delta, rf = parse_config(cfg)
        units = [sub] if method == 'multiGroupVI' else [sub[sub.group == gname] for gname in sub.group]   # cvi: one fit per study
        for u in units:
            if method == 'multiGroupVI':
                n_dims = u.n_dims.sum(skipna=False); au = u.active_units.sum(skipna=False); tk = u.total_kl.sum(skipna=False); grp = 'both'
            else:
                n_dims, au, tk, grp = u.n_dims.iloc[0], u.active_units.iloc[0], u.total_kl.iloc[0], u.group.iloc[0]
            rows.append(dict(experiment=experiment, method=method, config=cfg, scenario=scen, delta=delta, responder_fraction=rf,
                             free_bits=NAN, seed=int(seed), shared_fraction=NAN, active_units=au, n_dims=n_dims,
                             kl_per_dim=tk / n_dims if n_dims else NAN, total_kl=tk, probe_balanced_accuracy=u.probe.mean(),
                             test_mse=u.mse.mean(), recovery_top50_enrichment=u.enrich.mean(), recovery_auroc=u.auroc.mean(),
                             recovery_ap=u.ap.mean(), fit_id=f'{experiment}:{method}:{cfg}:{grp}:s{seed}', group=grp, source=source,
                             recovery_auroc_all_axes=NAN, recovery_auroc_shared_axis=NAN,
                             epochs_run=NAN, best_epoch=NAN, max_epochs=NAN, converged=NAN))
    return rows


def load_archive_external(root, gs, notes):
    a = root / 'analysis' / 'peerreview_20261004'
    if gs is None or not (a / 'R4_metrics_per_fit.csv').exists():
        notes.append(f'archived R4 metrics or gene sets not found under {a}; archived contrastiveVI / multiGroupVI fits are missing.')
        return []
    rd = lambda p: pd.read_csv(p, keep_default_na=False, na_values=[''])
    met = rd(a / 'R4_metrics_per_fit.csv'); met = met[met.model.isin(['contrastiveVI', 'multiGroupVI'])].rename(columns={'model': 'method'})
    vec = vector_scores(root / 'runs' / 'peerreview_20261004' / 'r4', gs)
    if not vec:
        notes.append('archived per-gene IG files (runs/peerreview_20261004/r4/*/*_ig.csv.gz) not found: no archived contrastiveVI / multiGroupVI recovery.')
    return external_rows(met, vec, 'archive_r4', 'archive')


def load_exp3(exp_root, gs, notes):
    recs = jobs_of(exp_root, 'exp3')
    if not recs:
        notes.append(f'no exp3 results under {exp_root / "exp3"}.')
        return []
    met_rows = []
    for r in recs:
        m = r['metrics']; method = {'contrastivevi': 'contrastiveVI', 'multigroupvi': 'multiGroupVI'}[r['method']]
        if method == 'contrastiveVI':
            met_rows.append(dict(method=method, config=r['config'], seed=r['seed'], group=STUDY_GROUP[r['study']],
                                 active_units=m['salient_active_units'], n_dims=8, total_kl_nats=m['salient_total_kl_nats'],
                                 probe_within_study=m['probe_within_study_salient_ungated'], canonical_log_mse=m['canonical_log_mse']))
        else:
            for label in GROUPS:
                met_rows.append(dict(method=method, config=r['config'], seed=r['seed'], group=label,
                                     active_units=m[f'{label}_private_active_units'], n_dims=10,
                                     total_kl_nats=m[f'{label}_private_total_kl_nats'],
                                     probe_within_study=m[f'probe_within_{label}_study_private_ungated'], canonical_log_mse=m['canonical_log_mse']))
    vec = vector_scores(exp_root / 'exp3' / ('smoke' if SMOKE else ''), gs) if gs is not None else {}
    if not vec:
        notes.append('exp3 per-gene IG files not found next to the JSONs: no exp3 recovery.')
    return external_rows(pd.DataFrame(met_rows), vec, 'exp3', 'new')


def load_exp4(exp_root, notes):
    recs = jobs_of(exp_root, 'exp4')
    if not recs:
        notes.append(f'no exp4 results under {exp_root / "exp4"}.')
        return []
    primary = {'scvi': ('scVI', 'abs_lfc_mean'), 'pca_wilcoxon': ('PCA32+Wilcoxon', 'wilcoxon_pca32_recon')}
    sizes = {'S': 40, 'A': 40, 'B': 40}
    rows = []
    for r in recs:
        name, ranking = primary[r['method']]
        sel = [x for x in r['recovery'] if x['ranking'] == ranking]
        if not sel:
            continue
        cfg = r['config']
        enr = mean_or_nan([recovery_from_hits(x, cfg, x['group'], sizes) for x in sel])
        scen, delta, rf = parse_config(cfg)
        rows.append(dict(experiment='exp4', method=name, config=cfg, scenario=scen, delta=delta, responder_fraction=rf, free_bits=NAN,
                         seed=int(r['seed']), shared_fraction=NAN, active_units=NAN, n_dims=NAN, kl_per_dim=NAN, total_kl=NAN,
                         probe_balanced_accuracy=NAN, test_mse=NAN, recovery_top50_enrichment=enr,
                         recovery_auroc=mean_or_nan([x['auroc'] for x in sel]), recovery_ap=mean_or_nan([x['average_precision'] for x in sel]),
                         fit_id=f'exp4:{name}:{cfg}:s{r["seed"]}', group='both', source='new', recovery_auroc_all_axes=NAN,
                         recovery_auroc_shared_axis=NAN, epochs_run=NAN, best_epoch=NAN, max_epochs=r.get('epochs', NAN), converged=NAN))
    return rows


# ------------------------------------------------------------------------------------------ statistics
def spearman(x, y):
    d = pd.DataFrame({'x': x, 'y': y}).dropna()
    if len(d) < 3 or d.x.nunique() < 2 or d.y.nunique() < 2:
        return NAN, NAN, len(d)
    rho, p = stats.spearmanr(d.x, d.y)
    return float(rho), float(p), len(d)


def holm(ps):
    order = np.argsort(ps); m = len(ps); adj = np.empty(m); run = 0.0
    for rank, i in enumerate(order):
        run = max(run, (m - rank) * ps[i]); adj[i] = min(1.0, run)
    return adj


def dose_response(master):
    e1 = master[master.experiment == 'exp1']
    inj = e1[e1.scenario != 'null']
    rows = []
    specs = [('recovery AUROC (group axes)', inj, 'recovery_auroc'), ('recovery top-50 enrichment', inj, 'recovery_top50_enrichment'),
             ('recovery AP', inj, 'recovery_ap'), ('recovery AUROC (all three axes)', inj, 'recovery_auroc_all_axes')]
    for label, d, col in specs:
        for pred in ('free_bits', 'active_units', 'kl_per_dim'):
            rho, p, n = spearman(d[pred], d[col]); rows.append(dict(outcome=label, predictor=pred, scope='exp1 injected configs pooled', rho=rho, P=p, n=n))
    for cfg, d in inj.groupby('config'):
        for pred in ('free_bits', 'active_units', 'kl_per_dim'):
            rho, p, n = spearman(d[pred], d.recovery_auroc); rows.append(dict(outcome='recovery AUROC (group axes)', predictor=pred, scope=f'exp1 {cfg}', rho=rho, P=p, n=n))
    sf = []
    for scope, d in (('exp1 all configs', e1), ('exp1 injected configs', inj)):
        rho, p, n = spearman(d.free_bits, d.shared_fraction); sf.append(dict(outcome='shared_fraction', predictor='free_bits', scope=scope, rho=rho, P=p, n=n))
    for cfg, d in e1.groupby('config'):
        rho, p, n = spearman(d.free_bits, d.shared_fraction); sf.append(dict(outcome='shared_fraction', predictor='free_bits', scope=f'exp1 {cfg}', rho=rho, P=p, n=n))
    return pd.DataFrame(rows), pd.DataFrame(sf)


def null_calibration(master):
    null = master[master.experiment == 'exp2'].shared_fraction.dropna().to_numpy()
    inj = master[(master.experiment == 'archive_semisynthetic') & (master.scenario != 'null')].shared_fraction.dropna().to_numpy()
    out = dict(n_null=len(null), n_inj=len(inj))
    if len(null):
        out.update(null_mean=float(null.mean()), null_sd=float(null.std(ddof=1)) if len(null) > 1 else NAN,
                   null_p05=float(np.quantile(null, .05)), null_p95=float(np.quantile(null, .95)))
    if len(null) and len(inj):
        u, p = stats.mannwhitneyu(inj, null, alternative='two-sided')
        out.update(inj_mean=float(inj.mean()), inj_median=float(np.median(inj)), U=float(u), P=float(p), auc=float(u / (len(inj) * len(null))),
                   frac_inj_above_null_p95=float((inj > out['null_p95']).mean()))
    return out


def cross_method(master, deltas, metric):
    d = master[(master.scenario != 'null') & master.delta.isin(deltas) & (master.responder_fraction == 1.0)].copy()
    d = d[d[metric].notna()]
    unit = d.groupby(['method', 'config', 'seed'], as_index=False)[metric].mean()      # groups of one fit already averaged; cvi studies averaged here
    methods = sorted(unit.method.unique())
    res = dict(n={m: int((unit.method == m).sum()) for m in methods})
    groups = [unit[unit.method == m][metric].to_numpy() for m in methods]
    if len(groups) >= 2 and all(len(g) >= 2 for g in groups):
        h, p = stats.kruskal(*groups); res.update(H=float(h), P=float(p), df=len(groups) - 1)
    ref = 'MC-ContrastiveVI'; rows = []
    if ref in methods:
        a = unit[unit.method == ref].set_index(['config', 'seed'])[metric]
        for m in methods:
            if m == ref:
                continue
            b = unit[unit.method == m].set_index(['config', 'seed'])[metric]
            j = pd.concat([a.rename('ref'), b.rename('oth')], axis=1, join='inner').dropna()
            if len(j) >= 6 and (j.oth - j.ref).abs().sum() > 0:
                w, p = stats.wilcoxon(j.oth, j.ref); rows.append(dict(method=m, n_pairs=len(j), median_diff_vs_MC=float((j.oth - j.ref).median()), W=float(w), P=float(p)))
            else:
                rows.append(dict(method=m, n_pairs=len(j), median_diff_vs_MC=float((j.oth - j.ref).median()) if len(j) else NAN, W=NAN, P=NAN))
        ps = [r['P'] for r in rows if not np.isnan(r['P'])]
        if ps:
            adj = iter(holm(np.array(ps)))
            for r in rows:
                r['P_holm'] = float(next(adj)) if not np.isnan(r['P']) else NAN
    res['wilcoxon'] = pd.DataFrame(rows)
    res['means'] = unit.groupby('method')[metric].agg(['mean', 'std', 'count']).reset_index()
    return res


def exp6_summary(out_dir, exp2_mean):
    p = out_dir / 'fsh_noise_null.csv'
    if not p.exists():
        return None
    df = pd.read_csv(p); main = df[(df.section == 'main') & (df.f_sh_kind == 'pooled')]
    if main.empty:
        return None
    r = main.iloc[0]
    out = dict(reps=int(r.reps), mean=float(r['mean']), sd=float(r.sd), q025=float(r.q025), q975=float(r.q975), median=float(r['median']),
               analytic=float(r.analytic_expected), obs_mean=float(np.mean(OBSERVED_NULL)), n_dox=int(r.n_dox_cells), n_cis=int(r.n_cis_cells),
               p_ge={o: float(r[f'p_ge_{o}']) for o in OBSERVED_NULL}, exp2_mean=exp2_mean)
    out['diff_obs'] = out['obs_mean'] - out['mean']
    out['z_obs'] = out['diff_obs'] / out['sd'] if out['sd'] > 0 else NAN
    sweep = df[(df.section == 'sweep_total16') & (df.f_sh_kind == 'pooled')][['n_shared_dims', 'n_drug_dims_per_block', 'mean']]
    out['sweep'] = sweep
    return out


# ------------------------------------------------------------------------------------------ STATS.md
def paired_group_tests():
    """One-sided paired Wilcoxon (method > MC-ContrastiveVI) on matching (config, seed, group), archived fits, delta 1 and 2, rf 1.0."""
    g = pd.DataFrame(GROUP_ROWS)
    rows = []
    if g.empty:
        return pd.DataFrame(rows)
    g = g[g.delta.isin([1.0, 2.0]) & (g.responder_fraction == 1.0)]
    mc = g[g.method == 'MC-ContrastiveVI'].set_index(['config', 'seed', 'group'])
    for method in ('contrastiveVI', 'multiGroupVI'):
        o = g[g.method == method].set_index(['config', 'seed', 'group'])
        for metric, col in (('auroc', 'auroc'), ('enrichment', 'enrichment')):
            j = pd.concat([o[col].rename('o'), mc[col].rename('mc')], axis=1, join='inner').dropna()
            if len(j) < 6:
                continue
            w, p = stats.wilcoxon(j.o, j.mc, alternative='greater')
            rows.append(dict(method=method, metric=metric, n_groups=len(j), median_method=float(j.o.median()), median_MC=float(j.mc.median()),
                             median_diff=float((j.o - j.mc).median()), n_ties=int((j.o == j.mc).sum()), W=float(w), P_one_sided=float(p),
                             P_expected=EXPECTED_P[(method, metric)]))
    return pd.DataFrame(rows)


def mc_variant_table():
    v = pd.DataFrame(MC_VARIANTS)
    if v.empty:
        return pd.DataFrame()
    v = v[v.delta.isin([1.0, 2.0]) & (v.responder_fraction == 1.0)]
    t = v.groupby('variant').agg(n_groups=('auroc', 'count'), auroc_median=('auroc', 'median'), enrichment_median=('enrichment', 'median'),
                                 auroc_min=('auroc', 'min'), auroc_max=('auroc', 'max')).reset_index()
    order = ['shared only', 'drug only', 'mean', 'max', 'rank-average', 'norm-concat']
    return t.set_index('variant').loc[[o for o in order if o in set(t.variant)]].reset_index()


def build_stats(master, flags, notes, out_dir, exp_root, args):
    L = ['# Statistics', '', f'Generated by `experiments/collect.py` from `master.csv` ({len(master)} fits). All tests are two-sided and '
         'uncorrected unless stated. Numbers come straight from the files; no interpretation is built into the tables.', '']
    counts = master.groupby(['experiment', 'method']).size().rename('n_fits').reset_index()
    L += ['## 0. Fits included', '', md_table(counts)]

    dr, sf = dose_response(master)
    L += ['## 1. Within-model dose-response (exp1 only, architecture fixed)', '',
          'Spearman correlation across exp1 fits (free_bits in {0, .05, .1, .25, .5, 1} x configs x seeds). Recovery needs injected '
          'truth, so the `null` configuration drops out of the recovery rows. active_units = shared + dox + cis units with posterior-mean '
          'variance > 0.01 (of 16); kl_per_dim = summed block KL / 16.', '', md_table(dr), '',
          '### shared_fraction against free_bits', '', md_table(sf)]
    pooled = dr[(dr.scope == 'exp1 injected configs pooled') & (dr.outcome == 'recovery AUROC (group axes)') & (dr.predictor == 'free_bits')]
    s_all = sf[sf.scope == 'exp1 all configs']
    if len(pooled) and len(s_all):
        rec_move = bool(pooled.P.iloc[0] < 0.05); sh_flat = bool(s_all.P.iloc[0] >= 0.05)
        L += ['', f'Mechanical check of the stated prediction (recovery moves with free_bits, shared_fraction does not; alpha = 0.05, '
              f'pooled rows): recovery rho = {pooled.rho.iloc[0]:.3f} (P = {fmt_p(pooled.P.iloc[0])}, n = {int(pooled.n.iloc[0])}); '
              f'shared_fraction rho = {s_all.rho.iloc[0]:.3f} (P = {fmt_p(s_all.P.iloc[0])}, n = {int(s_all.n.iloc[0])}). '
              f'Recovery moves: {"yes" if rec_move else "no"}; shared_fraction flat: {"yes" if sh_flat else "no"}. '
              'A non-significant P with this n is absence of evidence, not evidence of flatness (see caveats).', '']
    elif not len(pooled):
        L += ['', '_exp1 results not found or without recovery; section 1 is empty._', '']

    e7 = master[master.experiment == 'exp7']
    L += ['### 1b. Epoch ablation (exp7: free_bits = 0, shared_d1_rf1.0, patience 30)', '']
    if len(e7):
        g = e7.groupby('max_epochs')
        t = g.agg(n=('seed', 'count'), n_converged=('converged', 'sum'), epochs_run=('epochs_run', 'mean'), best_epoch=('best_epoch', 'mean'),
                  active_units=('active_units', 'mean'), kl_per_dim=('kl_per_dim', 'mean'), shared_fraction=('shared_fraction', 'mean'),
                  recovery_auroc=('recovery_auroc', 'mean'), recovery_top50_enrichment=('recovery_top50_enrichment', 'mean')).reset_index()
        L += [md_table(t), 'Means over seeds. `active_units` is shared + dox + cis (of 16); `kl_per_dim` is total KL / 16 over all three blocks (the shared block '
              'alone per dimension is in exp7/epoch_ablation_summary.csv). n_converged counts fits where early stopping fired before the cap.', '']
    else:
        L += ['_exp7 results not found._', '']

    nc = null_calibration(master)
    L += ['## 2. Null calibration of shared_fraction', '']
    if 'null_mean' in nc:
        L += [f'exp2 (n = {nc["n_null"]} uninjected fits, default settings): mean {nc["null_mean"]:.4f}, sd {nc["null_sd"]:.4f}, '
              f'5th-95th percentile [{nc["null_p05"]:.4f}, {nc["null_p95"]:.4f}].', '']
    else:
        L += ['_exp2 results not found._', '']
    if 'U' in nc:
        L += [f'Mann-Whitney, {nc["n_inj"]} injected archived fits (semi-synthetic grid) vs the {nc["n_null"]} null fits: U = {nc["U"]:.1f}, '
              f'P = {fmt_p(nc["P"])}, probability of superiority (U / n1 n2) = {nc["auc"]:.3f}; injected mean {nc["inj_mean"]:.4f}, median '
              f'{nc["inj_median"]:.4f}; {100 * nc["frac_inj_above_null_p95"]:.1f}% of injected fits lie above the null 95th percentile.', '']
        if nc['n_inj'] != 72:
            L += [f'_Note: {nc["n_inj"]} injected fits were found, not the expected 72 (24 configurations x 3 seeds)._', '']
    else:
        L += ['_The injected-vs-null test needs both the archived injected fits and exp2; one is missing._', '']

    e6 = exp6_summary(out_dir, nc.get('null_mean', NAN))
    L += ['## 3. Simulated F_sh noise null (exp6)', '']
    if e6:
        L += [f'{e6["reps"]} draws, 8 shared + 2 x 4 gated drug dims, equal variance, {e6["n_dox"]} + {e6["n_cis"]} treated cells (assumed counts): '
              f'mean {e6["mean"]:.4f} (analytic 8/(8+4) = {e6["analytic"]:.4f}), sd {e6["sd"]:.4f}, 2.5-97.5% [{e6["q025"]:.4f}, {e6["q975"]:.4f}].', '',
              f'Observed null fits ({", ".join(map(str, OBSERVED_NULL))}; mean {e6["obs_mean"]:.4f}) are {e6["diff_obs"]:+.4f} from the simulated mean '
              f'({e6["z_obs"]:+.1f} simulated sd). P(noise >= observed): ' + ', '.join(f'{o}: {fmt_p(p)}' for o, p in e6['p_ge'].items()) + '.']
        if not np.isnan(e6['exp2_mean']):
            L += ['', f'The exp2 null mean is {e6["exp2_mean"]:.4f}, {e6["exp2_mean"] - e6["mean"]:+.4f} from the simulated mean.']
        L += ['', 'F_sh of pure noise as the shared:drug dimension split changes (16 dims in total):', '', md_table(e6['sweep']), '']
    else:
        L += ['_exp6 output (`fsh_noise_null.csv`) not found._', '']

    pt = paired_group_tests(); mv = mc_variant_table()
    L += ['## 4. Cross-method recovery, AUROC and top-50 enrichment together', '',
          '### 4a. Archived fits at study-group level (the regression-gate population)', '',
          'Per study group, truth = injected S u the group\'s own set; both IG baselines averaged within a fit; delta in {1, 2}, '
          'responder fraction 1.0. MC-ContrastiveVI group ranking = mean of its shared axis and the group\'s own drug axis. '
          'Paired one-sided Wilcoxon (method > MC-ContrastiveVI) on matching (configuration, seed, group); P_expected is the value '
          'supplied as the rough expectation.', '', md_table(pt)]
    if len(pt):
        L += ['Ties (zero differences) are dropped by the default Wilcoxon zero method; the contrastiveVI enrichment test has ties '
              '(see n_ties).', '']
    L += ['### 4b. MC-ContrastiveVI readout robustness', '',
          'The same MC-ContrastiveVI fits scored six ways (norm-concat = elementwise sqrt(shared^2 + drug^2) of the baseline-averaged '
          'mean|IG|, an approximation of the concatenated-latent readout of scripts/44; rank-average = mean of the two axes\' ranks). '
          'Expectation supplied: AUROC between 0.46 and 0.49, enrichment 0.75.', '', md_table(mv)]
    if len(mv):
        okv = bool(mv.auroc_median.between(0.46, 0.49).all() and (mv.enrichment_median.sub(0.75).abs() < 0.006).all())
        L += [f'All six readouts within the expected range: {"yes" if okv else "NO"}.', '']
    L += ['### 4c. All fits, fit-level (new experiments included)', '',
          'Unit = one (method, configuration, seed) with groups averaged. Configurations: non-null scenarios, responder fraction 1.0. '
          'Kruskal-Wallis across methods, then a paired two-sided Wilcoxon signed-rank of each method against MC-ContrastiveVI on matching '
          '(configuration, seed); Holm-adjusted over the methods of that table. AUROC and top-50 enrichment side by side.', '']
    for deltas, title in (([1.0, 2.0], 'delta in {1, 2}'), ([0.25, 0.5], 'delta in {0.25, 0.5}')):
        L += [f'#### {title}', '']
        for metric, label in (('recovery_auroc', 'AUROC'), ('recovery_top50_enrichment', 'top-50 enrichment')):
            cm = cross_method(master, deltas, metric)
            L += [f'**{label}**', '', md_table(cm['means'])]
            if 'H' in cm:
                L += [f'Kruskal-Wallis: H = {cm["H"]:.3f}, df = {cm["df"]}, P = {fmt_p(cm["P"])}; n per method: ' +
                      ', '.join(f'{k} {v}' for k, v in cm['n'].items()) + '.', '']
            else:
                L += ['Kruskal-Wallis not computable (fewer than 2 methods with >= 2 units).', '']
            L += [md_table(cm['wilcoxon'])]

    caveats = []
    caveats += notes
    exp_n = master.groupby('experiment').size().to_dict()
    for exp, want in (('exp1', 54), ('exp2', 30), ('exp7', 12), ('exp3', 54), ('exp4', 36)):
        have = exp_n.get(exp, 0)
        if have < want:
            caveats.append(f'{exp}: {have} of {want} planned fits/jobs are in master.csv (exp3 counts one row per contrastiveVI study and one per multiGroupVI fit).')
    fl = pd.DataFrame(flags, columns=['experiment', 'not_converged', 'ig_quadrature'])
    loud = []
    if len(fl):
        e1 = fl[fl.experiment == 'exp1']
        if len(e1) and e1.not_converged.sum():
            cap = int(master[master.experiment == 'exp1'].max_epochs.max())
            loud.append(f'**WARNING: {int(e1.not_converged.sum())} of {len(e1)} exp1 fits did NOT converge (converged = False: they ran the full {cap}-epoch cap '
                        'without early stopping firing). For those fits the free-bits dose-response is confounded with training budget; see section 1b.**')
        for exp, g in fl.groupby('experiment'):
            if exp != 'exp1' and g.not_converged.sum():
                label = 'archived semi-synthetic' if exp == 'archive_semisynthetic' else exp
                loud.append(f'**{int(g.not_converged.sum())} of {len(g)} {label} fits did not converge (ran the full epoch cap; early stopping never fired).**')
            if g.ig_quadrature.sum():
                caveats.append(f'{exp}: {int(g.ig_quadrature.sum())} of {len(g)} fits have an IG completeness error >= 0.01 even at 4,096 quadrature points.')
    caveats.append('exp3 contrastiveVI / multiGroupVI and exp4 scVI train for a fixed number of epochs with no early stopping; convergence is not assessed.')
    n_dr = int(pooled.n.iloc[0]) if len(pooled) else 0
    if n_dr:
        caveats.append(f'Section 1 pools n = {n_dr} injected fits ({n_dr // 6} per free_bits level, 3 seeds x 2 configurations) in one Spearman test; '
                       f'with 3 seeds per cell the per-configuration rows (n = 18) detect only |rho| of about {np.tanh(2.8 / np.sqrt(18 - 3)):.2f} or more at 80% power. '
                       'Seeds within a configuration share data, and the two configurations differ in effect size, so pooled P-values overstate independence.')
    caveats.append('free_bits = 0 in exp1 repeats archived (config, seed) fits; the two differ only by non-bitwise-reproducible training, so they are not independent replicates and are never pooled in the tests above.')
    caveats.append('PCA32+Wilcoxon is deterministic, so its three seeds are identical; scVI and PCA are scored on training cells with a differential-expression ranking, not IG on test cells; their AUROC is against injected genes only.')
    caveats.append('Recovery is defined per drug group; the MC-ContrastiveVI shared axis (where free bits act most directly) is excluded from `recovery_auroc` and reported in `recovery_auroc_shared_axis`.')
    caveats.append('Probe accuracy is a different task for MC-ContrastiveVI (3-class condition, ungated drug block) than for contrastiveVI / multiGroupVI (within-study treated vs control); do not compare across methods.')
    caveats.append('exp6 uses assumed treated test-cell counts (358 + 357) and iid Gaussian noise; the analytic dimension-split result does not depend on them, the spread does.')
    caveats.append('Many tests are reported without family-wise correction except the Holm adjustment inside each section-4 table.')
    for d in sorted(out_dir.glob('**/progress.log')) + sorted(exp_root.glob('*/progress.log')):
        txt = d.read_text().splitlines()
        fails = [t for t in txt if '\tFAIL\t' in t]
        if fails:
            caveats.append(f'{len(fails)} failed job(s) logged in {d}.')
    L += ['## 5. Caveats: underpowered or not converged', ''] + [f'- {c}' for c in loud + list(dict.fromkeys(caveats))] + ['']
    if loud:
        L[2:2] = [c for c in loud if c.startswith('**WARNING')] + ['']
    return '\n'.join(L)


# ------------------------------------------------------------------------------------------ main
def build_master(args):
    notes = []
    gsp = (args.out / 'smoke_data' / 'D1_gene_sets.csv') if SMOKE else args.root / 'analysis' / 'peerreview_20261003' / 'D1_gene_sets.csv'
    gs = pd.read_csv(gsp) if gsp.exists() else None
    if gs is None:
        notes.append(f'{gsp} not found: AUROC/AP of archived fits cannot be computed and archived R4 enrichment uses 40-gene sets.')
    sizes = {s: int((gs.set == s).sum()) for s in ('S', 'A', 'B')} if gs is not None else {'S': 40, 'A': 40, 'B': 40}
    rows, flags = [], []
    if not SMOKE:
        r, f = load_archive_mccvi(args.root, gs, notes); rows += r; flags += f
    r, f = load_new_mccvi(args.exp_root, notes); rows += r; flags += f
    if not SMOKE:
        rows += load_archive_external(args.root, gs, notes)
    rows += load_exp3(args.exp_root, gs, notes)
    rows += load_exp4(args.exp_root, notes)
    master = pd.DataFrame(rows, columns=COLUMNS)
    return master, flags, notes


def regression_check():
    """Compare group-level medians of the archived fits with the expected values. Returns (ok, table)."""
    rows = []
    g = pd.DataFrame(GROUP_ROWS)
    if not g.empty:
        g = g[g.delta.isin([1.0, 2.0]) & (g.responder_fraction == 1.0)]
    for method, (e_auc, e_enr) in EXPECTED.items():
        d = g[g.method == method] if not g.empty else g
        auc = float(d.auroc.median()) if len(d) else NAN
        enr = float(d.enrichment.median()) if len(d) else NAN
        ok = bool(len(d)) and abs(round(auc, 3) - e_auc) <= 0.0011 and abs(round(enr, 2) - e_enr) <= 0.006
        rows.append(dict(method=method, n_group_rows=len(d), auroc_median=round(auc, 3), auroc_expected=e_auc,
                         enrichment_median=round(enr, 3), enrichment_expected=e_enr, match=ok))
    t = pd.DataFrame(rows)
    return bool(t.match.all()), t


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out', type=Path, default=R / 'analysis' / 'experiments')
    ap.add_argument('--root', type=Path, default=R, help='repo root holding the archived runs/ and analysis/')
    ap.add_argument('--exp-root', type=Path, default=R / 'runs' / 'experiments')
    ap.add_argument('--smoke', action='store_true', help='read the outputs of run_all.py --smoke (synthetic, no archived fits)')
    ap.add_argument('--no-regression-gate', action='store_true',
                    help='write STATS.md even if the scoring does not reproduce the expected archived medians (not for publication)')
    args = ap.parse_args(argv)
    global SMOKE
    SMOKE = args.smoke
    if SMOKE:
        args.out = args.out / 'smoke'
    args.out.mkdir(parents=True, exist_ok=True)
    master, flags, notes = build_master(args)
    master.to_csv(args.out / 'master.csv', index=False)
    print(f'master.csv: {len(master)} fits')
    print(master.groupby(['experiment', 'method']).size().to_string())
    for n in notes:
        print('note:', n)
    if not SMOKE:
        ok, table = regression_check()
        print('\nregression check (archived fits, delta 1 and 2, responder fraction 1.0, group-level medians, S u own set, baselines averaged):')
        print(table.to_string(index=False))
        if not ok and not args.no_regression_gate:
            print('\nREGRESSION CHECK FAILED: the scoring does not reproduce the expected archived medians. STATS.md was NOT written. '
                  f'master.csv was written to {args.out / "master.csv"}. Re-run with --no-regression-gate only to inspect.', file=sys.stderr)
            return 2
    (args.out / 'STATS.md').write_text(build_stats(master, flags, notes, args.out, args.exp_root, args))
    print(f'wrote {args.out / "master.csv"} and {args.out / "STATS.md"}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
