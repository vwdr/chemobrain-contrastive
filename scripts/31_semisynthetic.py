"""Phase 4: semi-synthetic positive control for MC-ContrastiveVI (binomial thinning of control cells).

Design (fixed in advance, CLAUDE_CODE_PROMPT.md Phase 4):
  * base cells: nonrescue controls only (GSE216146 PBS = source code PN; GSE271055 control = CNT);
    within each library a 50/50 random split into pseudo-control / pseudo-treated (seed 20261003);
    pseudo-treated cells take the cisplatin code (GSE216146) or doxorubicin code (GSE271055);
  * genes: the canonical 1,500 genes;
  * gene sets S, A, B (40 genes each, disjoint), 4 genes per decile of mean log-normalized base-cell
    expression per set, 2 up + 2 down per decile per set (seed 20261003);
  * injection: binomial thinning (p = 2**-delta) of raw counts; up genes thinned in all cells of the
    relevant study(ies) except responders, down genes thinned in responders; S in both studies, A in
    GSE216146, B in GSE271055; responders = random subset of pseudo-treated cells of the study at the
    responder fraction (seed 20261003); thinning draws use default_rng(20261003) per configuration;
    canonical normalization (1e4 over all genes, log1p);
  * scenarios: null; shared (S); specific (A+B); mixed (S+A+B); delta in {0.25,0.5,1,2} x responder
    fraction in {1.0, 0.3}; seeds 0-2;
  * split: canonical procedure (80/10/10 stratified by study x label, 1,500 training cells/stratum, seed 1729);
  * label-exchange reference: null and shared_d1_rf1.0, 100 within-library permutations of the
    pseudo labels (seed 20261003), seed-0 fits, pooled F_sh.

Usage:
  python scripts/31_semisynthetic.py --prepare
  python scripts/31_semisynthetic.py --fit --config shared_d1_rf1.0 --seed 0 [--threads=3]
  python scripts/31_semisynthetic.py --perm --config null --perm-index 0
  python scripts/31_semisynthetic.py --list-jobs grid|perm
  python scripts/31_semisynthetic.py --score
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.peerreview.threads import configure, strip_threads_arg  # noqa: E402
N_THREADS = configure(3)

import argparse  # noqa: E402
import itertools  # noqa: E402
import json  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import torch  # noqa: E402

from src.peerreview import core, fitjob, semisynth  # noqa: E402

torch.set_num_threads(N_THREADS)
SEED = 20261003
RUN = core.PR_RUN / 'semisynthetic'; OUT = core.PR_OUT
DELTAS = [0.25, 0.5, 1.0, 2.0]; FRACTIONS = [1.0, 0.3]; SCENARIOS = ['shared', 'specific', 'mixed']
SETS_BY_SCENARIO = {'null': [], 'shared': ['S'], 'specific': ['A', 'B'], 'mixed': ['S', 'A', 'B']}
N_PERM = 100
PERM_CONFIGS = ['null', 'shared_d1_rf1.0']


def config_name(scen, delta=None, rf=None):
    return 'null' if scen == 'null' else f'{scen}_d{delta:g}_rf{rf:.1f}'


def parse_config(name):
    if name == 'null':
        return 'null', 0.0, None
    scen, d, r = name.split('_')
    return scen, float(d[1:]), float(r[2:])


def all_configs(reductions=()):
    cfgs = ['null']
    for scen, delta, rf in itertools.product(SCENARIOS, DELTAS, FRACTIONS):
        if rf == 0.3 and ('drop_rf03_all' in reductions or ('drop_rf03_extreme' in reductions and delta in (0.25, 2.0))):
            continue
        cfgs.append(config_name(scen, delta, rf))
    return cfgs


# ------------------------------------------------------------------ preparation
def prepare():
    import anndata as ad
    import scanpy as sc
    RUN.mkdir(parents=True, exist_ok=True)
    if (RUN / 'base.npz').exists():
        print('already prepared'); return
    a = ad.read_h5ad(core.RAW_H5AD)
    o = a.obs
    base = np.where(~o.rescue.astype(bool).to_numpy() & (((o.study == 'GSE216146') & (o.source_code == 'PN')) |
                                                        ((o.study == 'GSE271055') & (o.source_code == 'CNT'))).to_numpy())[0]
    genes = pd.read_csv(core.CORR_OUT / 'benchmark_gene_universe.csv')['gene'].astype(str).tolist()
    sub = a[base]
    C = sub[:, genes].X.toarray().astype(np.float32)
    L = np.asarray(sub.X.sum(1)).ravel().astype(np.float64)
    obs = sub.obs[['study', 'sample_id', 'source_code', 'source_cell_type']].astype(str).reset_index()
    # normalization check against scanpy on the uninjected base cells
    chk = sub.copy(); sc.pp.normalize_total(chk, target_sum=1e4); sc.pp.log1p(chk)
    X0_scanpy = chk[:, genes].X.toarray().astype(np.float32)
    X0 = semisynth.log_normalize(C, L)
    norm_diff = float(np.abs(X0 - X0_scanpy).max())
    # 50/50 split within library
    rng = np.random.default_rng(SEED); treated = np.zeros(len(base), dtype=bool)
    for lib in sorted(obs.sample_id.unique()):
        idx = np.where(obs.sample_id.to_numpy() == lib)[0]
        perm = rng.permutation(idx); treated[perm[: len(idx) // 2]] = True
    b = (obs.study.to_numpy() == 'GSE271055').astype(np.int64)
    d = np.where(treated, np.where(b == 1, 1, 2), 0).astype(np.int64)
    # gene sets
    mean_expr = X0.mean(0)
    dec = pd.qcut(pd.Series(mean_expr), 10, labels=False).to_numpy()
    rng_g = np.random.default_rng(SEED); rows = []
    for k in range(10):
        gi = np.where(dec == k)[0]
        perm = rng_g.permutation(gi)
        for si, set_name in enumerate(['S', 'A', 'B']):
            chosen = perm[si * 4:(si + 1) * 4]
            for j, g in enumerate(chosen):
                rows.append(dict(set=set_name, gene=genes[g], gene_index=int(g), decile=int(k), direction='up' if j < 2 else 'down',
                                 mean_log_norm_expression_base=float(mean_expr[g])))
    gs = pd.DataFrame(rows)
    assert gs.gene.is_unique and (gs.groupby('set').size() == 40).all()
    assert (gs.groupby(['set', 'direction']).size() == 20).all()
    gs.to_csv(OUT / 'D1_gene_sets.csv', index=False)
    # responders per fraction
    resp = {}
    for rf in FRACTIONS:
        rng_r = np.random.default_rng(SEED); m = np.zeros(len(base), dtype=bool)
        for study_b in (0, 1):
            pt = np.where(treated & (b == study_b))[0]
            m[rng_r.choice(pt, int(round(rf * len(pt))), replace=False)] = True
        resp[rf] = m
    # canonical split on (study x pseudo-label)
    strata = (obs.study.to_numpy().astype(str) + '_' + np.where(treated, 'treated', 'control'))
    tr, va, te = core.canonical_split(strata, np.arange(len(base)))
    np.savez_compressed(RUN / 'base.npz', C=C, L=L, base_index=base, d=d, b=b, treated=treated, train=tr, val=va, test=te,
                        resp_1p0=resp[1.0], resp_0p3=resp[0.3], genes=np.asarray(genes, dtype='U'))
    obs.assign(pseudo_treated=treated, d=d, responder_rf1p0=resp[1.0], responder_rf0p3=resp[0.3]).to_csv(RUN / 'base_cells.csv', index=False)
    split = np.full(len(base), 'unused', dtype=object); split[tr] = 'train'; split[va] = 'validation'; split[te] = 'test'
    cnt = (obs.assign(pseudo=np.where(treated, 'pseudo_treated', 'pseudo_control'), split=split)
           .groupby(['study', 'sample_id', 'pseudo', 'split']).size().rename('n_cells').reset_index())
    cnt.to_csv(OUT / 'D1_base_cell_counts.csv', index=False)
    rcnt = []
    for rf in FRACTIONS:
        for study_b, study in ((0, 'GSE216146'), (1, 'GSE271055')):
            pt = treated & (b == study_b)
            rcnt.append(dict(responder_fraction=rf, study=study, n_pseudo_treated=int(pt.sum()), n_responders=int((resp[rf] & pt).sum()),
                             n_responders_test=int((resp[rf] & pt)[te].sum()), n_pseudo_treated_test=int(pt[te].sum())))
    pd.DataFrame(rcnt).to_csv(OUT / 'D1_responder_counts.csv', index=False)
    core.write_json(OUT / 'D1_design.json', {
        'seed': SEED, 'n_base_cells': int(len(base)), 'n_train': int(len(tr)), 'n_val': int(len(va)), 'n_test': int(len(te)),
        'max_abs_diff_manual_vs_scanpy_normalization_base': norm_diff, 'deltas': DELTAS, 'responder_fractions': FRACTIONS,
        'scenarios': SETS_BY_SCENARIO, 'n_permutations': N_PERM, 'permutation_configs': PERM_CONFIGS,
        'rng_streams': 'separate numpy default_rng(20261003) generators for: library split, gene-set draw, responders '
                       '(re-created per responder fraction), thinning draws (re-created per configuration), label permutations',
        'decile_definition': 'pandas.qcut of mean log-normalized expression (canonical normalization) over base cells, 10 bins'})
    print('prepared', len(base), len(tr), len(va), len(te), 'norm diff', norm_diff, flush=True)


def load_base():
    z = np.load(RUN / 'base.npz'); return {k: z[k] for k in z.files}


def build_inputs(cfg, base=None):
    base = base or load_base()
    scen, delta, rf = parse_config(cfg)
    gs = pd.read_csv(OUT / 'D1_gene_sets.csv')
    C = base['C']; b = base['b']
    progs = []
    if scen != 'null':
        resp = base['resp_1p0'] if rf == 1.0 else base['resp_0p3']
        for s in SETS_BY_SCENARIO[scen]:
            g = gs[gs.set == s]
            mask = np.ones(len(b), bool) if s == 'S' else (b == (0 if s == 'A' else 1))
            progs.append(dict(name=s, up_cols=g[g.direction == 'up'].gene_index.to_numpy(),
                              down_cols=g[g.direction == 'down'].gene_index.to_numpy(), study_mask=mask, responders=resp & mask,
                              delta=delta))
    Cn = semisynth.inject(C, progs, np.random.default_rng(SEED)) if progs else C
    Lnew = base['L'] - (C.astype(np.float64) - Cn.astype(np.float64)).sum(1)
    X = semisynth.log_normalize(Cn, Lnew)
    return dict(X=X, d=base['d'], b=b, train=base['train'], val=base['val'], test=base['test'], genes=base['genes']), progs


def true_shared_fraction(cfg, base):
    scen, delta, rf = parse_config(cfg)
    if scen == 'null':
        return dict(true_shared_energy_fraction=np.nan, true_shared_trace_fraction=np.nan)
    gs = pd.read_csv(OUT / 'D1_gene_sets.csv'); te = base['test']; d = base['d']; b = base['b']
    T = te[d[te] > 0]; resp = (base['resp_1p0'] if rf == 1.0 else base['resp_0p3'])
    def e_block(s):
        if s not in SETS_BY_SCENARIO[scen]:
            return np.zeros((len(T), 40))
        g = gs[gs.set == s]; sign = np.where(g.direction.to_numpy() == 'up', 1.0, -1.0) * delta * np.log(2)
        active = resp[T] & (np.ones(len(T), bool) if s == 'S' else (b[T] == (0 if s == 'A' else 1)))
        return active[:, None] * sign[None, :]
    eS, eA, eB = e_block('S'), e_block('A'), e_block('B')
    num = (eS ** 2).sum(); den = num + (eA ** 2).sum() + (eB ** 2).sum()
    tS = np.var(eS, axis=0).sum(); tD = np.var(np.hstack([eA, eB]), axis=0).sum()
    return dict(true_shared_energy_fraction=float(num / den) if den > 0 else np.nan,
                true_shared_trace_fraction=float(tS / (tS + tD)) if (tS + tD) > 0 else np.nan)


def fit(cfg, seed):
    inp, _ = build_inputs(cfg)
    fitjob.run_fit(inp, seed, RUN / cfg, f'{cfg}_s{seed}', save_checkpoint=True, save_latents=True,
                   extra={'config': cfg, 'threads': N_THREADS})


def perm_labels(base, k):
    rng = np.random.default_rng(SEED); libs = pd.read_csv(RUN / 'base_cells.csv').sample_id.to_numpy()
    d0 = base['d']; out = None
    for i in range(k + 1):
        d = d0.copy()
        for lib in sorted(np.unique(libs)):
            idx = np.where(libs == lib)[0]; d[idx] = rng.permutation(d0[idx])
        out = d
    return out


def perm(cfg, k):
    base = load_base(); inp, _ = build_inputs(cfg, base)
    inp['d'] = perm_labels(base, k)
    fitjob.run_fit(inp, 0, RUN / f'perm_{cfg}', f'perm{k:03d}', do_ig=False, save_checkpoint=False, save_latents=False,
                   extra={'config': cfg, 'perm_index': k, 'threads': N_THREADS})


# ------------------------------------------------------------------ scoring
def ig_scores(top_by_axis, gs):
    S = set(gs[gs.set == 'S'].gene); A = set(gs[gs.set == 'A'].gene); B = set(gs[gs.set == 'B'].gene)
    sh, dox, cis = top_by_axis['shared'], top_by_axis['doxorubicin'], top_by_axis['cisplatin']
    r = {}
    for nm, lst, tgt in [('shared_vs_S', sh, S), ('cis_vs_A', cis, A), ('dox_vs_B', dox, B)]:
        k = len(set(lst) & tgt); r[nm + '_hits'] = k
        r[nm + '_precision'] = k / len(lst) if len(lst) else np.nan; r[nm + '_recall'] = k / len(tgt)
    r['shared_top_n_AB_genes'] = len(set(sh) & (A | B)); r['cis_top_n_S_genes'] = len(set(cis) & S)
    r['dox_top_n_S_genes'] = len(set(dox) & S); r['cis_top_n_B_genes'] = len(set(cis) & B); r['dox_top_n_A_genes'] = len(set(dox) & A)
    r['list_size_shared'] = len(sh); r['list_size_dox'] = len(dox); r['list_size_cis'] = len(cis)
    return r


def score(reductions=()):
    base = load_base(); gs = pd.read_csv(OUT / 'D1_gene_sets.csv'); te = base['test']; genes = base['genes']
    fits = []; igrows = []; cons_rows = []; cka_rows = []
    for cfg in all_configs(reductions):
        tsf = true_shared_fraction(cfg, base)
        names = [f'{cfg}_s{s}' for s in (0, 1, 2)]
        if not all((RUN / cfg / f'{n}_done.json').exists() for n in names):
            print('incomplete', cfg); continue
        for s, n in zip((0, 1, 2), names):
            m = json.loads((RUN / cfg / f'{n}_done.json').read_text()); m.update(tsf); fits.append(m)
            df = pd.read_csv(RUN / cfg / f'{n}_ig_summary.csv.gz')
            for base_name, g in df.groupby('baseline'):
                top = {core.AXIS_NAME[ax]: list(core.top_genes(gg.mean_abs_ig.to_numpy(), gg.gene.to_numpy())) for ax, gg in g.groupby('axis')}
                igrows.append(dict(config=cfg, seed=s, baseline=base_name, **ig_scores(top, gs)))
        rankings, gorder = fitjob.load_ig_rankings(RUN / cfg, names)
        cons, _ = core.consensus_lists(rankings, gorder)
        cons_rows.append(dict(config=cfg, **ig_scores(cons, gs), shared_list='|'.join(cons['shared']),
                              dox_list='|'.join(cons['doxorubicin']), cis_list='|'.join(cons['cisplatin'])))
        lats = {s: np.load(RUN / cfg / f'{cfg}_s{s}_latents.npz') for s in (0, 1, 2)}
        for a_, c_ in itertools.combinations((0, 1, 2), 2):
            for key in ('bg', 'shared', 'drug'):
                cka_rows.append(dict(config=cfg, seed_a=a_, seed_b=c_, block=key, test_cka=core.cka(lats[a_][key][te], lats[c_][key][te])))
    fits = pd.DataFrame(fits); igr = pd.DataFrame(igrows); cons = pd.DataFrame(cons_rows); ck = pd.DataFrame(cka_rows)
    for col, fn in [('scenario', lambda c: parse_config(c)[0]), ('delta', lambda c: parse_config(c)[1]), ('responder_fraction', lambda c: parse_config(c)[2])]:
        for df in (fits, igr, cons, ck):
            df.insert(1, col, df.config.map(fn))
    fits.to_csv(OUT / 'D1_fit_metrics.csv', index=False); igr.to_csv(OUT / 'D1_ig_per_fit.csv', index=False)
    cons.to_csv(OUT / 'D1_ig_consensus.csv', index=False); ck.to_csv(OUT / 'D1_cross_seed_cka.csv', index=False)
    keys = ['config', 'scenario', 'delta', 'responder_fraction']
    summ = fits.groupby(keys, dropna=False, sort=False).agg(
        true_shared_energy_fraction=('true_shared_energy_fraction', 'first'), true_shared_trace_fraction=('true_shared_trace_fraction', 'first'),
        F_sh_mean=('shared_fraction', 'mean'), F_sh_sd=('shared_fraction', 'std'), F_sh_min=('shared_fraction', 'min'), F_sh_max=('shared_fraction', 'max'),
        F_sh_dox_mean=('shared_fraction_doxorubicin', 'mean'), F_sh_cis_mean=('shared_fraction_cisplatin', 'mean'),
        test_mse_mean=('test_mse', 'mean'), shared_AU_mean=('shared_active_units_var_gt_0p01', 'mean'),
        dox_AU_mean=('dox_active_units_var_gt_0p01', 'mean'), cis_AU_mean=('cis_active_units_var_gt_0p01', 'mean'),
        shared_KL_mean=('shared_mean_total_kl_nats', 'mean'), dox_KL_mean=('dox_mean_total_kl_nats', 'mean'),
        cis_KL_mean=('cis_mean_total_kl_nats', 'mean')).reset_index()
    summ['F_sh_range'] = summ.F_sh_max - summ.F_sh_min
    ckw = ck.groupby(['config', 'block']).test_cka.mean().unstack().add_prefix('mean_cross_seed_cka_').reset_index()
    summ = summ.merge(ckw, on='config', how='left')
    igs = igr.groupby(['config', 'baseline']).mean(numeric_only=True).reset_index()
    summ.to_csv(OUT / 'D1_config_summary.csv', index=False); igs.to_csv(OUT / 'D1_ig_per_config_mean.csv', index=False)
    # permutations
    prow = []
    for cfg in PERM_CONFIGS:
        pd_ = RUN / f'perm_{cfg}'
        vals = [json.loads(p.read_text())['shared_fraction'] for p in sorted(pd_.glob('perm*_done.json'))] if pd_.exists() else []
        obs = fits[(fits.config == cfg) & (fits.seed == 0)].shared_fraction
        if len(vals) and len(obs):
            o = float(obs.iloc[0]); v = np.asarray(vals)
            prow.append(dict(config=cfg, observed_seed0_F_sh=o, n_permutations=len(v), n_perm_ge_observed=int((v >= o).sum()),
                             p_value=(1 + int((v >= o).sum())) / (1 + len(v)), perm_F_sh_mean=float(v.mean()), perm_F_sh_sd=float(v.std(ddof=1)),
                             perm_F_sh_min=float(v.min()), perm_F_sh_q025=float(np.quantile(v, .025)), perm_F_sh_median=float(np.median(v)),
                             perm_F_sh_q975=float(np.quantile(v, .975)), perm_F_sh_max=float(v.max())))
            pd.DataFrame({'perm_index': range(len(v)), 'F_sh': v}).to_csv(OUT / f'D1_permutation_values_{cfg}.csv', index=False)
    pd.DataFrame(prow).to_csv(OUT / 'D1_permutation_summary.csv', index=False)
    pd.set_option('display.width', 300); pd.set_option('display.max_rows', 200); pd.set_option('display.max_columns', 40)
    print(summ.to_string(index=False)); print(cons.drop(columns=['shared_list', 'dox_list', 'cis_list']).to_string(index=False))
    print(pd.DataFrame(prow).to_string(index=False))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--prepare', action='store_true'); ap.add_argument('--fit', action='store_true'); ap.add_argument('--perm', action='store_true')
    ap.add_argument('--config'); ap.add_argument('--seed', type=int, default=0); ap.add_argument('--perm-index', type=int, default=0)
    ap.add_argument('--list-jobs'); ap.add_argument('--score', action='store_true'); ap.add_argument('--reductions', default='')
    a = ap.parse_args(strip_threads_arg(sys.argv[1:])); red = tuple(x for x in a.reductions.split(',') if x)
    if a.prepare:
        prepare()
    if a.fit:
        fit(a.config, a.seed)
    if a.perm:
        perm(a.config, a.perm_index)
    if a.list_jobs == 'grid':
        for c in all_configs(red):
            for s in (0, 1, 2):
                print(f'fit {c} {s}')
    elif a.list_jobs == 'perm':
        for c in PERM_CONFIGS:
            for k in range(N_PERM):
                print(f'perm {c} {k}')
    if a.score:
        score(red)
