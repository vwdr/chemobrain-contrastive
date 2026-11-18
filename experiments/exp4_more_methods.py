"""Exp 4: two more methods on the same semi-synthetic recovery benchmark: scVI and PCA(32)+Wilcoxon.

delta in {1, 2} x scenario {shared, specific, mixed} (responder fraction 1.0) x 3 seeds, i.e. the injected inputs of
scripts/41 (reused from runs/peerreview_20261004/r4/ when present, otherwise built by scripts/41's own main()).

Scoring is the benchmark of scripts/42-44: for each study group (cisplatin = GSE216146 -> S u A, doxorubicin =
GSE271055 -> S u B) a per-gene importance vector is ranked and scored with src.peerreview.ig_generic.recovery
(top-50 hits vs the union, random expectation, per-set hits) plus AUROC / average precision against the same union.
Neither method has an integrated-gradients readout, so each uses its natural gene-level ranking, computed on the
group's training cells (the cells the model was fitted on):
  * scVI (n_hidden 128, n_latent 32, 2 layers, NB, study as batch, as in scripts/20 and 29; 100 epochs): model.
    differential_expression pseudo-treated vs pseudo-control within the study; rankings |lfc_mean| and |bayes_factor|.
  * PCA(32)+Wilcoxon: PCA(32) fitted on all training cells (log-normalized X), cells of the study reconstructed
    through the 32 components, Wilcoxon rank-sum (scanpy) pseudo-treated vs pseudo-control, ranked by |z|. The same
    Wilcoxon on the unreconstructed X is stored as a reference ranking ('wilcoxon_raw').
PCA and Wilcoxon are deterministic, so their three seeds are identical (they are kept so the table is rectangular).

Usage:
  python experiments/exp4_more_methods.py --out runs/experiments/exp4 [--seeds 3] [--resume] [--threads=3]
  python experiments/exp4_more_methods.py --smoke
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from experiments import common  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

DELTAS = [1.0, 2.0]
SCENARIOS = ['shared', 'specific', 'mixed']
METHODS = ['scvi', 'pca_wilcoxon']
GROUPS = (('cisplatin', 0, ('S', 'A')), ('doxorubicin', 1, ('S', 'B')))   # label, study code b, target sets


def wilcoxon_scores(X, treated, genes):
    """|z| of a Wilcoxon rank-sum test (scanpy) pseudo-treated vs pseudo-control for every gene."""
    import anndata as ad
    import scanpy as sc
    a = ad.AnnData(X=np.asarray(X, dtype=np.float32), var=pd.DataFrame(index=list(genes)))
    a.obs['pseudo'] = pd.Categorical(np.where(treated, 'treated', 'control'), categories=['control', 'treated'])
    sc.tl.rank_genes_groups(a, groupby='pseudo', groups=['treated'], reference='control', method='wilcoxon',
                            n_genes=a.n_vars, use_raw=False)
    df = sc.get.rank_genes_groups_df(a, group='treated').set_index('names')
    return df.scores.reindex(list(genes)).abs().to_numpy()


def score_rows(rankings, genes, label, sets, SETS, method, cfg, seed, recovery, truth):
    """top-50 hits as in scripts/44 (union of the group's sets, whether or not injected); AUROC/AP against the
    sets actually injected in this scenario (``truth``, None when the scenario injects nothing for the group)."""
    rows = []
    targets = {s: SETS[s] for s in sets}
    for name, vec in rankings.items():
        rec, top = recovery(np.asarray(vec, dtype=float), genes, targets)
        m = common.rank_metrics(vec, genes, truth) if truth else {}
        rows.append(dict(method=method, config=cfg, seed=seed, group=label, ranking=name, **rec,
                         n_injected_positive=m.get('n_positive', 0), auroc=m.get('auroc', float('nan')),
                         average_precision=m.get('average_precision', float('nan')),
                         average_precision_baseline=m.get('average_precision_baseline', float('nan'))))
    return rows


def main(argv=None):
    args = common.parse(common.make_parser(__doc__, default_seeds=3, default_out='runs/experiments/exp4'), argv)
    H = common.Harness(); H.prepare(args); ss = H.ss
    from src.peerreview.ig_generic import recovery
    gs = H.gene_sets(); SETS = {s: set(gs[gs.set == s].gene) for s in ('S', 'A', 'B')}
    jobs = [(f'{m}_{ss.config_name(sc_, dl, 1.0)}_s{seed}', dict(method=m, cfg=ss.config_name(sc_, dl, 1.0), seed=seed))
            for m in METHODS for sc_ in SCENARIOS for dl in DELTAS for seed in range(args.seeds)]
    jobs = common.select_jobs(jobs, args)
    runner = common.Runner('exp4', args)
    inputs_dir = args.out / 'inputs'

    def fn(job_id, method, cfg, seed):
        import time
        t0 = time.time()
        z = np.load(H.ensure_inputs(cfg, args, inputs_dir) / f'inputs_{cfg}.npz')
        C, X, d, b, tr, genes, treated = z['C'], z['X'], z['d'], z['b'], z['train'], z['genes'], z['treated']
        scen, delta, rf = ss.parse_config(cfg)
        rows = []
        truth = common.injected_targets(cfg, gs, ss)['with_shared']
        if method == 'pca_wilcoxon':
            from sklearn.decomposition import PCA
            k = min(32, X.shape[1], len(tr))
            pca = PCA(k, svd_solver='full', random_state=seed).fit(X[tr])
            for label, sb, sets in GROUPS:
                trs = tr[b[tr] == sb]
                recon = pca.inverse_transform(pca.transform(X[trs]))
                rows += score_rows({'wilcoxon_pca32_recon': wilcoxon_scores(recon, treated[trs], genes),
                                    'wilcoxon_raw': wilcoxon_scores(X[trs], treated[trs], genes)},
                                   genes, label, sets, SETS, method, cfg, seed, recovery, truth.get(label))
            extra = dict(n_components=k, explained_variance_ratio_sum=float(pca.explained_variance_ratio_.sum()))
        else:
            import anndata as ad
            import scvi
            import torch
            scvi.settings.num_threads = common.N_THREADS
            scvi.settings.seed = seed; torch.manual_seed(seed); np.random.seed(seed)
            obs = pd.DataFrame({'study': np.where(b[tr] == 1, 'GSE271055', 'GSE216146'),
                                'pseudo': np.where(treated[tr], 'treated', 'control')}, index=[f'c{i}' for i in tr])
            atr = ad.AnnData(X=np.rint(C[tr]).astype(np.int64), obs=obs, var=pd.DataFrame(index=list(genes)))
            scvi.model.SCVI.setup_anndata(atr, batch_key='study')
            model = scvi.model.SCVI(atr, n_hidden=128, n_latent=32, n_layers=2, dropout_rate=0.1, dispersion='gene',
                                    gene_likelihood='nb', use_observed_lib_size=True)
            model.train(max_epochs=args.epochs, train_size=1.0, validation_size=0.0, early_stopping=False,
                        accelerator='cpu', devices=1, batch_size=256, enable_progress_bar=False)
            for label, sb, sets in GROUPS:
                in_study = (b[tr] == sb)
                i1 = in_study & treated[tr]; i2 = in_study & ~treated[tr]
                de = model.differential_expression(adata=atr, idx1=i1, idx2=i2, mode='change', delta=0.25,
                                                   batch_correction=False, silent=True).reindex(list(genes))
                rows += score_rows({'abs_lfc_mean': de.lfc_mean.abs().to_numpy(),
                                    'abs_bayes_factor': de.bayes_factor.abs().to_numpy()},
                                   genes, label, sets, SETS, method, cfg, seed, recovery, truth.get(label))
            extra = dict(epochs=args.epochs)
        return dict(method=method, config=cfg, scenario=scen, delta=delta, responder_fraction=rf, seed=seed,
                    n_train_cells=int(len(tr)), smoke=bool(args.smoke), seconds=time.time() - t0, recovery=rows, **extra)

    failed = runner.run(jobs, fn, resume=args.resume)
    common.write_summary(runner, [j for j, _ in jobs])
    flat = [dict(r, scenario=res['scenario'], delta=res['delta']) for res in runner.load_results([j for j, _ in jobs])
            for r in res.get('recovery', [])]
    if flat:
        pd.DataFrame(flat).to_csv(args.out / 'recovery_per_fit.csv', index=False)
    return runner.stats


if __name__ == '__main__':
    sys.exit(1 if main()['failed'] else 0)
