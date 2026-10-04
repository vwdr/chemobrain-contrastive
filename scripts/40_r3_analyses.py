"""R3 analyses on the fresh canonical fits and the R3 comparison fits (outputs: analysis/peerreview_20261004/R3_*).

Parts (each writes its own tables):
  comparison   : held-out MSE, background study / source-cell-type probes, NB shared fraction; NB latents are computed
                 from the NB checkpoints (06_validation.means logic) and saved next to them.
  consensus    : Figure-5 data for the fresh consensus lists (top 12 per block; per-seed mean |IG|, squared-norm target,
                 averaged over zero and control-median baselines; from R3_09_interpretation/complete_corrected_attributions.csv).
  enrichment   : 19_gene_pattern_audit procedure on the fresh consensus lists (and the frozen lists for reference),
                 overlap with the regenerated direction-consistent cisplatin pseudobulk set and the frozen scDisInFact top-100.
  pseudobulk   : committed vs regenerated pseudobulk quantities reported in the manuscript.
Requires: scripts/36 (comparison fits), scripts/39 --script 09/11/16 outputs.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.peerreview.threads import configure  # noqa: E402
configure(2)

import argparse  # noqa: E402
import json  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import torch  # noqa: E402
from sklearn.linear_model import LogisticRegression  # noqa: E402
from sklearn.metrics import balanced_accuracy_score  # noqa: E402

from src.peerreview import core  # noqa: E402

torch.set_num_threads(2)
R = core.R
OUT = R / 'analysis' / 'peerreview_20261004'
COMP = R / 'runs' / 'peerreview_20261004' / 'comparison'
D11 = OUT / 'R3_11_diagnostics'; D09 = OUT / 'R3_09_interpretation'
PB_NEW = R / 'runs' / 'peerreview_20261004' / 'pseudobulk_task4'
IN = core.PR_RUN / 'canonical_inputs'


def comparison():
    z = np.load(IN / 'inputs.npz'); d, b, tr, te = z['d'], z['b'], z['train'], z['test']
    cells = pd.read_csv(IN / 'cells.csv', index_col=0); ct = cells.source_cell_type.astype(str).to_numpy()
    abl = pd.read_csv(D11 / 'ablation_diagnostics.csv')
    rows = []
    for mode in ('full', 'no_hsic', 'no_gating', 'gaussian_vae'):
        for s in (0, 1, 2):
            p = (core.PR_RUN / 'canonical' / f'full_{s}_metrics.json') if mode == 'full' else (COMP / f'{mode}_{s}_metrics.json')
            m = json.loads(p.read_text())
            a = abl[(abl['mode'] == mode) & (abl.seed == s)].iloc[0]
            rows.append(dict(model=mode, seed=s, test_mse=m['test_mse'], bg_study_balanced_accuracy=m['bg_study_balanced_accuracy'],
                             bg_source_cell_type_balanced_accuracy=a.background_cell_type_balanced_accuracy,
                             shared_fraction=m.get('shared_fraction'), best_epoch=m['best_epoch'],
                             source='runs/peerreview_20261003/canonical' if mode == 'full' else 'runs/peerreview_20261004/comparison'))
    # NB: latents from checkpoint (06 means logic, gated as 'full')
    Xt = torch.from_numpy(z['X']); dt = torch.from_numpy(d); bt = torch.from_numpy(b)
    for s in (0, 1, 2):
        m = json.loads((COMP / f'negative_binomial_{s}_metrics.json').read_text())
        c = torch.load(COMP / f'negative_binomial_{s}.pt', map_location='cpu', weights_only=True)
        model = core.mm.MCContrastiveVI(core.mm.MCContrastiveVIConfig(**c['config'])); model.load_state_dict(c['state']); model.eval()
        lat = {k: [] for k in ('bg', 'shared', 'drug', 'raw_shared', 'raw_drug')}
        with torch.no_grad():
            for ix in np.array_split(np.arange(len(Xt)), int(np.ceil(len(Xt) / 512))):
                bg, sh, dr, e = core.means(model, Xt[ix], dt[ix])
                lat['bg'].append(bg.numpy()); lat['shared'].append(sh.numpy()); lat['drug'].append(dr.numpy())
                lat['raw_shared'].append(e['mu_shared'].numpy()); lat['raw_drug'].append(torch.cat(e['mu_drug_list'], 1).numpy())
        lat = {k: np.concatenate(v) for k, v in lat.items()}
        np.savez_compressed(COMP / f'negative_binomial_{s}_latents.npz', **lat)
        st = LogisticRegression(max_iter=500, class_weight='balanced').fit(lat['bg'][tr], b[tr])
        aa = tr[b[tr] == 0]; bb = te[b[te] == 0]
        cc = LogisticRegression(max_iter=800, class_weight='balanced').fit(lat['bg'][aa], ct[aa])
        rows.append(dict(model='negative_binomial', seed=s, test_mse=m['test_mse'],
                         bg_study_balanced_accuracy=float(balanced_accuracy_score(b[te], st.predict(lat['bg'][te]))),
                         bg_source_cell_type_balanced_accuracy=float(balanced_accuracy_score(ct[bb], cc.predict(lat['bg'][bb]))),
                         shared_fraction=m['shared_fraction'], best_epoch=m['best_epoch'], source='runs/peerreview_20261004/comparison'))
    pm = json.loads((COMP / 'pca_metrics.json').read_text())
    rows.append(dict(model='pca32', seed=None, test_mse=pm['test_mse'], bg_study_balanced_accuracy=pm['study_balanced_accuracy'],
                     bg_source_cell_type_balanced_accuracy=None, shared_fraction=None, best_epoch=None, source='runs/peerreview_20261004/comparison'))
    df = pd.DataFrame(rows); df.to_csv(OUT / 'R3_comparison_models.csv', index=False)
    summ = df.groupby('model', sort=False).agg(n=('test_mse', 'size'), test_mse_min=('test_mse', 'min'), test_mse_max=('test_mse', 'max'),
                                                bg_study_min=('bg_study_balanced_accuracy', 'min'), bg_study_max=('bg_study_balanced_accuracy', 'max'),
                                                bg_celltype_min=('bg_source_cell_type_balanced_accuracy', 'min'),
                                                bg_celltype_max=('bg_source_cell_type_balanced_accuracy', 'max'),
                                                shared_fraction_min=('shared_fraction', 'min'), shared_fraction_max=('shared_fraction', 'max')).reset_index()
    summ.to_csv(OUT / 'R3_comparison_models_summary.csv', index=False)
    print(summ.to_string(index=False))


def consensus():
    cons = json.loads((core.PR_OUT / 'A1_fresh_consensus.json').read_text())
    att = pd.read_csv(D09 / 'complete_corrected_attributions.csv')
    att = att[att.target == 'norm_squared']
    axis_of = {'shared': 'shared', 'doxorubicin': 'drug0', 'cisplatin': 'drug1'}
    rows = []
    for block, genes in cons.items():
        g = att[(att.axis == axis_of[block]) & att.gene.isin(genes[:12])]
        per = g.groupby(['gene', 'seed']).mean_abs_ig.mean().reset_index()   # mean over the two baselines
        for rank, gene in enumerate(genes[:12], 1):
            for r in per[per.gene == gene].itertuples():
                rows.append(dict(block=block, consensus_rank=rank, gene=gene, seed=int(r.seed), mean_abs_ig_avg_two_baselines=float(r.mean_abs_ig),
                                 target='norm_squared', baselines='zero|control_median'))
    df = pd.DataFrame(rows); df.to_csv(OUT / 'R3_consensus_top12_figure5_data.csv', index=False)
    print(df.groupby('block').gene.nunique())


def enrichment():
    g19 = core.load_script('gene_pattern_audit_19', '19_gene_pattern_audit.py')
    g19.OUT = PB_NEW   # derive_cisplatin_de_support reads the regenerated paired results here and writes its gene table here
    support = g19.derive_cisplatin_de_support()
    universe = set(pd.read_csv(core.CORR_OUT / 'benchmark_gene_universe.csv')['gene'].astype(str))
    support_u = support & universe
    lists = {('fresh', k): v for k, v in json.loads((core.PR_OUT / 'A1_fresh_consensus.json').read_text()).items()}
    lists.update({('frozen', k): v for k, v in g19.CONSENSUS.items()})
    cross = []; allr = []
    for (ver, axis), genes in lists.items():
        sel = set(genes)
        k_de, p_de, g_de = g19.overlap_test(sel, support_u, len(universe))
        k_sc, p_sc, g_sc = g19.overlap_test(sel, g19.SCD_TOP100, len(universe))
        cross.append(dict(consensus_version=ver, axis=axis, n_consensus_genes=len(sel), n_support_genes_in_universe=len(support_u),
                          n_overlap_cisplatin_DE=k_de, cisplatin_DE_overlap_genes=g_de, hypergeom_p_cisplatin_DE=p_de,
                          n_overlap_scdisinfact_top100=k_sc, scdisinfact_overlap_genes=g_sc, hypergeom_p_scdisinfact_top100=p_sc))
        for lib, fn in g19.LIBRARIES.items():
            for r in g19.enrichment(sel, universe, lib, g19.EVID / fn):
                r.update(consensus_version=ver, axis=axis); allr.append(r)
    cross = pd.DataFrame(cross)
    from statsmodels.stats.multitest import multipletests
    for ver in ('fresh', 'frozen'):
        m = cross.consensus_version == ver
        cross.loc[m, 'q_cisplatin_DE_bh'] = multipletests(cross.loc[m, 'hypergeom_p_cisplatin_DE'], method='fdr_bh')[1]
        cross.loc[m, 'q_scdisinfact_top100_bh'] = multipletests(cross.loc[m, 'hypergeom_p_scdisinfact_top100'], method='fdr_bh')[1]
    cross.to_csv(OUT / 'R3_enrichment_cross_checks.csv', index=False)
    en = pd.DataFrame(allr); en.to_csv(OUT / 'R3_enrichment_all.csv.gz', index=False, compression='gzip')
    summ = (en.assign(sig=en.q_value < 0.05).groupby(['consensus_version', 'axis', 'library'], as_index=False)
            .agg(n_tested=('term', 'size'), n_fdr_005=('sig', 'sum'), minimum_q=('q_value', 'min')))
    summ.to_csv(OUT / 'R3_enrichment_summary.csv', index=False)
    top = en.sort_values(['consensus_version', 'axis', 'library', 'q_value', 'p_value', 'term']).groupby(
        ['consensus_version', 'axis', 'library'], group_keys=False).head(5)
    top.to_csv(OUT / 'R3_enrichment_top_terms.csv', index=False)
    named = en[en.term.isin(['REACTOME_IMMUNE_SYSTEM', 'REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX', 'ZHANG_UTERUS_C5_MACROPHAGE'])]
    named.to_csv(OUT / 'R3_enrichment_named_terms.csv', index=False)
    pd.set_option('display.width', 250)
    print(cross.drop(columns=['cisplatin_DE_overlap_genes', 'scdisinfact_overlap_genes']).to_string(index=False)); print(summ.to_string(index=False))
    print(named[['consensus_version', 'axis', 'term', 'overlap', 'fold_enrichment', 'q_value']].to_string(index=False))


def pseudobulk():
    rows = []
    for label, folder in (('committed', core.CORR_OUT), ('regenerated', PB_NEW)):
        s = pd.read_csv(folder / 'pseudobulk_task4_summary.csv'); rb = pd.read_csv(folder / 'pseudobulk_task4_robustness.csv')
        for c, g in s.groupby('contrast'):
            rows.append(dict(version=label, quantity=f'n_BH_findings_{c}', value=int(g.n_wald_padj_lt_0p05.sum())))
        rows.append(dict(version=label, quantity='pooled_BH_findings_all_contrasts', value=int(s.n_wald_padj_lt_0p05.sum())))
        rows.append(dict(version=label, quantity='n_findings_all3_pairs_same_direction', value=int(rb.n_deseq2_hits_all3_pairs_same_direction.sum())))
        rows.append(dict(version=label, quantity='n_welch_fdr_005_total', value=int(rb.n_welch_fdr_005.sum())))
        hit = rb[rb.n_deseq2_fdr_005 > 0]
        rows.append(dict(version=label, quantity='spearman_deseq2_vs_paired_min_rows_with_findings', value=float(hit.effect_spearman_deseq2_vs_paired_logcpm.min())))
        rows.append(dict(version=label, quantity='spearman_deseq2_vs_paired_max_rows_with_findings', value=float(hit.effect_spearman_deseq2_vs_paired_logcpm.max())))
        rows.append(dict(version=label, quantity='spearman_deseq2_vs_paired_min_all_rows', value=float(rb.effect_spearman_deseq2_vs_paired_logcpm.min())))
        rows.append(dict(version=label, quantity='spearman_deseq2_vs_paired_max_all_rows', value=float(rb.effect_spearman_deseq2_vs_paired_logcpm.max())))
        rows.append(dict(version=label, quantity='n_cell_type_contrasts_analyzed', value=int(len(s))))
    df = pd.DataFrame(rows)
    w = df.pivot(index='quantity', columns='version', values='value').reset_index()
    w['changed'] = [('yes' if not np.isclose(float(a), float(b)) else 'no') for a, b in zip(w.committed, w.regenerated)]
    w.to_csv(OUT / 'R3_pseudobulk_committed_vs_regenerated.csv', index=False)
    s0 = pd.read_csv(core.CORR_OUT / 'pseudobulk_task4_summary.csv'); s1 = pd.read_csv(PB_NEW / 'pseudobulk_task4_summary.csv')
    m = s0.merge(s1, on=['cell_type', 'contrast'], suffixes=('_committed', '_regenerated'), how='outer')
    m[['cell_type', 'contrast', 'n_genes_tested_committed', 'n_genes_tested_regenerated', 'n_wald_padj_lt_0p05_committed',
       'n_wald_padj_lt_0p05_regenerated', 'n_welch_q_lt_0p05_committed', 'n_welch_q_lt_0p05_regenerated']].to_csv(
        OUT / 'R3_pseudobulk_per_celltype_contrast.csv', index=False)
    print(w.to_string(index=False))


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--part', required=True, choices=['comparison', 'consensus', 'enrichment', 'pseudobulk'])
    a = ap.parse_args([x for x in sys.argv[1:] if not x.startswith('--threads=')])
    globals()[a.part]()
