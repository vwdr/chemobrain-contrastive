"""R3.9: table of every numeric result reported in manuscript_final_text.txt with fresh (or committed) values.

Fresh values are read from result files of this run and the previous run; values that do not depend on the
archived full-model fits are taken from the committed files and marked as such in `notes`. Methods settings that
are numeric are listed with section 'methods (setting)' and the value used by the fresh fits.
Output: analysis/peerreview_20261004/R3_paper_numbers.csv
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

R = Path(__file__).resolve().parents[1]
C = R / 'analysis' / 'corrected_20260920'; P = R / 'analysis' / 'peerreview_20261003'; O = R / 'analysis' / 'peerreview_20261004'
TXT = (R.parent / 'manuscript_final_text.txt').read_text(encoding='utf-8')
REL = lambda p: str(Path(p).relative_to(R)).replace('\\', '/')


def rng(v, fmt='{:.4f}'):
    v = np.asarray(v, dtype=float); return f'{fmt.format(v.min())}–{fmt.format(v.max())}'


def pct(v):
    v = np.asarray(v, dtype=float) * 100; return f'{v.min():.1f}–{v.max():.1f}%'


rows = []


def add(section, quote, quantity, manuscript_value, fresh_value, changed, source_file, notes=''):
    q = ' '.join(quote.split())
    found = q in ' '.join(TXT.split())
    rows.append(dict(section=section, quote=quote, quantity=quantity, manuscript_value=manuscript_value, fresh_value=fresh_value,
                     changed=changed, source_file=source_file, notes=(notes + ('' if found else ' [quote not found verbatim in text extraction; check line breaks]')).strip()))


def main():
    # ------------ sources
    qc = pd.read_csv(C / 'cohort_qc.csv'); la = pd.read_csv(C / 'label_audit.csv')
    a1 = pd.read_csv(P / 'A1_fresh_metrics.csv'); comp = pd.read_csv(O / 'R3_comparison_models.csv', keep_default_na=False, na_values=[''])
    vr = pd.read_csv(O / 'R3_11_diagnostics' / 'variance_robustness.csv'); ab = pd.read_csv(O / 'R3_11_diagnostics' / 'ablation_diagnostics.csv')
    lp = pd.read_csv(O / 'R3_11_diagnostics' / 'latent_probes.csv'); t1 = pd.read_csv(C / 'task1_reproduction_combined.csv')
    t2 = pd.read_csv(C / 'task2_null_test_summary.csv').iloc[0]; t5 = pd.read_csv(C / 'task5_unseen_library_model_summary.csv').set_index('mode')
    sd = pd.read_csv(C / 'task6_scdisinfact_summary.csv').set_index('metric')
    acv = pd.read_csv(C / 'annotation_reference_cv.csv'); acn = pd.read_csv(C / 'annotation_counts.csv')
    st = pd.read_csv(O / 'R3_09_interpretation' / 'attribution_stability.csv'); st = st[st.comparison == 'seed']
    cons = json.loads((P / 'A1_fresh_consensus.json').read_text())
    es = pd.read_csv(O / 'R3_enrichment_summary.csv').set_index(['consensus_version', 'axis', 'library'])
    en = pd.read_csv(O / 'R3_enrichment_all.csv.gz'); cx = pd.read_csv(O / 'R3_enrichment_cross_checks.csv')
    pb = pd.read_csv(O / 'R3_pseudobulk_committed_vs_regenerated.csv').set_index('quantity')
    rx = pd.read_csv(O / 'R3_11_diagnostics' / 'rescue_exploratory_comparisons.csv').iloc[0]
    sig = pd.read_csv(O / 'R3_11_diagnostics' / 'rescue_signflip_distribution.csv')
    full = comp[comp.model == 'full']; nb = comp[comp.model == 'negative_binomial']; pca = comp[comp.model == 'pca32'].iloc[0]
    ccq = REL(C / 'cohort_qc.csv')
    COMMITTED = 'Committed value; analysis does not depend on the archived full-model fits.'

    # ------------ Abstract
    add('abstract', 'we analyzed 28,808 observations from the noninterven- tion groups', 'nonintervention observations', '28,808',
        f"{int(la[~la.rescue].n_cells.sum()):,}", 'no', REL(C / 'label_audit.csv'), COMMITTED + ' Recomputed from the recovered cohort (identical SHA-256).')
    add('abstract', 'The mean squared error (MSE) for MC- ContrastiveVI ranged from 0.1687 to 0.1688', 'MC-ContrastiveVI held-out MSE range',
        '0.1687–0.1688', rng(a1.test_mse), 'yes', REL(P / 'A1_fresh_metrics.csv'), 'Fresh seeds 0–2: ' + '/'.join(f'{x:.6f}' for x in a1.test_mse))
    add('abstract', 'while PCA reached 0.1552', 'PCA(32) held-out MSE', '0.1552', f'{pca.test_mse:.4f}', 'no', REL(O / 'R3_comparison_models.csv'), f'{pca.test_mse:.8f}')
    add('abstract', 'the shared treatment block accounted for 92.1% to 95.7% of the treatment-associated posterior-mean vari- ance',
        'pooled F_sh range (three primary fits)', '92.1–95.7%', pct(vr[vr.drug == 'pooled'].fraction), 'yes', REL(O / 'R3_11_diagnostics' / 'variance_robustness.csv'))
    add('abstract', 'the shared fraction across the nine additional fits ranged from 75.0% to 98.1%', 'F_sh range, nine Task 1 reruns',
        '75.0–98.1%', pct(t1.new_shared_fraction), 'no', REL(C / 'task1_reproduction_combined.csv'), COMMITTED + ' (reruns of the archived configuration on three hosted CPU environments; not rerun here).')
    add('abstract', 'label exchangeability across 16 replicate-preserving assignments', 'number of label assignments', '16', str(int(t2.n_assignments)), 'no', REL(C / 'task2_null_test_summary.csv'), COMMITTED)
    add('abstract', 'Twelve of the 16 partitions produced a shared fraction at least as large as the observed one (P = 0.750)', 'label-exchange upper-tail count and P',
        '12 of 16; P = 0.750', f'{int(t2.upper_tail_count)} of 16; P = {t2.exact_upper_tail_p:.3f}', 'no', REL(C / 'task2_null_test_summary.csv'), COMMITTED)

    # ------------ 2.1 cohort
    g = qc.groupby('study')[['before_qc', 'after_qc']].sum(); libs = qc.groupby('study').sample_id.nunique()
    add('2.1', 'GSE216146 contained 19,834 deposited whole-cell profiles from 12 libraries', 'GSE216146 deposited profiles / libraries', '19,834 / 12',
        f"{int(g.loc['GSE216146','before_qc']):,} / {libs['GSE216146']}", 'no', ccq, COMMITTED)
    add('2.1', 'GSE271055 contained 35,322 de- posited nuclei from three libraries', 'GSE271055 deposited nuclei / libraries', '35,322 / 3',
        f"{int(g.loc['GSE271055','before_qc']):,} / {libs['GSE271055']}", 'no', ccq, COMMITTED)
    add('2.1', 'Quality filtering re- tained 14,965 cells and 31,892 nuclei across 18,271 com- mon genes', 'retained cells / nuclei / genes', '14,965 / 31,892 / 18,271',
        f"{int(g.loc['GSE216146','after_qc']):,} / {int(g.loc['GSE271055','after_qc']):,} / 18,271", 'no', ccq, COMMITTED + ' Gene count from recovered_counts.h5ad (46,857 × 18,271).')
    add('2.1', 'Each GSE271055 library pooled hippocampal material from four mice', 'mice pooled per GSE271055 library', '4', '4', 'not-recomputable', 'data/evidence/GSE271055.soft', 'Source-study metadata; not computed.')
    add('2.1', 'The resulting analysis cohort contained 28,808 noninterven- tion observations', 'nonintervention observations', '28,808', f"{int(la[~la.rescue].n_cells.sum()):,}", 'no', REL(C / 'label_audit.csv'), COMMITTED)
    add('2.1', 'The remaining 18,049 intervention-arm observations were reserved for projection', 'intervention-arm observations', '18,049', f"{int(la[la.rescue].n_cells.sum()):,}", 'no', REL(C / 'label_audit.csv'), COMMITTED)
    arms = [('GSE216146', 'control', 'PBS', 'GSE216146 PBS 3 5,989 4,830 Whole cells'), ('GSE216146', 'control_GENUS', 'PBS + GENUS', 'PBS + GENUS 3 4,868 4,098 Whole cells'),
            ('GSE216146', 'cisplatin', 'Cisplatin', 'Cisplatin 3 5,749 4,005 Whole cells'), ('GSE216146', 'cisplatin_GENUS', 'Cisplatin + GENUS', 'Cisplatin + GENUS 3 3,228 2,032 Whole cells'),
            ('GSE271055', 'control', 'Control', 'GSE271055 Control 1 10,909 9,455 Nuclei'), ('GSE271055', 'doxorubicin', 'Doxorubicin', 'Doxorubicin 1 12,221 10,518 Nuclei'),
            ('GSE271055', 'doxorubicin_ACY1083', 'Doxorubicin + ACY-1083', 'Doxorubicin + ACY-1083 1 12,192 11,919 Nuclei')]
    tab = qc.groupby(['study', 'arm']).agg(libs=('sample_id', 'nunique'), dep=('before_qc', 'sum'), ret=('after_qc', 'sum'))
    for study, arm, lab, quote in arms:
        r = tab.loc[(study, arm)]
        add('Table 1 / Figure 2A', quote, f'{lab}: libraries / deposited / retained', 'see Table 1', f'{r.libs} / {r.dep:,} / {r.ret:,}', 'no', ccq,
            COMMITTED + ' Table 1 values (manuscript: PBS 3/5,989/4,830; PBS+GENUS 3/4,868/4,098; Cisplatin 3/5,749/4,005; Cisplatin+GENUS 3/3,228/2,032; Control 1/10,909/9,455; Doxorubicin 1/12,221/10,518; Doxorubicin+ACY-1083 1/12,192/11,919).')

    add('Figure 1A', 'GSE216146 12 whole-cell libraries 14,965 retained cells', 'GSE216146 libraries / retained cells', '12 / 14,965',
        f"{libs['GSE216146']} / {int(g.loc['GSE216146','after_qc']):,}", 'no', ccq, COMMITTED)
    add('Figure 1A', 'GSE271055 3 pooled-nucleus libraries 31,892 retained nuclei', 'GSE271055 libraries / retained nuclei', '3 / 31,892',
        f"{libs['GSE271055']} / {int(g.loc['GSE271055','after_qc']):,}", 'no', ccq, COMMITTED)
    # ------------ 2.2 reconstruction
    add('2.2', 'The primary benchmark used 1,500 genes selected from training cells and 32 total latent coordinates', 'genes / latent coordinates', '1,500 / 32', '1,500 / 32', 'no', REL(C / 'benchmark_gene_universe.csv'), 'Design setting; committed gene universe used by the fresh fits.')
    add('2.2', 'Three inde- pendently initialized MC-ContrastiveVI fits achieved held- out log-expression mean squared error of 0.1687–0.1688', 'MC-ContrastiveVI held-out MSE range',
        '0.1687–0.1688', rng(a1.test_mse), 'yes', REL(P / 'A1_fresh_metrics.csv'))
    add('2.2', 'PCA with 32 components achieved lower error at 0.1552', 'PCA(32) MSE', '0.1552', f'{pca.test_mse:.4f}', 'no', REL(O / 'R3_comparison_models.csv'))
    for mode, lab in (('no_hsic', 'no HSIC'), ('no_gating', 'no gating'), ('gaussian_vae', 'noncontrastive Gaussian VAE')):
        cm = pd.read_csv(C / 'model_benchmarks.csv'); cm = cm[cm['mode'] == mode].test_mse
        fm = comp[comp.model == mode].test_mse
        add('2.2', 'Removing the independence penalty or hard gating did not improve reconstruction materially, and a noncontrastive Gaussian variational autoencoder produced similar error',
            f'{lab} held-out MSE range (qualitative statement; values from Figure 3A data)', f'committed {rng(cm)}', rng(fm), 'yes', REL(O / 'R3_comparison_models.csv'),
            'Manuscript gives no number in text; committed range from analysis/corrected_20260920/model_benchmarks.csv.')
    add('2.2', 'The negative-binomial sensitivity model produced higher error of 0.2109–0.2123', 'NB sensitivity held-out MSE range', '0.2109–0.2123', rng(nb.test_mse), 'yes', REL(O / 'R3_comparison_models.csv'))
    add('2.2', 'The full model had mean balanced MSE of 0.0910 on vali- dation cells from represented libraries and 0.0940 on the excluded libraries',
        'Task 5 full model seen / unseen balanced MSE', '0.0910 / 0.0940', f"{t5.loc['full','mean_seen_library_validation_balanced_mse']:.4f} / {t5.loc['full','mean_unseen_library_balanced_mse']:.4f}", 'no', REL(C / 'task5_unseen_library_model_summary.csv'), COMMITTED)
    add('2.2', 'The noncontrastive Gaussian model produced 0.0926 on excluded libraries, and PCA produced 0.0860', 'Task 5 Gaussian VAE / PCA unseen balanced MSE', '0.0926 / 0.0860',
        f"{t5.loc['gaussian_vae','mean_unseen_library_balanced_mse']:.4f} / {t5.loc['pca32','mean_unseen_library_balanced_mse']:.4f}", 'no', REL(C / 'task5_unseen_library_model_summary.csv'), COMMITTED)
    add('2.2', 'scDisInFact produced mean held-out log-expression MSE of 0.2313', 'scDisInFact MSE', '0.2313', f"{sd.loc['test_log_expression_mse','mean']:.4f}", 'no', REL(C / 'task6_scdisinfact_summary.csv'), COMMITTED)
    add('2.2', 'Balanced study accuracy was 0.895 for its shared-bio factor and 0.981 for its condition-associated factor', 'scDisInFact study accuracy shared-bio / condition factor',
        '0.895 / 0.981', f"{sd.loc['shared_bio_study_balanced_accuracy','mean']:.3f} / {sd.loc['unshared_bio_study_balanced_accuracy','mean']:.3f}", 'no', REL(C / 'task6_scdisinfact_summary.csv'), COMMITTED)
    add('2.2', 'MC-ContrastiveVI background representations like- wise retained study information, with balanced accuracy of 0.726–0.856',
        'background study balanced accuracy range', '0.726–0.856', rng(full.bg_study_balanced_accuracy, '{:.3f}'), 'yes', REL(O / 'R3_comparison_models.csv'), 'From fresh metrics JSON (06 probes).')
    add('2.2', 'preserving GSE216146 source cell-type labels with balanced accuracy of 0.985–0.991', 'background source cell-type balanced accuracy range',
        '0.985–0.991', rng(lp[(lp.representation == 'background') & (lp.task == 'cis_source_cell_type')].balanced_accuracy, '{:.3f}'), 'yes', REL(O / 'R3_11_diagnostics' / 'latent_probes.csv'))
    add('Figure 3B caption', 'Chance is 0.5', 'chance level of study probe', '0.5', '0.5', 'no', '', 'Definition (two studies, balanced accuracy).')

    # ------------ 2.3 latent utilization
    add('2.3', 'the shared block accounted for 92.1–95.7% of treatment-associated posterior-mean vari- ance', 'pooled F_sh range', '92.1–95.7%',
        pct(vr[vr.drug == 'pooled'].fraction), 'yes', REL(O / 'R3_11_diagnostics' / 'variance_robustness.csv'),
        'Per drug: doxorubicin ' + pct(vr[vr.drug == 'doxorubicin'].fraction) + '; cisplatin ' + pct(vr[vr.drug == 'cisplatin'].fraction) + '; bootstrap intervals in the source file.')
    add('2.3', 'the three fits contained 0, 0 and 1 active shared coordinate out of 8', 'active shared coordinates per fit', '0, 0, 1',
        ', '.join(str(int(x)) for x in a1.shared_active_units_var_gt_0p01), 'yes', REL(P / 'A1_fresh_metrics.csv'))
    add('2.3', 'All four doxorubicin-specific and all four cisplatin-specific coordinates were inactive in every fit', 'active drug coordinates per fit', '0 (dox), 0 (cis)',
        ', '.join(f'{int(x)}/{int(y)}' for x, y in zip(a1.dox_active_units_var_gt_0p01, a1.cis_active_units_var_gt_0p01)), 'no', REL(P / 'A1_fresh_metrics.csv'))
    add('2.3', 'no individual treatment-associated coordinate had mean KL above 0.01 nats', 'coordinates with mean KL > 0.01', '0',
        str(int((a1.shared_dims_kl_gt_0p01 + a1.dox_dims_kl_gt_0p01 + a1.cis_dims_kl_gt_0p01).sum())), 'no', REL(P / 'A1_fresh_metrics.csv'),
        'Max per-dimension mean KL across fresh fits/blocks: ' + f"{a1[['shared_max_dim_kl_nats','dox_max_dim_kl_nats','cis_max_dim_kl_nats']].values.max():.4f} nats.")
    add('2.3', 'Nine additional fits performed with identical data, hyperpa- rameters and software settings across three hosted CPU environments produced pooled shared fractions from 75.0% to 98.1%',
        'F_sh range, nine reruns', '75.0–98.1%', pct(t1.new_shared_fraction), 'no', REL(C / 'task1_reproduction_combined.csv'), COMMITTED)
    add('2.3', 'while held-out MSE remained within 0.168516– 0.168816', 'MSE range, nine reruns', '0.168516–0.168816', rng(t1.new_test_mse, '{:.6f}'), 'no', REL(C / 'task1_reproduction_combined.csv'), COMMITTED)
    add('2.3', 'We evaluated all 16 assign- ments', 'label assignments', '16', str(int(t2.n_assignments)), 'no', REL(C / 'task2_null_test_summary.csv'), COMMITTED)
    add('2.3', 'within the three GSE216146 replicate pairs', 'GSE216146 replicate pairs', '3', '3', 'no', REL(C / 'task2_library_design.csv'), COMMITTED)
    add('2.3', 'The observed assignment had mean shared fraction 0.776', 'Task 2 observed mean F_sh', '0.776', f'{t2.observed_value:.3f}', 'no', REL(C / 'task2_null_test_summary.csv'), COMMITTED)
    add('2.3', 'Twelve of 16 assignments produced values at least as large, giving an upper-tail reference probability of P = 0.750', 'Task 2 upper tail', '12 of 16; P = 0.750',
        f'{int(t2.upper_tail_count)} of 16; P = {t2.exact_upper_tail_p:.3f}', 'no', REL(C / 'task2_null_test_summary.csv'), COMMITTED)
    add('2.3', 'The reference median was 0.888', 'Task 2 reference median', '0.888', f'{t2.reference_median:.3f}', 'no', REL(C / 'task2_null_test_summary.csv'), COMMITTED)
    add('2.3', 'The negative-binomial sensitivity model also assigned a high fraction to the shared block, 82.9–89.2%', 'NB F_sh range', '82.9–89.2%', pct(nb.shared_fraction), 'yes', REL(O / 'R3_comparison_models.csv'))

    # ------------ 2.4 annotation and attribution
    add('2.4', 'Leave-library-out balanced accuracy within the reference was 0.909–0.955', 'annotation leave-library-out balanced accuracy', '0.909–0.955',
        rng(acv.balanced_accuracy, '{:.3f}'), 'no', REL(C / 'annotation_reference_cv.csv'), COMMITTED + ' (08_annotation.py; independent of model fits).')
    d = acn[acn.study == 'GSE271055']
    add('2.4', '72.6% of GSE271055 nuclei failed the joint confidence and marker criterion', 'unassigned GSE271055 nuclei', '72.6%',
        f"{100 * d[d.annotation == 'uncertain'].n_cells.sum() / d.n_cells.sum():.1f}%", 'no', REL(C / 'annotation_counts.csv'), COMMITTED)
    s0 = st[(st.axis == 'shared') & (st.baseline == 'zero') & (st.target == 'norm_squared')].top50_jaccard
    add('2.4', 'Shared-block top-50 overlap across initializations ranged from 0.190 to 0.639 with the squared-norm target and zero baseline',
        'shared top-50 cross-seed Jaccard (norm squared, zero baseline)', '0.190–0.639', rng(s0, '{:.3f}'), 'yes', REL(O / 'R3_09_interpretation' / 'attribution_stability.csv'))
    add('2.4', 'consensus gene as one appearing in the top 50 in at least four of six initialization-by-baseline rankings', 'consensus rule', 'top 50, 4 of 6', 'top 50, 4 of 6', 'no', '', 'Rule; applied unchanged.')
    add('2.4', 'This yielded 34 shared, 20 doxorubicin-associated and 33 cisplatin-associated genes', 'consensus list sizes', '34 / 20 / 33',
        f"{len(cons['shared'])} / {len(cons['doxorubicin'])} / {len(cons['cisplatin'])}", 'yes', REL(P / 'A1_fresh_consensus.json'))
    fig5 = {'shared': 'Mrc1, Pf4, Ms4a7, Cbr2, Cd163, F13a1, Gpx3, Dab2, mt-Co3, Cd52, Stab1, C5ar1',
            'doxorubicin': 'Snhg11, Kcnq1ot1, Meg3, Ttr, Gria2, Rtn1, Grin2b, C1ql3, Opcml, Pcdh9, Nrcam, Dclk1',
            'cisplatin': 'Ccn3, Col25a1, Igfbp6, Mgp, Foxc2, Rspo3, Col6a2, mt-Nd4, Ogn, Col6a1, Adamtsl3, Foxd1'}
    for blk in ('shared', 'doxorubicin', 'cisplatin'):
        add('Figure 5', 'Figure 5 shows the twelve highest-ranked consensus genes in each block', f'top-12 consensus genes, {blk}', fig5[blk],
            ', '.join(cons[blk][:12]), 'yes', REL(O / 'R3_consensus_top12_figure5_data.csv'),
            f'Overlap of top-12 sets: {len(set(fig5[blk].split(", ")) & set(cons[blk][:12]))}/12. Per-seed mean |IG| (two-baseline average) in the source file.')
    for blk, ax in (('shared', 'shared'), ('doxorubicin', 'drug0'), ('cisplatin', 'drug1')):
        for tgt in ('coordinate_sum', 'norm_squared'):
            v = st[(st.axis == ax) & (st.target == tgt)].top50_jaccard
            add('Figure 4B', 'Pairwise top-50 Jaccard overlap for integrated-gradient rankings across model initializations', f'cross-seed top-50 Jaccard, {blk}, {tgt} (both baselines)',
                'figure only', rng(v, '{:.3f}') + f' (mean {v.mean():.3f})', 'not-recomputable', REL(O / 'R3_09_interpretation' / 'attribution_stability.csv'),
                'Figure values not stated in text; committed values in analysis/corrected_20260920/attribution_stability.csv.')
    for drug in ('pooled', 'doxorubicin', 'cisplatin'):
        v = vr[vr.drug == drug]
        add('Figure 4A', 'Intervals are conditional cell-bootstrap percentiles', f'F_sh with bootstrap 95% interval, {drug}', 'figure only',
            '; '.join(f's{int(r.seed)} {r.fraction:.3f} [{r.cell_bootstrap_low:.3f}, {r.cell_bootstrap_high:.3f}]' for r in v.itertuples()), 'not-recomputable',
            REL(O / 'R3_11_diagnostics' / 'variance_robustness.csv'), 'Figure values; committed values in analysis/corrected_20260920/variance_robustness.csv.')

    # ------------ 2.5 enrichment and pseudobulk
    def n(ver, axis, lib):
        return int(es.loc[(ver, axis, lib), 'n_fdr_005'])
    def term(ver, axis, t):
        r = en[(en.consensus_version == ver) & (en.axis == axis) & (en.term == t)].iloc[0]; return r
    add('2.5', 'M8 signatures and 11 Gene Ontology terms passed FDR 0.05', 'shared consensus: M8 / GO BP terms at FDR 0.05', '46 / 11',
        f"{n('fresh','shared','M8_cell_type')} / {n('fresh','shared','GO_BP')}", 'yes', REL(O / 'R3_enrichment_summary.csv'), f"'Forty-six' precedes this phrase across a page break in the text extraction. Frozen lists recomputed: {n('frozen','shared','M8_cell_type')} / {n('frozen','shared','GO_BP')}.")
    t = term('fresh', 'shared', 'ZHANG_UTERUS_C5_MACROPHAGE'); tf = term('frozen', 'shared', 'ZHANG_UTERUS_C5_MACROPHAGE')
    top_fresh = en[(en.consensus_version == 'fresh') & (en.axis == 'shared') & (en.library == 'M8_cell_type')].sort_values(['q_value', 'p_value']).iloc[0].term
    add('2.5', 'The strongest M8 result was a macrophage signature containing 14 of 34 shared consensus genes with adjusted P = 4.18 × 10−10',
        'top M8 term (ZHANG_UTERUS_C5_MACROPHAGE): overlap / q', '14 of 34; 4.18e-10', f'{int(t.overlap)} of {len(cons["shared"])}; {t.q_value:.3g}', 'yes', REL(O / 'R3_enrichment_all.csv.gz'),
        f'Top fresh M8 term: {top_fresh}. Frozen list recomputed: {int(tf.overlap)} of 34; q = {tf.q_value:.3g}.')
    t = term('fresh', 'shared', 'REACTOME_IMMUNE_SYSTEM'); tf = term('frozen', 'shared', 'REACTOME_IMMUNE_SYSTEM')
    add('2.5', 'Reactome immune system did not pass the same threshold after consensus filtering, with adjusted P = 0.092', 'Reactome immune system q (shared)', '0.092',
        f'{t.q_value:.3g} (overlap {int(t.overlap)})', 'yes', REL(O / 'R3_enrichment_all.csv.gz'), f'Frozen list recomputed: q = {tf.q_value:.3f}.')
    add('2.5', 'Seventy-seven M8 signatures and five Reactome terms passed FDR 0.05', 'cisplatin consensus: M8 / Reactome terms at FDR 0.05', '77 / 5',
        f"{n('fresh','cisplatin','M8_cell_type')} / {n('fresh','cisplatin','Reactome')}", 'yes', REL(O / 'R3_enrichment_summary.csv'),
        f"Frozen lists recomputed: {n('frozen','cisplatin','M8_cell_type')} / {n('frozen','cisplatin','Reactome')}.")
    t = term('fresh', 'cisplatin', 'REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX'); tf = term('frozen', 'cisplatin', 'REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX')
    add('2.5', 'Re- actome extracellular-matrix degradation contained five of 33 cisplatin consensus genes, corresponding to 7.58-fold enrichment and adjusted P = 0.034',
        'Reactome ECM degradation (cisplatin): overlap / fold / q', '5 of 33; 7.58; 0.034', f'{int(t.overlap)} of {len(cons["cisplatin"])}; {t.fold_enrichment:.2f}; {t.q_value:.3g}', 'yes',
        REL(O / 'R3_enrichment_all.csv.gz'), f'Frozen list recomputed: {int(tf.overlap)} of 33; {tf.fold_enrichment:.2f}; q = {tf.q_value:.3f}.')
    dx = sum(n('fresh', 'doxorubicin', l) for l in ('Reactome', 'GO_BP', 'M8_cell_type'))
    add('2.5', 'No Reactome, Gene Ontology or M8 term passed FDR 0.05 for the 20-gene dox- orubicin consensus', 'doxorubicin consensus: terms at FDR 0.05 (size)', '0 (20 genes)',
        f"{dx} ({len(cons['doxorubicin'])} genes): Reactome {n('fresh','doxorubicin','Reactome')}, GO {n('fresh','doxorubicin','GO_BP')}, M8 {n('fresh','doxorubicin','M8_cell_type')}", 'yes',
        REL(O / 'R3_enrichment_summary.csv'), 'Fresh doxorubicin Reactome terms at FDR 0.05: REACTOME_NEURONAL_SYSTEM, REACTOME_NEUREXINS_AND_NEUROLIGINS (see R3_enrichment_top_terms.csv). Frozen list recomputed: 0.')
    def pbv(qn, fmt='{:.0f}'):
        r = pb.loc[qn]; return fmt.format(float(r.committed)), fmt.format(float(r.regenerated)), r.changed
    for qn, quote, lab in (('n_BH_findings_cisplatin_vs_control', 'identified 297 within-contrast BH findings for cisplatin versus control', 'cisplatin vs control BH findings'),
                           ('n_BH_findings_cisplatin_GENUS_vs_cisplatin', '55 for cisplatin plus GENUS versus cisplatin', 'cisplatin+GENUS vs cisplatin BH findings'),
                           ('n_BH_findings_control_GENUS_vs_control', '2 for control plus GENUS versus control', 'control+GENUS vs control BH findings'),
                           ('pooled_BH_findings_all_contrasts', 'the pooled count of 354 findings', 'pooled BH findings'),
                           ('n_findings_all3_pairs_same_direction', 'Of these findings, 353 had the same effect direction in all three replicate pairs', 'findings with same direction in all 3 pairs'),
                           ('n_welch_fdr_005_total', 'Welch sensitivity analysis produced only 1 FDR finding', 'Welch FDR findings')):
        c_, r_, ch = pbv(qn)
        add('2.5', quote, lab, c_, r_, ch, REL(O / 'R3_pseudobulk_committed_vs_regenerated.csv'), 'Regenerated with scripts/16 logic (runs/peerreview_20261004/pseudobulk_task4); manuscript value = committed file value.')
    c_, r_, _ = pbv('spearman_deseq2_vs_paired_min_rows_with_findings', '{:.3f}'); c2, r2, _ = pbv('spearman_deseq2_vs_paired_max_rows_with_findings', '{:.3f}')
    add('2.5', 'Spearman ρ from 0.932 to 0.997', 'Spearman range DESeq2 vs paired effects', f'{c_}–{c2}', f'{r_}–{r2}', 'no' if (c_, c2) == (r_, r2) else 'yes',
        REL(O / 'R3_pseudobulk_committed_vs_regenerated.csv'), 'Range over cell-type contrasts with ≥ 1 BH finding (identical over all 16 contrasts).')
    add('2.5', 'The exact paired sign-flip reference has a minimum two-sided probability of 0.25', 'minimum two-sided sign-flip P', '0.25', '0.25', 'no', '', 'Design property (2^3 assignments).')
    cf = cx[(cx.consensus_version == 'frozen') & (cx.axis == 'cisplatin')].iloc[0]; cfr = cx[(cx.consensus_version == 'fresh') & (cx.axis == 'cisplatin')].iloc[0]
    add('2.5', 'Five cisplatin consensus genes, Apoe, Mbp, Mt1, Tmsb10 and Ttr, overlapped the direction-consistent cisplatin pseu- dobulk support set',
        'cisplatin consensus ∩ direction-consistent pseudobulk set', '5: Apoe, Mbp, Mt1, Tmsb10, Ttr', f"fresh list {int(cfr.n_overlap_cisplatin_DE)}: {cfr.cisplatin_DE_overlap_genes.replace('|', ', ')}",
        'yes', REL(O / 'R3_enrichment_cross_checks.csv'), f"Frozen list vs regenerated support set: {int(cf.n_overlap_cisplatin_DE)}: {cf.cisplatin_DE_overlap_genes.replace('|', ', ')}; regenerated support set = {int(cf.n_support_genes_in_universe)} universe genes.")

    # ------------ 2.6 rescue
    add('2.6', 'For cisplatin, the mean rescue-minus- treatment difference was 0.0003', 'cisplatin rescue-minus-treatment mean distance difference', '0.0003',
        f'{rx.mean_difference:.4f}', 'yes', REL(O / 'R3_11_diagnostics' / 'rescue_exploratory_comparisons.csv'), f'{rx.mean_difference:.6f}; per-pair values in rescue_replicate_pairs.csv.')
    add('2.6', 'The three GEO replicate pairs permit eight exact within-pair sign assignments', 'sign assignments', '8', str(len(sig)), 'no', REL(O / 'R3_11_diagnostics' / 'rescue_signflip_distribution.csv'))
    add('2.6', 'All eight produced an absolute mean difference at least as large as observed, giving P = 1.000', 'assignments ≥ observed; P', '8; 1.000',
        f'{int(sig.absolute_ge_observed.sum())}; {rx.exact_two_sided_p:.3f}', 'no', REL(O / 'R3_11_diagnostics' / 'rescue_exploratory_comparisons.csv'))
    add('2.6', 'The minimum attainable two-sided probability is 0.25', 'minimum attainable P', '0.25', f'{rx.minimum_attainable_two_sided_p:.2f}', 'no', REL(O / 'R3_11_diagnostics' / 'rescue_exploratory_comparisons.csv'))
    add('Figure 6', 'Each point is the mean distance for one deposited library, averaged across model initializations', 'per-library mean distances (figure data)', 'figure only',
        'see source file', 'not-recomputable', REL(O / 'R3_11_diagnostics' / 'rescue_library_projection.csv'), 'Figure data recomputed on the fresh fits (per seed); not stated numerically in text.')

    # ------------ methods: stated results and settings
    add('3.1', 'Validation and test sets each contained 2,881 observations', 'validation / test cells', '2,881 / 2,881', '2,881 / 2,881', 'no', 'runs/peerreview_20261003/canonical_inputs/inputs.npz', 'Regenerated split identical to committed (A0).')
    add('3.1', 'yielding 6,000 training observations', 'training cells', '6,000', '6,000', 'no', 'runs/peerreview_20261003/canonical_inputs/inputs.npz', '')
    add('3.1', 'We sampled 1,500 training observations from each of four study-by-treatment strata', 'training cells per stratum', '1,500', '1,500', 'no', 'runs/peerreview_20261003/canonical_inputs/inputs.npz', 'Setting.')
    for quote, qn, val in (('split with seed 1729 into an 80% candidate training pool and 10% each validation and test sets', 'split seed / proportions', '1729; 80/10/10'),
                           ('We selected 1,500 highly variable genes', 'HVGs', '1,500'),
                           ('Counts were normalized to 10,000 per observation', 'normalization target', '10,000'),
                           ('we retained observations with at least 500 detected genes, genes detected in at least ten observa- tions, and observations with mitochondrial count fraction below 10%', 'QC thresholds', '500 genes; 10 observations; 10% mito'),
                           ('Latent dimensions were 16 for background, 8 for the shared block and 4 for each treatment-specific block', 'latent dimensions', '16 / 8 / 4+4'),
                           ('Each encoder used two hidden layers with 128 units', 'encoder width/depth', '2 × 128'),
                           ('dropout 0.1', 'dropout', '0.1'),
                           ('We used β = 1 and (λ1, λ2, λ3) = (10, 10, 5)', 'KL weight; HSIC weights', '1; (10, 10, 5)'),
                           ('minibatches of approximately 256 observations', 'batch size', '256'),
                           ('gradient clipping at norm 5, and a maximum of 100 epochs', 'gradient clip; max epochs', '5; 100'),
                           ('with patience 15', 'early-stopping patience', '15'),
                           ('A coordinate was classified as active for this diagnostic when posterior-mean variance exceeded 0.01', 'active-unit threshold', '0.01'),
                           ('Mean KL above 0.01 nats was recorded separately', 'KL threshold', '0.01 nats'),
                           ('PCA used 32 components', 'PCA components', '32'),
                           ('The remaining observations supplied 1,500 training cells and 100 validation cells from each study-by-treatment stratum', 'Task 5 cells per stratum', '1,500 / 100'),
                           ('with latent dimensions 8 and 2, batch size 64, learning rate 5 × 10−4 , 50 epochs and three initializations', 'scDisInFact settings', '8,2; 64; 5e-4; 50; 3'),
                           ('Predictions were retained only when the maximum class probability was at least 0.8', 'annotation probability threshold', '0.8'),
                           ('Gene sets containing 10–500 modeled genes were tested', 'gene-set size limits', '10–500'),
                           ('Cell-type contrasts were analyzed only when all six required library aggregates contained at least 20 cells', 'pseudobulk min cells', '20'),
                           ('Genes required at least ten pseudobulk counts in at least three of six libraries', 'pseudobulk gene filter', '10 counts in 3 of 6'),
                           ('Py- DESeq2 0.5.4', 'PyDESeq2 version', '0.5.4'),
                           ('Nine additional training runs were performed across three independent CPU environments', 'Task 1 reruns', '9 runs; 3 environments')):
        add('methods (setting)', quote, qn, val, val, 'no', 'src/peerreview/core.py; scripts/06,09,11,16,19 (unmodified logic)', 'Setting used unchanged by the fresh analyses (or committed analysis, for Task 1/5/6/annotation).')
    add('methods (setting)', 'Models were trained with AdamW at learning rate 10−3 and weight decay 10−5', 'learning rate; weight decay', '1e-3; 1e-5', '1e-3; 1e-5', 'no', 'src/peerreview/core.py', 'Setting.')
    add('methods (setting)', 'we performed three leave-one-replicate-pair-out analyses within GSE216146', 'Task 5 folds', '3', '3', 'no', REL(C / 'task5_fold_design.csv'), COMMITTED)
    add('methods (setting)', 'fit with three initializations per fold', 'Task 5 initializations per fold', '3', '3', 'no', REL(C / 'task5_unseen_library_metrics.csv'), COMMITTED)
    add('discussion', 'Values above 90% initially sug- gested that most treatment-associated variation occupied the shared block', 'pooled F_sh above 90% in primary fits',
        '> 90%', pct(vr[vr.drug == 'pooled'].fraction), 'no', REL(O / 'R3_11_diagnostics' / 'variance_robustness.csv'), 'All fresh pooled values exceed 90%.')
    add('3.4', 'This produced 16 replicate-', 'label assignments', '16', '16', 'no', REL(C / 'task2_null_test_summary.csv'), COMMITTED + " Phrase continues 'preserving assignments' after a page number in the text extraction.")
    add('3.4', 'Each assignment was fit with three initializations', 'initializations per assignment', '3', '3', 'no', REL(C / 'task2_permutation_fits.csv'), COMMITTED)

    df = pd.DataFrame(rows); df.to_csv(O / 'R3_paper_numbers.csv', index=False)
    print(df.changed.value_counts()); print((df.notes.str.contains('quote not found')).sum(), 'quotes not found verbatim')
    print(df[df.notes.str.contains('quote not found')][['section', 'quote']].to_string())


if __name__ == '__main__':
    main()
