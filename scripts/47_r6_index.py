"""R6: INDEX.md listing every result file of both runs with a one-line description and the producing script.

Files without an entry in DESC are listed under 'Undescribed files' so that no file is silently omitted.
"""
import sys
from pathlib import Path

R = Path(__file__).resolve().parents[1]
D3 = R / 'analysis' / 'peerreview_20261003'; D4 = R / 'analysis' / 'peerreview_20261004'

S23, S24, S25, S26, S27, S28 = ('scripts/23_prepare_and_ved_check.py', 'scripts/24_canonical_fits.py', 'scripts/25_shuffled_label_attribution.py',
                                'scripts/26_c1_marker_list.py', 'scripts/27_c1_marker_comparison.py', 'scripts/28_c5_control_identity.py')
S29, S30, S31, S32, S33, S34, S35 = ('scripts/29_baselines_extended.py', 'scripts/30_calibration_pseudobulk.py', 'scripts/31_semisynthetic.py',
                                    'scripts/32_c6_ttr_distribution.py', 'scripts/33_d1_report_tables.py', 'scripts/34_c4_multigroupvi.py',
                                    'scripts/35_d1_signal_check.py')
S37, S38, S39, S40, S44, S45, S46, S47 = ('scripts/37_r2_baseline_headline.py', 'scripts/38_r1_injected_signal.py', 'scripts/39_r3_rerun_canonical_scripts.py',
                                         'scripts/40_r3_analyses.py', 'scripts/44_r4_mccvi_and_collect.py', 'scripts/45_r3_paper_numbers.py',
                                         'scripts/46_r5_archive.py', 'scripts/47_r6_index.py')
MAN = 'written by hand from result files'

DESC = {
    # ---------------- run 1 (peerreview_20261003)
    'A0_cells_comparison.csv': ('ved_archive cells.csv vs regenerated cell metadata, per column', S23),
    'A0_full0_latents_comparison.csv': ('max abs difference of latents recomputed from ved full_0.pt vs ved full_0_latents.npz', S23),
    'A0_full0_metrics_comparison.csv': ('metrics of ved full_0.pt vs committed full_0_metrics.json', S23),
    'A0_inputs_comparison.csv': ('ved inputs.npz vs regenerated canonical inputs, per array', S23),
    'A0_latent_probes_seed0_comparison.csv': ('seed-0 latent_probes rows recomputed from ved full_0 vs committed', S23),
    'A0_ved_archive_check.md': ('note: content check of the ved_archive files', MAN),
    'A1_consensus_jaccard.csv': ('fresh vs frozen consensus lists: sizes, overlap, Jaccard, gene lists', S24),
    'A1_cross_seed_cka.csv': ('cross-seed linear CKA per block, committed vs fresh', S24),
    'A1_fresh_canonical_fits.md': ('note: fresh canonical fits vs committed values', MAN),
    'A1_fresh_consensus.json': ('fresh consensus gene lists per block (ordered)', S24),
    'A1_fresh_consensus_genes.csv': ('fresh consensus genes with rank, top-50 count and frozen-list membership', S24),
    'A1_fresh_metrics.csv': ('per-seed metrics of the fresh canonical fits (MSE, F_sh, active units, KL, probes)', S24),
    'A1_fresh_ranking_top50_jaccard.csv': ('pairwise top-50 Jaccard among the six fresh seed x baseline rankings', S24),
    'A1_latent_probes_fresh.csv': ('latent_probes.csv rows (11_diagnostics logic) for the fresh fits', S24),
    'A1_latent_usage_per_dimension.csv': ('per-dimension KL and posterior-mean variance, fresh fits', S24),
    'A1_metric_comparison.csv': ('fresh vs committed metric per seed', S24),
    'A1_seed0_vs_ved.csv': ('fresh seed 0 vs ved full_0: metrics and test-cell CKA', S24),
    'C1_consensus_gene_annotation.csv': ('marker-set membership of each consensus gene (frozen and fresh lists)', S27),
    'C1_marker_comparison.md': ('note: consensus lists vs published marker lists', MAN),
    'C1_marker_list.csv': ('frozen marker list (cell_type, gene, source) from published sources', S26),
    'C1_marker_overlap.csv': ('overlap, hypergeometric P and BH q of consensus lists vs marker sets', S27),
    'C1_marker_set_sizes.csv': ('marker set sizes listed and in the 1,500-gene universe', S27),
    'C1_marker_sources.json': ('marker source files: SHA-256, sizes, symbol-matching rule, sources not used', S26),
    'C2_enrichment_all.csv.gz': ('all enrichment tests (shuffled, frozen, fresh lists)', S25),
    'C2_enrichment_summary.csv': ('number of FDR<0.05 terms per list x block x library', S25),
    'C2_label_shuffle_design.csv': ('label-shuffle design counts per split x study', S25),
    'C2_overlap.csv': ('shuffled vs real consensus overlaps (counts, Jaccard, hypergeometric P)', S25),
    'C2_paper_top_m8_term_lookup.csv': ('ZHANG_UTERUS_C5_MACROPHAGE results in every list', S25),
    'C2_shuffled_consensus.json': ('shuffled-label consensus lists', S25),
    'C2_shuffled_consensus_genes.csv': ('shuffled-label consensus genes with ranks', S25),
    'C2_shuffled_fit_metrics.csv': ('per-seed metrics of the shuffled-label fits', S25),
    'C2_shuffled_label_attribution.md': ('note: shuffled-label attribution', MAN),
    'C2_top_enrichment_terms.csv': ('top 5 terms per list x block x library', S25),
    'C3_baselines.md': ('note: scVI and contrastiveVI baselines', MAN),
    'C3_contrastivevi_metrics.csv': ('contrastiveVI per pair x seed: MSE both scales, probes, salient utilization', S29),
    'C3_contrastivevi_salient_usage_per_dimension.csv': ('contrastiveVI salient per-dimension KL/variance', S29),
    'C3_cross_seed_cka.csv': ('cross-seed CKA for scVI, contrastiveVI and MC-ContrastiveVI', S29),
    'C3_design.json': ('C3 settings and conversions', S29),
    'C3_mccvi_pca_both_scales.csv': ('MC-ContrastiveVI and PCA MSE on both scales', S29),
    'C3_mse_both_scales.csv': ('all C3 MSE rows on both scales', S29),
    'C3_mse_summary.csv': ('C3 MSE mean/SD per method x test-cell set', S29),
    'C3_scvi_metrics.csv': ('scVI per seed: MSE both scales, pooled and within-study probes', S29),
    'C4_mgvi_pip_freeze.txt': ('pip freeze of .venv-mgvi', 'pip freeze'),
    'C4_multigroupvi.md': ('note: multiGroupVI baseline', MAN),
    'C4_multigroupvi_cross_seed_cka.csv': ('multiGroupVI cross-seed CKA', S34),
    'C4_multigroupvi_metrics.csv': ('multiGroupVI per seed: MSE, probes, utilization', S34),
    'C4_multigroupvi_usage_per_dimension.csv': ('multiGroupVI per-dimension KL/variance', S34),
    'C5_block_norms_by_celltype.csv': ('ungated block norms per cell type, arm, scope', S28),
    'C5_block_norms_control_vs_cisplatin.csv': ('ungated block norms control vs cisplatin (wide)', S28),
    'C5_cell_counts.csv': ('cell counts per arm x split x cell type (GSE216146)', S28),
    'C5_celltype_probe_accuracy.csv': ('source cell-type probe balanced accuracy per block x arm x seed', S28),
    'C5_control_cell_identity.md': ('note: control-cell identity check', MAN),
    'C6_ttr_by_celltype_arm.csv': ('Ttr statistics per source cell type x arm', S32),
    'C6_ttr_choroid_plexus_vs_other.csv': ('Ttr statistics, Choroid Plexus vs other cell types', S32),
    'C6_ttr_distribution.md': ('note: Ttr distribution', MAN),
    'D1_base_cell_counts.csv': ('semi-synthetic base cells per library x pseudo-label x split', S31),
    'D1_calibration.json': ('effect-size calibration from the Task 4 pseudobulk analysis', S30),
    'D1_calibration_supported_rows.csv': ('direction-consistent BH-significant cisplatin rows used for calibration', S30),
    'D1_config_summary.csv': ('per-configuration summary (F_sh, true shared fraction, MSE, AU, KL, CKA)', S31),
    'D1_cross_seed_cka.csv': ('cross-seed CKA per configuration and block', S31),
    'D1_design.json': ('semi-synthetic design record', S31),
    'D1_fit_metrics.csv': ('per-fit metrics of the 75 grid fits', S31),
    'D1_gene_sets.csv': ('injected gene sets S, A, B with deciles and directions', S31),
    'D1_ig_consensus.csv': ('4-of-6 consensus recovery per configuration', S31),
    'D1_ig_per_config_mean.csv': ('per-fit IG recovery averaged per configuration x baseline', S31),
    'D1_ig_per_fit.csv': ('per-fit IG top-50 recovery', S31),
    'D1_injection_check.csv': ('realized log2 fold changes of the injection (data level)', 'runs/peerreview_20261003/scratch/check_injection.py (not committed)'),
    'D1_permutation_summary.csv': ('label-exchange reference summary', S31),
    'D1_permutation_values_null.csv': ('permuted F_sh values, null configuration', S31),
    'D1_permutation_values_shared_d1_rf1.0.csv': ('permuted F_sh values, shared_d1_rf1.0', S31),
    'D1_responder_counts.csv': ('responder counts per fraction x study', S31),
    'D1_semisynthetic.md': ('note: semi-synthetic positive control', MAN + ' and ' + S33),
    'D1_signal_check_seed0.csv': ('injected-set ranks by input mean difference vs IG ranks (seed 0)', S35),
    'PROGRESS.md': ('step log with times, commands and deviations', MAN),
    'SUMMARY.md': ('run summary', MAN),
    'environment.json': ('Python/platform/CPU and pip freeze of .venv', 'inline Python in run 1'),
    # ---------------- run A (peerreview_20261004)
    'R1_injected_signal_location.md': ('note: location of the injected semi-synthetic signal', MAN),
    'R1_probes_per_fit.csv': ('within-study probes per grid fit x study x representation x label', S38),
    'R1_reconstruction_per_fit.csv': ('observed vs reconstructed injected effects per fit x set x study x direction', S38),
    'R1_summary.csv': ('R1 probe and reconstruction summaries stacked', S38),
    'R1_summary_probes.csv': ('probe means/SDs per configuration', S38),
    'R1_summary_reconstruction.csv': ('reconstruction means/SDs per configuration', S38),
    'R2_baseline_headline.csv': ('one row per model: MSE, ungated probes, utilization, CKA, gated multiGroupVI probe', S37),
    'R2_baseline_headline.md': ('note: baseline headline table', MAN + ' and ' + S37),
    'R3_comparison_models.csv': ('per-fit MSE, background probes and F_sh for full, ablations, NB, PCA', S40),
    'R3_comparison_models_summary.csv': ('ranges per comparison model', S40),
    'R3_consensus_top12_figure5_data.csv': ('Figure-5 data: top-12 fresh consensus genes, per-seed mean |IG|', S40),
    'R3_consistent_numbers.md': ('note: consistent numbers from the fresh fits', MAN),
    'R3_enrichment_all.csv.gz': ('all enrichment tests, fresh and frozen consensus lists', S40),
    'R3_enrichment_cross_checks.csv': ('consensus overlap with regenerated cisplatin pseudobulk set and scDisInFact top-100', S40),
    'R3_enrichment_named_terms.csv': ('manuscript-named terms in every list', S40),
    'R3_enrichment_summary.csv': ('FDR<0.05 term counts per list x block x library', S40),
    'R3_enrichment_top_terms.csv': ('top 5 terms per list x block x library', S40),
    'R3_paper_numbers.csv': ('every manuscript number with fresh/committed value and source file', S45),
    'R3_pseudobulk_committed_vs_regenerated.csv': ('manuscript pseudobulk quantities, committed vs regenerated', S40),
    'R3_pseudobulk_per_celltype_contrast.csv': ('pseudobulk counts per cell type x contrast, committed vs regenerated', S40),
    'R3_09_interpretation/attribution_stability.csv': ('09 logic on fresh fits: seed/baseline/bootstrap top-50 Jaccard and Spearman', S39 + ' --script 09'),
    'R3_09_interpretation/complete_corrected_attributions.csv': ('09 logic: per-gene mean and mean |IG| per seed x axis x baseline x target', S39 + ' --script 09'),
    'R3_09_interpretation/ig_completeness.csv': ('09 logic: completeness errors and quadrature orders', S39 + ' --script 09'),
    'R3_09_interpretation/pathway_enrichment.csv': ('09 logic: enrichment of the top-50 mean-|IG| genes', S39 + ' --script 09'),
    'R3_11_diagnostics/ablation_diagnostics.csv': ('11 logic: background cell-type accuracy and block CKA per model x seed', S39 + ' --script 11'),
    'R3_11_diagnostics/cross_seed_latent_stability.csv': ('11 logic: cross-seed CKA per block', S39 + ' --script 11'),
    'R3_11_diagnostics/latent_dependence.csv': ('11 logic: between-block CKA with 100-permutation null', S39 + ' --script 11'),
    'R3_11_diagnostics/latent_probes.csv': ('11 logic: cell-type and within-study treatment probes', S39 + ' --script 11'),
    'R3_11_diagnostics/rescue_exploratory_comparisons.csv': ('11 logic: rescue-minus-treatment comparisons and sign-flip P', S39 + ' --script 11'),
    'R3_11_diagnostics/rescue_library_projection.csv': ('11 logic: per-library distances to the training-control centroid per seed', S39 + ' --script 11'),
    'R3_11_diagnostics/rescue_replicate_pairs.csv': ('11 logic: seed-averaged per-pair rescue differences', S39 + ' --script 11'),
    'R3_11_diagnostics/rescue_seedwise_signflip.csv': ('11 logic: per-seed sign-flip results', S39 + ' --script 11'),
    'R3_11_diagnostics/rescue_signflip_distribution.csv': ('11 logic: the 8 sign assignments', S39 + ' --script 11'),
    'R3_11_diagnostics/variance_robustness.csv': ('11 logic: F_sh with 200-resample cell-bootstrap intervals', S39 + ' --script 11'),
    'R4_baseline_positive_control.md': ('note: baseline positive control', MAN),
    'R4_comparison_table.csv': ('R4 comparison table: model x configuration x group', S44 + ' --collect'),
    'R4_contrastivevi_fits.csv': ('R4 contrastiveVI per-fit metrics', S44 + ' --collect'),
    'R4_ig_recovery_per_fit.csv': ('R4 attribution recovery per fit x group x baseline (all models)', S44 + ' --collect'),
    'R4_ig_recovery_summary.csv': ('R4 attribution recovery summaries', S44 + ' --collect'),
    'R4_metrics_per_fit.csv': ('R4 MSE, probes, active units, KL per model x fit x group', S44 + ' --collect'),
    'R4_metrics_summary.csv': ('R4 metric summaries', S44 + ' --collect'),
    'R4_multigroupvi_fits.csv': ('R4 multiGroupVI per-fit metrics', S44 + ' --collect'),
    'MANIFEST.csv': ('copy of the archive manifest (file, bytes, SHA-256, producer, command, commit)', S46),
    'ARCHIVE_ZIPS.csv': ('archive zip files with bytes and SHA-256', S46),
    'EXCLUDED_INPUT_HASHES.csv': ('hashes of inputs not archived (raw GEO, recovered cohort, MSigDB)', S46),
    'REPRODUCE.md': ('ordered reproduction commands for both runs', MAN),
    'INDEX.md': ('this index', S47),
    'environment_notes.md': ('environment notes for Run A', MAN),
}


def main():
    L = ['# Result index (both peer-review runs)', '', 'Generated by `scripts/47_r6_index.py`. One line per result file: description; producing script.', '']
    missing = []
    for title, D in (('analysis/peerreview_20261003 (previous run)', D3), ('analysis/peerreview_20261004 (Run A)', D4)):
        L += [f'## {title}', '', '| file | description | produced by |', '|---|---|---|']
        for p in sorted(x for x in D.rglob('*') if x.is_file()):
            rel = p.relative_to(D).as_posix()
            if rel in DESC:
                d, s = DESC[rel]; L.append(f'| `{rel}` | {d} | {s} |')
            else:
                missing.append(f'{title}: {rel}')
        L.append('')
    if missing:
        L += ['## Undescribed files', ''] + [f'- {m}' for m in missing]
    (D4 / 'INDEX.md').write_text('\n'.join(L) + '\n', encoding='utf-8')
    print('undescribed:', missing)


if __name__ == '__main__':
    sys.exit(main())
