"""R2: baseline headline table assembled from existing result files only (no refits).

Sources (read-only): analysis/peerreview_20261003/{A1_*, C3_*, C4_*} and analysis/corrected_20260920/task6_*.
Within-study treatment probes are taken on ungated representations only. multiGroupVI probes on gated
(label-masked) group-specific latents are listed separately in `gated_label_masked_probe`.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

R = Path(__file__).resolve().parents[1]
P = R / 'analysis' / 'peerreview_20261003'; C = R / 'analysis' / 'corrected_20260920'; OUT = R / 'analysis' / 'peerreview_20261004'


def ms(v):
    v = np.asarray([x for x in v if x is not None and not (isinstance(x, float) and np.isnan(x))], dtype=float)
    if len(v) == 0:
        return np.nan, np.nan, ''
    return float(v.mean()), float(v.std(ddof=1)) if len(v) > 1 else np.nan, '/'.join(f'{x:.4f}' for x in v)


def row(model, scope, **kw):
    r = dict(model=model, scope=scope); r.update(kw); return r


def main():
    rows = []
    # ---------------- MC-ContrastiveVI fresh fits
    a1 = pd.read_csv(P / 'A1_fresh_metrics.csv'); pr = pd.read_csv(P / 'A1_latent_probes_fresh.csv')
    c3m = pd.read_csv(P / 'C3_mse_both_scales.csv'); ck1 = pd.read_csv(P / 'A1_cross_seed_cka.csv')
    def probe(rep, task):
        return pr[(pr.representation == rep) & (pr.task == task)].sort_values('seed').balanced_accuracy.tolist()
    mc = c3m[(c3m.method == 'MC-ContrastiveVI (fresh)') & (c3m.study_or_pair == 'combined')]
    au = [f'sh {a}/dox {b}/cis {c}' for a, b, c in zip(a1.shared_active_units_var_gt_0p01, a1.dox_active_units_var_gt_0p01, a1.cis_active_units_var_gt_0p01)]
    kl = [f'sh {a:.4f}/dox {b:.4f}/cis {c:.4f}' for a, b, c in zip(a1.shared_mean_total_kl_nats, a1.dox_mean_total_kl_nats, a1.cis_mean_total_kl_nats)]
    m_, s_, l_ = ms(a1.test_mse); mm_, _, _ = ms(mc.modeled_gene_log_mse)
    p1 = ms(probe('drug_ungated', 'cisplatin_treatment')); p2 = ms(probe('drug_ungated', 'doxorubicin_treatment'))
    q1 = ms(probe('shared_ungated', 'cisplatin_treatment')); q2 = ms(probe('shared_ungated', 'doxorubicin_treatment'))
    cka = ck1.groupby('space').fresh_test_cka.agg(['mean', 'min', 'max'])
    rows.append(row('MC-ContrastiveVI (fresh fits)', 'all canonical test cells', n_seeds=3, canonical_mse_mean=m_, canonical_mse_sd=s_, canonical_mse_per_seed=l_,
                    modeled_gene_mse_mean=mm_, primary_representation='drug block, ungated posterior means (8 dims)',
                    probe_cis_within_GSE216146_mean=p1[0], probe_cis_within_GSE216146_per_seed=p1[2],
                    probe_dox_within_GSE271055_mean=p2[0], probe_dox_within_GSE271055_per_seed=p2[2],
                    secondary_representation='shared block, ungated posterior means (8 dims)',
                    secondary_probe_cis_mean=q1[0], secondary_probe_cis_per_seed=q1[2], secondary_probe_dox_mean=q2[0], secondary_probe_dox_per_seed=q2[2],
                    treatment_block_active_units_per_seed='; '.join(au), treatment_block_total_kl_nats_per_seed='; '.join(kl),
                    cross_seed_cka=f"shared(gated) {cka.loc['shared','mean']:.3f} [{cka.loc['shared','min']:.3f}-{cka.loc['shared','max']:.3f}]; drug(gated) {cka.loc['drug','mean']:.3f} [{cka.loc['drug','min']:.3f}-{cka.loc['drug','max']:.3f}]",
                    gated_label_masked_probe='', source_files='A1_fresh_metrics.csv; A1_latent_probes_fresh.csv; A1_cross_seed_cka.csv; C3_mse_both_scales.csv',
                    notes='Active units: Var[E q(z|x)] > 0.01 on treated test cells; CKA on gated posterior means as stored in A1_cross_seed_cka.csv.'))
    # ---------------- PCA(32)
    pca = c3m[(c3m.method == 'PCA(32) fit on canonical scale') & (c3m.study_or_pair == 'combined')].iloc[0]
    pcam = c3m[(c3m.method == 'PCA(32) fit on modeled-gene scale') & (c3m.study_or_pair == 'combined')].iloc[0]
    rows.append(row('PCA(32)', 'all canonical test cells', n_seeds=1, canonical_mse_mean=pca.canonical_log_mse, canonical_mse_per_seed=f'{pca.canonical_log_mse:.4f}',
                    modeled_gene_mse_mean=pca.modeled_gene_log_mse, primary_representation='not available (no within-study probe in existing files)',
                    treatment_block_active_units_per_seed='not applicable', treatment_block_total_kl_nats_per_seed='not applicable', cross_seed_cka='not applicable (deterministic)',
                    gated_label_masked_probe='', source_files='C3_mse_both_scales.csv',
                    notes=f'Modeled-gene MSE of PCA(32) fitted directly on modeled-gene-scale data: {pcam.modeled_gene_log_mse:.4f}.'))
    # ---------------- scVI
    sv = pd.read_csv(P / 'C3_scvi_metrics.csv'); ck3 = pd.read_csv(P / 'C3_cross_seed_cka.csv')
    k = ck3[(ck3.method == 'scVI')].linear_cka
    a, b, c = ms(sv.canonical_log_mse), ms(sv.modeled_gene_log_mse), ms(sv.treatment_probe_cisplatin_within_GSE216146_balanced_accuracy)
    d = ms(sv.treatment_probe_doxorubicin_within_GSE271055_balanced_accuracy)
    rows.append(row('scVI', 'all canonical test cells', n_seeds=3, canonical_mse_mean=a[0], canonical_mse_sd=a[1], canonical_mse_per_seed=a[2], modeled_gene_mse_mean=b[0],
                    primary_representation='latent (32 dims; single block)', probe_cis_within_GSE216146_mean=c[0], probe_cis_within_GSE216146_per_seed=c[2],
                    probe_dox_within_GSE271055_mean=d[0], probe_dox_within_GSE271055_per_seed=d[2],
                    treatment_block_active_units_per_seed='not computed in existing files', treatment_block_total_kl_nats_per_seed='not computed in existing files',
                    cross_seed_cka=f'latent {k.mean():.3f} [{k.min():.3f}-{k.max():.3f}]', gated_label_masked_probe='', source_files='C3_scvi_metrics.csv; C3_cross_seed_cka.csv',
                    notes='No treatment-specific block; the probe uses the whole latent.'))
    # ---------------- contrastiveVI per pair
    cv = pd.read_csv(P / 'C3_contrastivevi_metrics.csv')
    for pair, study_col in (('cisplatin', 'cis'), ('doxorubicin', 'dox')):
        g = cv[cv.study_or_pair == pair].sort_values('seed'); kk = ck3[(ck3.method == 'contrastiveVI') & (ck3.study_or_pair == pair)]
        a, b = ms(g.canonical_log_mse), ms(g.modeled_gene_log_mse); sp = ms(g.salient_condition_balanced_accuracy); bp = ms(g.background_condition_balanced_accuracy)
        ks = kk[kk.representation == 'salient'].linear_cka; kt = kk[kk.representation == 'salient_target_cells_only'].linear_cka
        r = row(f'contrastiveVI ({pair} pair)', f'{pair}-pair canonical test cells', n_seeds=3, canonical_mse_mean=a[0], canonical_mse_sd=a[1], canonical_mse_per_seed=a[2],
                modeled_gene_mse_mean=b[0], primary_representation='salient (8 dims, ungated posterior means)',
                secondary_representation='background (16 dims)',
                treatment_block_active_units_per_seed='; '.join(f'{x}/8' for x in g.salient_target_test_active_units_var_gt_0p01),
                treatment_block_total_kl_nats_per_seed='; '.join(f'{x:.3f}' for x in g.salient_target_test_total_kl_nats),
                cross_seed_cka=f'salient {ks.mean():.3f} [{ks.min():.3f}-{ks.max():.3f}]; salient target cells {kt.mean():.3f} [{kt.min():.3f}-{kt.max():.3f}]',
                gated_label_masked_probe='', source_files='C3_contrastivevi_metrics.csv; C3_cross_seed_cka.csv',
                notes='Pairwise model fitted within one study; active units/KL on held-out target (treated) cells.')
        if pair == 'cisplatin':
            r.update(probe_cis_within_GSE216146_mean=sp[0], probe_cis_within_GSE216146_per_seed=sp[2], secondary_probe_cis_mean=bp[0], secondary_probe_cis_per_seed=bp[2])
        else:
            r.update(probe_dox_within_GSE271055_mean=sp[0], probe_dox_within_GSE271055_per_seed=sp[2], secondary_probe_dox_mean=bp[0], secondary_probe_dox_per_seed=bp[2])
        rows.append(r)
    # ---------------- multiGroupVI
    mg = pd.read_csv(P / 'C4_multigroupvi_metrics.csv').sort_values('seed'); ck4 = pd.read_csv(P / 'C4_multigroupvi_cross_seed_cka.csv')
    a, b = ms(mg.canonical_log_mse), ms(mg.modeled_gene_log_mse)
    c = ms(mg.private_ungated_treatment_probe_cisplatin_within_GSE216146_balanced_accuracy); d = ms(mg.private_ungated_treatment_probe_doxorubicin_within_GSE271055_balanced_accuracy)
    e = ms(mg.shared_treatment_probe_cisplatin_within_GSE216146_balanced_accuracy); f = ms(mg.shared_treatment_probe_doxorubicin_within_GSE271055_balanced_accuracy)
    gc = ms(mg.private_gated_treatment_probe_cisplatin_within_GSE216146_balanced_accuracy); gd = ms(mg.private_gated_treatment_probe_doxorubicin_within_GSE271055_balanced_accuracy)
    au = [f'cis {x}/10, dox {y}/10' for x, y in zip(mg['private_cisplatin__cisplatin_test__active_units_var_gt_0p01'], mg['private_doxorubicin__doxorubicin_test__active_units_var_gt_0p01'])]
    kl = [f'cis {x:.3f}, dox {y:.3f}' for x, y in zip(mg['private_cisplatin__cisplatin_test__total_kl_nats'], mg['private_doxorubicin__doxorubicin_test__total_kl_nats'])]
    k1 = ck4[ck4.representation == 'private_ungated'].linear_cka; k2 = ck4[ck4.representation == 'shared'].linear_cka
    rows.append(row('multiGroupVI', 'all canonical test cells', n_seeds=3, canonical_mse_mean=a[0], canonical_mse_sd=a[1], canonical_mse_per_seed=a[2], modeled_gene_mse_mean=b[0],
                    primary_representation='group-specific latents, ungated posterior means (3 x 10 dims)',
                    probe_cis_within_GSE216146_mean=c[0], probe_cis_within_GSE216146_per_seed=c[2], probe_dox_within_GSE271055_mean=d[0], probe_dox_within_GSE271055_per_seed=d[2],
                    secondary_representation='shared latent (10 dims)', secondary_probe_cis_mean=e[0], secondary_probe_cis_per_seed=e[2],
                    secondary_probe_dox_mean=f[0], secondary_probe_dox_per_seed=f[2],
                    treatment_block_active_units_per_seed='; '.join(au), treatment_block_total_kl_nats_per_seed='; '.join(kl),
                    cross_seed_cka=f'group-specific ungated {k1.mean():.3f} [{k1.min():.3f}-{k1.max():.3f}]; shared {k2.mean():.3f} [{k2.min():.3f}-{k2.max():.3f}]',
                    gated_label_masked_probe=f'cis {gc[2]} (mean {gc[0]:.4f}); dox {gd[2]} (mean {gd[0]:.4f})',
                    source_files='C4_multigroupvi_metrics.csv; C4_multigroupvi_cross_seed_cka.csv',
                    notes='gated_label_masked_probe: group-specific latents multiplied by the group-label mask; gating by group label makes these values structurally label-informative. Active units/KL on held-out cells of the block\'s own group.'))
    # ---------------- scDisInFact (committed Task 6)
    sd = pd.read_csv(C / 'task6_scdisinfact_metrics.csv').sort_values('seed'); sk = pd.read_csv(C / 'task6_scdisinfact_seed_cka.csv')
    a = ms(sd.test_log_expression_mse); u = ms(sd.unshared_bio_condition_balanced_accuracy)
    k1 = sk[sk.representation == 'unshared_bio'].linear_cka; k2 = sk[sk.representation == 'shared_bio'].linear_cka
    rows.append(row('scDisInFact (committed Task 6)', 'all canonical test cells', n_seeds=3, canonical_mse_mean=a[0], canonical_mse_sd=a[1], canonical_mse_per_seed=a[2],
                    modeled_gene_mse_mean=np.nan, primary_representation='unshared-bio (condition) latent; within-study probes not available in committed files',
                    treatment_block_active_units_per_seed='not available in committed files', treatment_block_total_kl_nats_per_seed='not available in committed files',
                    cross_seed_cka=f'unshared_bio {k1.mean():.3f} [{k1.min():.3f}-{k1.max():.3f}]; shared_bio {k2.mean():.3f} [{k2.min():.3f}-{k2.max():.3f}]',
                    gated_label_masked_probe='', source_files='analysis/corrected_20260920/task6_scdisinfact_metrics.csv; task6_scdisinfact_seed_cka.csv',
                    notes=f'Only pooled 3-class condition probes exist: unshared_bio {u[2]} (mean {u[0]:.4f}). Canonical-scale MSE from predicted raw counts with full-library normalization (task6 design).'))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / 'R2_baseline_headline.csv', index=False)
    print(df.T.to_string())


if __name__ == '__main__':
    sys.exit(main())
