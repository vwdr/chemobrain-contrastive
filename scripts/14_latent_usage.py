"""Held-out latent-utilization diagnostics for corrected full MC-ContrastiveVI fits.

Computes per-dimension KL divergence and posterior-mean variance on the treated
held-out test cells for the canonical full-model seeds. Active units follow the
common criterion Var_x[E_q(z_j|x)] > 0.01. A KL > 0.01 nats indicator is also
reported as a sensitivity diagnostic, but is not labeled as the active-unit
definition.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch

R = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(R))

import src.models.mc_contrastive_vi as mm

RUN = R / "runs" / "corrected_20260920"
OUT = R / "analysis" / "corrected_20260920"
SEEDS = (0, 1, 2)
AU_VAR_THRESHOLD = 0.01
KL_SENSITIVITY_THRESHOLD = 0.01
BATCH_SIZE = 512


def per_dim_kl(mu: np.ndarray, logvar: np.ndarray) -> np.ndarray:
    """KL(q(z_j|x)||N(0,1)) for every cell x latent dimension j."""
    return 0.5 * (mu**2 + np.exp(logvar) - 1.0 - logvar)


def load_posteriors(model, X: torch.Tensor, d: torch.Tensor):
    buckets = {
        "shared_mu": [],
        "shared_lv": [],
        "dox_mu": [],
        "dox_lv": [],
        "cis_mu": [],
        "cis_lv": [],
    }
    model.eval()
    with torch.no_grad():
        for ix in np.array_split(np.arange(len(X)), int(np.ceil(len(X) / BATCH_SIZE))):
            e = model.encode(X[ix], d[ix])
            buckets["shared_mu"].append(e["mu_shared"].cpu().numpy())
            buckets["shared_lv"].append(e["lv_shared"].cpu().numpy())
            buckets["dox_mu"].append(e["mu_drug_list"][0].cpu().numpy())
            buckets["dox_lv"].append(e["lv_drug_list"][0].cpu().numpy())
            buckets["cis_mu"].append(e["mu_drug_list"][1].cpu().numpy())
            buckets["cis_lv"].append(e["lv_drug_list"][1].cpu().numpy())
    return {k: np.concatenate(v, axis=0) for k, v in buckets.items()}


def variance_fraction(shared_mu: np.ndarray, active_drug_mu: np.ndarray) -> float:
    vs = np.var(shared_mu, axis=0).sum()
    vd = np.var(active_drug_mu, axis=0).sum()
    return float(vs / (vs + vd))


def append_dim_rows(rows, seed, block, condition, mu, lv):
    kl = per_dim_kl(mu, lv)
    mean_kl = kl.mean(axis=0)
    mu_var = np.var(mu, axis=0)
    for j in range(mu.shape[1]):
        rows.append(
            {
                "seed": seed,
                "block": block,
                "condition": condition,
                "dimension": j,
                "n_cells": len(mu),
                "mean_kl_nats": float(mean_kl[j]),
                "posterior_mean_variance": float(mu_var[j]),
                "active_unit_var_gt_0p01": bool(mu_var[j] > AU_VAR_THRESHOLD),
                "kl_gt_0p01": bool(mean_kl[j] > KL_SENSITIVITY_THRESHOLD),
            }
        )
    return mean_kl, mu_var


def main():
    z = np.load(RUN / "inputs.npz")
    X = torch.from_numpy(z["X"])
    d = torch.from_numpy(z["d"])
    test = z["test"]

    treated = test[z["d"][test] > 0]
    dox = test[z["d"][test] == 1]
    cis = test[z["d"][test] == 2]

    assert len(treated) == len(dox) + len(cis)
    assert len(dox) > 0 and len(cis) > 0

    dim_rows = []
    summary_rows = []

    for seed in SEEDS:
        ckpt_path = RUN / f"full_{seed}.pt"
        metrics_path = RUN / f"full_{seed}_metrics.json"
        assert ckpt_path.exists(), f"Missing checkpoint: {ckpt_path}"

        checkpoint = torch.load(ckpt_path, map_location="cpu")
        cfg = mm.MCContrastiveVIConfig(**checkpoint["config"])
        model = mm.MCContrastiveVI(cfg)
        model.load_state_dict(checkpoint["state"])
        post = load_posteriors(model, X, d)

        sh_mu = post["shared_mu"][treated]
        sh_lv = post["shared_lv"][treated]
        dox_mu = post["dox_mu"][dox]
        dox_lv = post["dox_lv"][dox]
        cis_mu = post["cis_mu"][cis]
        cis_lv = post["cis_lv"][cis]

        sh_kl_dim, sh_var_dim = append_dim_rows(
            dim_rows, seed, "shared", "all_treated", sh_mu, sh_lv
        )
        dox_kl_dim, dox_var_dim = append_dim_rows(
            dim_rows, seed, "drug_specific", "doxorubicin", dox_mu, dox_lv
        )
        cis_kl_dim, cis_var_dim = append_dim_rows(
            dim_rows, seed, "drug_specific", "cisplatin", cis_mu, cis_lv
        )

        # Recreate the gated active drug matrix used by the published F_sh metric.
        active_drug = np.zeros((len(treated), cfg.latent_dim_drug * 2), dtype=np.float32)
        treated_d = z["d"][treated]
        dox_pos = np.where(treated_d == 1)[0]
        cis_pos = np.where(treated_d == 2)[0]
        active_drug[dox_pos, : cfg.latent_dim_drug] = post["dox_mu"][treated[dox_pos]]
        active_drug[cis_pos, cfg.latent_dim_drug :] = post["cis_mu"][treated[cis_pos]]

        pooled_fsh = variance_fraction(sh_mu, active_drug)
        dox_fsh = variance_fraction(
            post["shared_mu"][dox],
            post["dox_mu"][dox],
        )
        cis_fsh = variance_fraction(
            post["shared_mu"][cis],
            post["cis_mu"][cis],
        )

        shared_kl_total = float(sh_kl_dim.sum())
        dox_kl_total = float(dox_kl_dim.sum())
        cis_kl_total = float(cis_kl_dim.sum())

        objective_ratio = shared_kl_total / (
            shared_kl_total + dox_kl_total + cis_kl_total
        )

        # Observation-weighted active-drug KL on the treated test cells.
        sh_cell_total = per_dim_kl(sh_mu, sh_lv).sum(axis=1)
        active_drug_cell_total = np.empty(len(treated), dtype=np.float64)
        active_drug_cell_total[dox_pos] = per_dim_kl(
            post["dox_mu"][treated[dox_pos]], post["dox_lv"][treated[dox_pos]]
        ).sum(axis=1)
        active_drug_cell_total[cis_pos] = per_dim_kl(
            post["cis_mu"][treated[cis_pos]], post["cis_lv"][treated[cis_pos]]
        ).sum(axis=1)
        obs_weighted_ratio = float(
            sh_cell_total.mean()
            / (sh_cell_total.mean() + active_drug_cell_total.mean())
        )

        committed = json.loads(metrics_path.read_text())
        committed_fsh = float(committed["shared_fraction"])
        delta = pooled_fsh - committed_fsh
        if abs(delta) > 1e-6:
            raise RuntimeError(
                f"Seed {seed}: recomputed F_sh={pooled_fsh:.9f} does not match "
                f"metrics JSON {committed_fsh:.9f} (delta={delta:.3g})"
            )

        summary_rows.append(
            {
                "seed": seed,
                "n_treated_test": len(treated),
                "n_doxorubicin_test": len(dox),
                "n_cisplatin_test": len(cis),
                "shared_fraction_recomputed": pooled_fsh,
                "shared_fraction_metrics_json": committed_fsh,
                "doxorubicin_shared_fraction": dox_fsh,
                "cisplatin_shared_fraction": cis_fsh,
                "shared_active_units_var_gt_0p01": int(
                    (sh_var_dim > AU_VAR_THRESHOLD).sum()
                ),
                "shared_total_dims": len(sh_var_dim),
                "dox_active_units_var_gt_0p01": int(
                    (dox_var_dim > AU_VAR_THRESHOLD).sum()
                ),
                "dox_total_dims": len(dox_var_dim),
                "cis_active_units_var_gt_0p01": int(
                    (cis_var_dim > AU_VAR_THRESHOLD).sum()
                ),
                "cis_total_dims": len(cis_var_dim),
                "drug_active_units_var_gt_0p01_total": int(
                    (dox_var_dim > AU_VAR_THRESHOLD).sum()
                    + (cis_var_dim > AU_VAR_THRESHOLD).sum()
                ),
                "drug_total_dims": len(dox_var_dim) + len(cis_var_dim),
                "shared_dims_kl_gt_0p01": int(
                    (sh_kl_dim > KL_SENSITIVITY_THRESHOLD).sum()
                ),
                "dox_dims_kl_gt_0p01": int(
                    (dox_kl_dim > KL_SENSITIVITY_THRESHOLD).sum()
                ),
                "cis_dims_kl_gt_0p01": int(
                    (cis_kl_dim > KL_SENSITIVITY_THRESHOLD).sum()
                ),
                "shared_mean_total_kl_nats": shared_kl_total,
                "dox_mean_total_kl_nats": dox_kl_total,
                "cis_mean_total_kl_nats": cis_kl_total,
                "objective_matched_shared_kl_ratio": float(objective_ratio),
                "observation_weighted_shared_kl_ratio": obs_weighted_ratio,
            }
        )

    dim_df = pd.DataFrame(dim_rows)
    summary_df = pd.DataFrame(summary_rows)
    OUT.mkdir(parents=True, exist_ok=True)
    dim_df.to_csv(OUT / "latent_usage_per_dimension.csv", index=False)
    summary_df.to_csv(OUT / "latent_usage_summary.csv", index=False)

    thresholds = {
        "evaluation_population": "held-out treated test cells",
        "active_unit_definition": "Var_x[E_q(z_j|x)] > 0.01",
        "active_unit_variance_threshold": AU_VAR_THRESHOLD,
        "kl_sensitivity_definition": "mean per-dimension KL > 0.01 nats",
        "kl_sensitivity_threshold_nats": KL_SENSITIVITY_THRESHOLD,
        "note": (
            "KL > 0.01 is reported only as a sensitivity diagnostic; "
            "the active-unit count uses posterior-mean variance."
        ),
    }
    (OUT / "latent_usage_thresholds.json").write_text(
        json.dumps(thresholds, indent=2) + "\n"
    )

    pd.set_option("display.max_columns", None)
    print("\n===== LATENT USAGE SUMMARY =====")
    print(summary_df.to_string(index=False))
    print("\nWrote:")
    print(OUT / "latent_usage_per_dimension.csv")
    print(OUT / "latent_usage_summary.csv")
    print(OUT / "latent_usage_thresholds.json")


if __name__ == "__main__":
    main()
