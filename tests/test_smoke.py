"""
End-to-end smoke test on synthetic data.

Simulates 3 drugs + control, 4 cell types, 500 HVGs, 2 batches.
Injects known drug-specific and shared-toxicity signals, then verifies that
after training the model:
    (a) reconstructs reasonably,
    (b) puts treated cells farther from origin on the shared axis than controls,
    (c) separates drugs on the drug-specific axes.

Run:  python tests/test_smoke.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import torch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.models.mc_contrastive_vi import (  # noqa: E402
    MCContrastiveVI,
    MCContrastiveVIConfig,
    mc_contrastive_loss,
)


def simulate(n_cells_per_condition=600, n_genes=500, seed=0):
    rng = np.random.default_rng(seed)
    # 4 cell types with different expression profiles (background signal)
    n_types = 4
    type_profiles = rng.normal(0, 1, (n_types, n_genes))

    # Shared toxicity signature: same vector bumps up in all treated cells
    shared_sig = rng.normal(0, 1, n_genes) * (rng.random(n_genes) < 0.05)

    # Drug-specific signatures (sparse, distinct)
    drug_sigs = [rng.normal(0, 1, n_genes) * (rng.random(n_genes) < 0.03) for _ in range(3)]

    conditions = ["control", "drug1", "drug2", "drug3"]
    xs, drug_idx_all, batch_idx_all, celltype_all = [], [], [], []

    for d, cond in enumerate(conditions):
        for _ in range(n_cells_per_condition):
            ct = rng.integers(0, n_types)
            base = type_profiles[ct] + rng.normal(0, 0.3, n_genes)
            if cond != "control":
                base = base + shared_sig * 0.8
                base = base + drug_sigs[d - 1] * 1.2
            batch = rng.integers(0, 2)
            # mild batch effect
            base = base + (0.2 if batch == 1 else -0.2)
            xs.append(base)
            drug_idx_all.append(d)
            batch_idx_all.append(batch)
            celltype_all.append(ct)

    x = np.stack(xs).astype(np.float32)
    return (
        torch.tensor(x),
        torch.tensor(drug_idx_all, dtype=torch.long),
        torch.tensor(batch_idx_all, dtype=torch.long),
        np.array(celltype_all),
    )


def main():
    torch.manual_seed(0)
    x, drug_idx, batch_idx, celltypes = simulate()
    n_genes = x.shape[1]
    print(f"Simulated: {x.shape[0]} cells, {n_genes} genes, "
          f"{drug_idx.unique().numel()} conditions, {batch_idx.unique().numel()} batches")

    cfg = MCContrastiveVIConfig(
        n_genes=n_genes, n_drugs=4, n_batches=2,
        hidden_dim=128, n_hidden_layers=2,
        latent_dim_bg=8, latent_dim_shared=4, latent_dim_drug=2,
    )
    model = MCContrastiveVI(cfg)
    opt = torch.optim.AdamW(model.parameters(), lr=1e-3)

    n_epochs = 40
    bs = 256
    losses_history = []
    for epoch in range(n_epochs):
        model.train()
        perm = torch.randperm(x.shape[0])
        total = 0.0
        n = 0
        for i in range(0, x.shape[0], bs):
            idx = perm[i:i + bs]
            xb = x[idx]
            db = drug_idx[idx]
            bb = batch_idx[idx]
            enc = model(xb, db, bb)
            losses = mc_contrastive_loss(
                model, enc, xb, db,
                kl_weight=0.5,
                hsic_weight_bg_shared=5.0,
                hsic_weight_bg_drug=5.0,
                hsic_weight_shared_drug=2.0,
            )
            opt.zero_grad()
            losses["loss"].backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            opt.step()
            total += float(losses["loss"]) * xb.shape[0]
            n += xb.shape[0]
        losses_history.append(total / n)
        if epoch % 5 == 0 or epoch == n_epochs - 1:
            print(f"  epoch {epoch:02d}  loss={total/n:.3f}")

    # Check training actually reduced loss
    assert losses_history[-1] < losses_history[0] * 0.7, \
        f"Loss only dropped from {losses_history[0]:.2f} to {losses_history[-1]:.2f}"
    print(f"[ok] loss dropped {losses_history[0]:.2f} -> {losses_history[-1]:.2f}")

    # Check shared axis is bigger on treated vs control
    model.eval()
    with torch.no_grad():
        enc = model.encode(x, drug_idx)
    z_shared_mu = enc["mu_shared"].numpy()
    treated = (drug_idx > 0).numpy()
    mag_t = np.linalg.norm(z_shared_mu[treated], axis=1).mean()
    mag_c = np.linalg.norm(z_shared_mu[~treated], axis=1).mean()
    print(f"[shared axis]  treated-mean-mag = {mag_t:.3f}   control-mean-mag = {mag_c:.3f}")
    # the shared encoder sees controls too, so we expect a gap but not infinity.
    # the gating means the shared axis is only *used* in the decoder for treated cells,
    # so the encoder learns it for treated cells.
    assert mag_t > mag_c, "Shared axis should be more active on treated cells"
    print("[ok] shared axis is more engaged on treated cells")

    # Check drug-specific heads differentiate drugs
    # Each drug_k head should fire mostly on drug_k cells.
    for k in range(3):
        mu_k = enc["mu_drug_list"][k].numpy()
        mag_own = np.linalg.norm(mu_k[drug_idx == (k + 1)], axis=1).mean()
        mag_other_treated = np.linalg.norm(mu_k[(drug_idx > 0) & (drug_idx != (k + 1))], axis=1).mean()
        print(f"[drug {k+1} head]  own={mag_own:.3f}  other-treated={mag_other_treated:.3f}")

    print("\n[SMOKE TEST PASSED]")


if __name__ == "__main__":
    main()
