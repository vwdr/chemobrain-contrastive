"""
Train MC-ContrastiveVI on the merged chemo-brain AnnData.

Usage:
    python scripts/02_train.py --config configs/default.yaml

Outputs to runs/<timestamp>/:
    - model.pt               (best checkpoint)
    - config.yaml            (frozen config)
    - train.log
    - latents.h5ad           (AnnData with obsm['z_bg'], obsm['z_shared'], obsm['z_drug'])
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import torch
import yaml
from torch.utils.data import DataLoader, TensorDataset, random_split

# make src importable
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.models.mc_contrastive_vi import (  # noqa: E402
    MCContrastiveVI,
    MCContrastiveVIConfig,
    mc_contrastive_loss,
)


def load_config(path: Path) -> dict:
    with open(path) as f:
        return yaml.safe_load(f)


def build_tensors(adata, cfg: dict):
    """Convert AnnData -> (x, drug_idx, batch_idx) float tensors."""
    import scipy.sparse as sp

    drugs = cfg["data"]["drugs"]
    drug_to_idx = {d: i for i, d in enumerate(drugs)}

    # Expression matrix
    X = adata.X
    if sp.issparse(X):
        X = X.toarray()
    x = torch.tensor(X, dtype=torch.float32)

    # Condition labels
    obs_drug = adata.obs[cfg["data"]["condition_key"]].astype(str).str.lower().str.strip()
    # map unknown drugs to error
    missing = set(obs_drug.unique()) - set(drug_to_idx)
    if missing:
        raise ValueError(
            f"Drugs present in data but not in config.data.drugs: {missing}. "
            f"Expected: {drugs}"
        )
    drug_idx = torch.tensor([drug_to_idx[d] for d in obs_drug.tolist()], dtype=torch.long)

    # Batch labels (dataset_id)
    obs_batch = adata.obs[cfg["data"]["batch_key"]].astype(str)
    batches = sorted(obs_batch.unique())
    batch_to_idx = {b: i for i, b in enumerate(batches)}
    batch_idx = torch.tensor([batch_to_idx[b] for b in obs_batch.tolist()], dtype=torch.long)

    return x, drug_idx, batch_idx, drug_to_idx, batch_to_idx


def pick_device(pref: str) -> torch.device:
    if pref == "cpu":
        return torch.device("cpu")
    if pref == "cuda":
        return torch.device("cuda" if torch.cuda.is_available() else "cpu")
    if pref == "mps":
        return torch.device("mps" if torch.backends.mps.is_available() else "cpu")
    # auto
    if torch.cuda.is_available():
        return torch.device("cuda")
    if torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", type=Path, default="configs/default.yaml")
    args = ap.parse_args()
    cfg = load_config(args.config)

    # ----- output dir -----
    run_dir = Path(cfg["project"]["output_dir"]) / time.strftime("%Y%m%d_%H%M%S")
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "config.yaml").write_text(yaml.dump(cfg, sort_keys=False))

    torch.manual_seed(cfg["project"]["seed"])
    np.random.seed(cfg["project"]["seed"])

    # ----- data -----
    import anndata as ad
    adata = ad.read_h5ad(cfg["data"]["adata_path"])
    x, drug_idx, batch_idx, drug_to_idx, batch_to_idx = build_tensors(adata, cfg)
    n_genes = x.shape[1]
    n_drugs = len(drug_to_idx)
    n_batches = len(batch_to_idx)
    print(f"[data] {x.shape[0]} cells, {n_genes} genes, {n_drugs} drugs, {n_batches} batches")

    ds = TensorDataset(x, drug_idx, batch_idx)
    val_size = int(cfg["training"]["val_frac"] * len(ds))
    train_size = len(ds) - val_size
    train_ds, val_ds = random_split(
        ds, [train_size, val_size],
        generator=torch.Generator().manual_seed(cfg["project"]["seed"]),
    )
    train_loader = DataLoader(train_ds, batch_size=cfg["training"]["batch_size"], shuffle=True)
    val_loader = DataLoader(val_ds, batch_size=cfg["training"]["batch_size"])

    # ----- model -----
    m = cfg["model"]
    model_cfg = MCContrastiveVIConfig(
        n_genes=n_genes, n_drugs=n_drugs, n_batches=n_batches,
        hidden_dim=m["hidden_dim"], n_hidden_layers=m["n_hidden_layers"],
        latent_dim_bg=m["latent_dim_bg"], latent_dim_shared=m["latent_dim_shared"],
        latent_dim_drug=m["latent_dim_drug"], dropout=m["dropout"],
        use_batch_onehot=m["use_batch_onehot"],
    )
    device = pick_device(cfg["training"]["device"])
    model = MCContrastiveVI(model_cfg).to(device)
    opt = torch.optim.AdamW(
        model.parameters(),
        lr=cfg["training"]["lr"],
        weight_decay=cfg["training"]["weight_decay"],
    )

    # ----- train -----
    best_val = float("inf")
    patience = cfg["training"]["early_stopping_patience"]
    bad_epochs = 0
    log_path = run_dir / "train.log"
    log_path.write_text("epoch\tsplit\tloss\trecon\tkl_bg\tkl_sh\tkl_dr\thsic_bg_sh\thsic_bg_dr\thsic_sh_dr\n")

    for epoch in range(cfg["training"]["epochs"]):
        model.train()
        train_metrics = {k: 0.0 for k in
                         ["loss", "recon", "kl_bg", "kl_sh", "kl_dr",
                          "hsic_bg_sh", "hsic_bg_dr", "hsic_sh_dr"]}
        n_train = 0
        for xb, db, bb in train_loader:
            xb, db, bb = xb.to(device), db.to(device), bb.to(device)
            enc = model(xb, db, bb)
            losses = mc_contrastive_loss(
                model, enc, xb, db,
                kl_weight=m["kl_weight"],
                hsic_weight_bg_shared=m["hsic_weight_bg_shared"],
                hsic_weight_bg_drug=m["hsic_weight_bg_drug"],
                hsic_weight_shared_drug=m["hsic_weight_shared_drug"],
            )
            opt.zero_grad()
            losses["loss"].backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            opt.step()
            bs = xb.shape[0]
            n_train += bs
            for k, v in losses.items():
                key = {"recon_nll": "recon", "kl_shared": "kl_sh", "kl_drug": "kl_dr",
                       "hsic_bg_shared": "hsic_bg_sh", "hsic_bg_drug": "hsic_bg_dr",
                       "hsic_shared_drug": "hsic_sh_dr", "loss": "loss", "kl_bg": "kl_bg"}.get(k)
                if key is not None:
                    train_metrics[key] += float(v) * bs

        train_metrics = {k: v / n_train for k, v in train_metrics.items()}

        # validation
        model.eval()
        val_loss, n_val = 0.0, 0
        with torch.no_grad():
            for xb, db, bb in val_loader:
                xb, db, bb = xb.to(device), db.to(device), bb.to(device)
                enc = model(xb, db, bb)
                losses = mc_contrastive_loss(
                    model, enc, xb, db,
                    kl_weight=m["kl_weight"],
                    hsic_weight_bg_shared=m["hsic_weight_bg_shared"],
                    hsic_weight_bg_drug=m["hsic_weight_bg_drug"],
                    hsic_weight_shared_drug=m["hsic_weight_shared_drug"],
                )
                bs = xb.shape[0]
                val_loss += float(losses["loss"]) * bs
                n_val += bs
        val_loss /= max(n_val, 1)

        print(f"[{epoch:03d}] train_loss={train_metrics['loss']:.3f}  val_loss={val_loss:.3f}")
        with open(log_path, "a") as f:
            f.write(f"{epoch}\ttrain\t" + "\t".join(
                f"{train_metrics[k]:.4f}" for k in
                ["loss", "recon", "kl_bg", "kl_sh", "kl_dr",
                 "hsic_bg_sh", "hsic_bg_dr", "hsic_sh_dr"]) + "\n")
            f.write(f"{epoch}\tval\t{val_loss:.4f}\n")

        if val_loss < best_val - 1e-3:
            best_val = val_loss
            bad_epochs = 0
            torch.save({
                "model_state": model.state_dict(),
                "model_cfg": model_cfg.__dict__,
                "drug_to_idx": drug_to_idx,
                "batch_to_idx": batch_to_idx,
                "epoch": epoch,
                "val_loss": val_loss,
            }, run_dir / "model.pt")
        else:
            bad_epochs += 1
            if bad_epochs >= patience:
                print(f"[early stop] no improvement for {patience} epochs")
                break

    # ----- export latents for analysis -----
    print("[export] computing latents for all cells ...")
    ckpt = torch.load(run_dir / "model.pt", map_location=device, weights_only=False)
    model.load_state_dict(ckpt["model_state"])
    model.eval()
    all_bg, all_sh, all_dr = [], [], []
    with torch.no_grad():
        loader = DataLoader(ds, batch_size=cfg["training"]["batch_size"])
        for xb, db, bb in loader:
            xb, db, bb = xb.to(device), db.to(device), bb.to(device)
            enc = model.encode(xb, db)
            # use means for reproducibility
            all_bg.append(enc["mu_bg"].cpu().numpy())
            sh_gate = (db > 0).float().unsqueeze(1).cpu().numpy()
            all_sh.append(enc["mu_shared"].cpu().numpy() * sh_gate)
            # drug-specific: gated means
            dr_parts = []
            for k, mu_dr in enumerate(enc["mu_drug_list"]):
                gate = (db == (k + 1)).float().unsqueeze(1).cpu().numpy()
                dr_parts.append(mu_dr.cpu().numpy() * gate)
            all_dr.append(np.concatenate(dr_parts, axis=1))

    adata.obsm["z_bg"] = np.concatenate(all_bg)
    adata.obsm["z_shared"] = np.concatenate(all_sh)
    adata.obsm["z_drug"] = np.concatenate(all_dr)
    adata.write_h5ad(run_dir / "latents.h5ad")

    (run_dir / "result.json").write_text(json.dumps({
        "best_val_loss": best_val, "best_epoch": ckpt["epoch"]
    }, indent=2))
    print(f"[done] artifacts in {run_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
