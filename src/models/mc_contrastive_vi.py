"""
Multi-Condition Contrastive Variational Inference (MC-ContrastiveVI).

Extends the contrastiveVI framework (Weinberger et al., Nature Methods 2023) from
a single treatment-vs-control setup to a multi-drug setup, learning:

    z_bg       : background latent  (shared biological variance across all cells)
    z_shared   : shared-toxicity latent  (common to ALL drugs vs control)
    z_drug[k]  : drug-specific latent for drug k  (active only when cell is in drug k)

Gating rule (hard):
    control cells:          z_shared = 0, z_drug[k] = 0 for all k
    drug-k cells:           z_shared = active, z_drug[k] = active, z_drug[j!=k] = 0

Independence is encouraged via HSIC penalties between z_bg and the salient latents,
and (weaker) between z_shared and z_drug[k].

This implementation uses a Gaussian likelihood on log-normalized HVG expression for
clarity. For a production run, swap `_decode` to a ZINB head and use raw counts
(see scvi-tools' contrastive_vi for a count-based reference implementation).

Author: <your name>
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, Tuple

import torch
import torch.nn as nn
import torch.nn.functional as F


# -----------------------------------------------------------------------------
# Config
# -----------------------------------------------------------------------------

@dataclass
class MCContrastiveVIConfig:
    n_genes: int
    n_drugs: int                  # includes control at index 0
    n_batches: int = 1            # dataset_id cardinality for batch correction
    hidden_dim: int = 256
    n_hidden_layers: int = 2
    latent_dim_bg: int = 16
    latent_dim_shared: int = 8
    latent_dim_drug: int = 4      # per-drug; total salient = latent_dim_drug * (n_drugs-1)
    dropout: float = 0.1
    use_batch_onehot: bool = True


# -----------------------------------------------------------------------------
# Building blocks
# -----------------------------------------------------------------------------

def _mlp(in_dim: int, hidden_dim: int, n_layers: int, dropout: float) -> nn.Sequential:
    layers = []
    d = in_dim
    for _ in range(n_layers):
        layers += [nn.Linear(d, hidden_dim), nn.LayerNorm(hidden_dim),
                   nn.ReLU(), nn.Dropout(dropout)]
        d = hidden_dim
    return nn.Sequential(*layers)


class GaussianLatentHead(nn.Module):
    """Two linear heads producing (mu, logvar) from a trunk feature."""

    def __init__(self, in_dim: int, out_dim: int):
        super().__init__()
        self.mu = nn.Linear(in_dim, out_dim)
        self.logvar = nn.Linear(in_dim, out_dim)

    def forward(self, h: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor]:
        return self.mu(h), self.logvar(h).clamp(min=-10.0, max=10.0)


def reparameterize(mu: torch.Tensor, logvar: torch.Tensor) -> torch.Tensor:
    std = torch.exp(0.5 * logvar)
    eps = torch.randn_like(std)
    return mu + eps * std


# -----------------------------------------------------------------------------
# The model
# -----------------------------------------------------------------------------

class MCContrastiveVI(nn.Module):
    """
    Multi-condition contrastive VI.

    Inputs
    ------
    x : (B, G)       log-normalized expression (HVGs)
    drug_idx : (B,)  long tensor, 0 == control, 1..K-1 == drugs
    batch_idx: (B,)  long tensor of dataset_id (for batch-aware decoding)

    Outputs (via .forward)
    ----------------------
    dict with keys:
        recon        : (B, G)
        z_bg_mu, z_bg_logvar, z_bg
        z_shared_mu, z_shared_logvar, z_shared_gated
        z_drug_mu, z_drug_logvar, z_drug_gated   # per-drug concatenated
    """

    def __init__(self, cfg: MCContrastiveVIConfig):
        super().__init__()
        self.cfg = cfg
        assert cfg.n_drugs >= 2, "Need at least control + 1 drug."

        # ---- Encoders ----
        # background encoder sees only expression
        self.bg_trunk = _mlp(cfg.n_genes, cfg.hidden_dim, cfg.n_hidden_layers, cfg.dropout)
        self.bg_head = GaussianLatentHead(cfg.hidden_dim, cfg.latent_dim_bg)

        # shared-toxicity encoder also sees expression only
        self.shared_trunk = _mlp(cfg.n_genes, cfg.hidden_dim, cfg.n_hidden_layers, cfg.dropout)
        self.shared_head = GaussianLatentHead(cfg.hidden_dim, cfg.latent_dim_shared)

        # ONE drug-specific encoder per drug (excluding control).
        # We could share a trunk; separate trunks give more capacity.
        n_drug_conditions = cfg.n_drugs - 1
        self.drug_trunks = nn.ModuleList([
            _mlp(cfg.n_genes, cfg.hidden_dim, cfg.n_hidden_layers, cfg.dropout)
            for _ in range(n_drug_conditions)
        ])
        self.drug_heads = nn.ModuleList([
            GaussianLatentHead(cfg.hidden_dim, cfg.latent_dim_drug)
            for _ in range(n_drug_conditions)
        ])

        # ---- Decoder ----
        total_latent = (
            cfg.latent_dim_bg
            + cfg.latent_dim_shared
            + cfg.latent_dim_drug * n_drug_conditions
        )
        decoder_in = total_latent + (cfg.n_batches if cfg.use_batch_onehot else 0)
        self.decoder_trunk = _mlp(decoder_in, cfg.hidden_dim, cfg.n_hidden_layers, cfg.dropout)
        self.decoder_out = nn.Linear(cfg.hidden_dim, cfg.n_genes)
        # log-sigma for Gaussian likelihood on log-normalized expression
        self.log_sigma = nn.Parameter(torch.zeros(cfg.n_genes))

    # ------------------------------------------------------------------ encode
    def encode(
        self, x: torch.Tensor, drug_idx: torch.Tensor
    ) -> Dict[str, torch.Tensor]:
        B = x.shape[0]
        K = self.cfg.n_drugs - 1  # number of non-control drugs
        device = x.device

        # background (always active)
        h_bg = self.bg_trunk(x)
        mu_bg, lv_bg = self.bg_head(h_bg)
        z_bg = reparameterize(mu_bg, lv_bg)

        # shared-toxicity (active iff drug != control)
        h_sh = self.shared_trunk(x)
        mu_sh, lv_sh = self.shared_head(h_sh)
        z_sh = reparameterize(mu_sh, lv_sh)
        is_treated = (drug_idx > 0).float().unsqueeze(1)  # (B, 1)
        z_sh_gated = z_sh * is_treated

        # drug-specific (one-hot gated)
        mu_dr_list, lv_dr_list, z_dr_list = [], [], []
        # drug_idx == k+1  ->  k-th non-control encoder active
        for k in range(K):
            h_dr = self.drug_trunks[k](x)
            mu_dr, lv_dr = self.drug_heads[k](h_dr)
            z_dr = reparameterize(mu_dr, lv_dr)
            gate = (drug_idx == (k + 1)).float().unsqueeze(1)
            z_dr_list.append(z_dr * gate)
            mu_dr_list.append(mu_dr)
            lv_dr_list.append(lv_dr)

        z_dr_cat = torch.cat(z_dr_list, dim=1)  # (B, latent_dim_drug * K)

        return {
            "z_bg": z_bg, "mu_bg": mu_bg, "lv_bg": lv_bg,
            "z_shared": z_sh, "z_shared_gated": z_sh_gated,
            "mu_shared": mu_sh, "lv_shared": lv_sh,
            "z_drug_cat": z_dr_cat,
            "mu_drug_list": mu_dr_list, "lv_drug_list": lv_dr_list,
        }

    # ------------------------------------------------------------------ decode
    def decode(
        self,
        z_bg: torch.Tensor,
        z_shared_gated: torch.Tensor,
        z_drug_cat: torch.Tensor,
        batch_idx: torch.Tensor | None = None,
    ) -> torch.Tensor:
        parts = [z_bg, z_shared_gated, z_drug_cat]
        if self.cfg.use_batch_onehot and batch_idx is not None:
            parts.append(F.one_hot(batch_idx, self.cfg.n_batches).float())
        z = torch.cat(parts, dim=1)
        return self.decoder_out(self.decoder_trunk(z))

    # ------------------------------------------------------------------ forward
    def forward(
        self,
        x: torch.Tensor,
        drug_idx: torch.Tensor,
        batch_idx: torch.Tensor | None = None,
    ) -> Dict[str, torch.Tensor]:
        enc = self.encode(x, drug_idx)
        recon = self.decode(
            enc["z_bg"], enc["z_shared_gated"], enc["z_drug_cat"], batch_idx
        )
        enc["recon"] = recon
        return enc


# -----------------------------------------------------------------------------
# Loss components
# -----------------------------------------------------------------------------

def gaussian_nll(x: torch.Tensor, x_hat: torch.Tensor, log_sigma: torch.Tensor) -> torch.Tensor:
    """Per-sample negative log-likelihood, summed over genes."""
    var = torch.exp(2 * log_sigma)
    nll = 0.5 * ((x - x_hat) ** 2 / var + 2 * log_sigma + math.log(2 * math.pi))
    return nll.sum(dim=1)


def kl_standard_normal(mu: torch.Tensor, logvar: torch.Tensor) -> torch.Tensor:
    """KL(N(mu, sigma) || N(0, 1)), summed across latent dims."""
    return 0.5 * (mu.pow(2) + logvar.exp() - 1 - logvar).sum(dim=1)


def _rbf_kernel(x: torch.Tensor, sigma: float | None = None) -> torch.Tensor:
    """Median-heuristic RBF kernel matrix."""
    x2 = (x ** 2).sum(dim=1, keepdim=True)
    d2 = x2 + x2.T - 2 * x @ x.T
    d2 = d2.clamp(min=0.0)
    if sigma is None:
        with torch.no_grad():
            med = torch.median(d2[d2 > 0]) if (d2 > 0).any() else torch.tensor(1.0, device=x.device)
            sigma = torch.sqrt(0.5 * med).item() if med > 0 else 1.0
    return torch.exp(-d2 / (2 * sigma ** 2 + 1e-8))


def hsic(x: torch.Tensor, y: torch.Tensor) -> torch.Tensor:
    """Biased HSIC estimator. Penalizes statistical dependence between x and y."""
    n = x.shape[0]
    if n < 4:
        return torch.tensor(0.0, device=x.device)
    Kx = _rbf_kernel(x)
    Ky = _rbf_kernel(y)
    H = torch.eye(n, device=x.device) - torch.full((n, n), 1.0 / n, device=x.device)
    return (Kx @ H @ Ky @ H).diagonal().sum() / (n - 1) ** 2


# -----------------------------------------------------------------------------
# Combined objective
# -----------------------------------------------------------------------------

def mc_contrastive_loss(
    model: MCContrastiveVI,
    enc: Dict[str, torch.Tensor],
    x: torch.Tensor,
    drug_idx: torch.Tensor,
    kl_weight: float,
    hsic_weight_bg_shared: float,
    hsic_weight_bg_drug: float,
    hsic_weight_shared_drug: float,
) -> Dict[str, torch.Tensor]:
    """
    Returns a dict with 'loss' (to backprop) plus individual components for logging.
    """
    recon = enc["recon"]
    recon_nll = gaussian_nll(x, recon, model.log_sigma).mean()

    # KL for bg (always)
    kl_bg = kl_standard_normal(enc["mu_bg"], enc["lv_bg"]).mean()

    # KL for shared only on treated cells (matches contrastiveVI rule)
    treated_mask = drug_idx > 0
    if treated_mask.any():
        kl_sh = kl_standard_normal(
            enc["mu_shared"][treated_mask], enc["lv_shared"][treated_mask]
        ).mean()
    else:
        kl_sh = torch.tensor(0.0, device=x.device)

    # KL for each drug head, only on that drug's cells
    kl_drug_total = torch.tensor(0.0, device=x.device)
    K = len(enc["mu_drug_list"])
    for k in range(K):
        mask = drug_idx == (k + 1)
        if mask.any():
            kl_drug_total = kl_drug_total + kl_standard_normal(
                enc["mu_drug_list"][k][mask], enc["lv_drug_list"][k][mask]
            ).mean()

    # HSIC disentanglement (computed on mini-batch)
    hsic_bg_shared = hsic(enc["z_bg"], enc["z_shared"])
    hsic_bg_drug = hsic(enc["z_bg"], enc["z_drug_cat"])
    hsic_shared_drug = hsic(enc["z_shared"], enc["z_drug_cat"])

    loss = (
        recon_nll
        + kl_weight * (kl_bg + kl_sh + kl_drug_total)
        + hsic_weight_bg_shared * hsic_bg_shared
        + hsic_weight_bg_drug * hsic_bg_drug
        + hsic_weight_shared_drug * hsic_shared_drug
    )

    return {
        "loss": loss,
        "recon_nll": recon_nll.detach(),
        "kl_bg": kl_bg.detach(),
        "kl_shared": kl_sh.detach(),
        "kl_drug": kl_drug_total.detach(),
        "hsic_bg_shared": hsic_bg_shared.detach(),
        "hsic_bg_drug": hsic_bg_drug.detach(),
        "hsic_shared_drug": hsic_shared_drug.detach(),
    }
