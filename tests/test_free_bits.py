"""
Free-bits / KL-utilization controls of MC-ContrastiveVI on tiny CPU tensors (no training).

Run:  python tests/test_free_bits.py      (or pytest tests/test_free_bits.py)
"""
from __future__ import annotations

import sys
from pathlib import Path

import torch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import src.models.mc_contrastive_vi as mm  # noqa: E402

WEIGHTS = dict(hsic_weight_bg_shared=10., hsic_weight_bg_drug=10., hsic_weight_shared_drug=5.)


def _setup(seed: int = 0):
    torch.manual_seed(seed)
    cfg = mm.MCContrastiveVIConfig(n_genes=12, n_drugs=3, n_batches=1, hidden_dim=8, n_hidden_layers=1,
                                   latent_dim_bg=2, latent_dim_shared=3, latent_dim_drug=2, dropout=0.0)
    model = mm.MCContrastiveVI(cfg)
    x = torch.randn(24, 12)
    d = torch.arange(24) % 3                     # controls, drug 1 and drug 2 all present
    b = torch.zeros(24, dtype=torch.long)
    enc = model(x, d, b)                         # one forward; every loss below reuses it
    return model, enc, x, d


def _reference_loss(model, enc, x, d, kl_weight):
    """The pre-change objective written out with kl_standard_normal."""
    recon_nll = mm.gaussian_nll(x, enc["recon"], model.log_sigma).mean()
    kl_bg = mm.kl_standard_normal(enc["mu_bg"], enc["lv_bg"]).mean()
    t = d > 0
    kl_sh = mm.kl_standard_normal(enc["mu_shared"][t], enc["lv_shared"][t]).mean()
    kl_dr = torch.tensor(0.0)
    for k in range(2):
        m = d == k + 1
        kl_dr = kl_dr + mm.kl_standard_normal(enc["mu_drug_list"][k][m], enc["lv_drug_list"][k][m]).mean()
    return (recon_nll + kl_weight * (kl_bg + kl_sh + kl_dr)
            + 10. * mm.hsic(enc["z_bg"], enc["z_shared"]) + 10. * mm.hsic(enc["z_bg"], enc["z_drug_cat"])
            + 5. * mm.hsic(enc["z_shared"], enc["z_drug_cat"]))


def test_kl_per_dim_sums_to_kl_standard_normal():
    mu, lv = torch.randn(5, 4), torch.randn(5, 4)
    kl = mm.kl_per_dim(mu, lv)
    assert kl.shape == (5, 4)
    assert torch.equal(kl.sum(dim=1), mm.kl_standard_normal(mu, lv))


def test_default_config_is_unchanged_behaviour():
    cfg = mm.MCContrastiveVIConfig(n_genes=4, n_drugs=2)
    assert (cfg.free_bits, cfg.kl_warmup_epochs, cfg.treatment_kl_beta) == (0.0, 0, 1.0)
    assert mm.kl_weight_at_epoch(cfg, 0) == 1.0 and mm.kl_weight_at_epoch(cfg, 50) == 1.0


def test_free_bits_zero_reproduces_existing_loss():
    model, enc, x, d = _setup()
    out = mm.mc_contrastive_loss(model, enc, x, d, 1.0, **WEIGHTS)
    ref = _reference_loss(model, enc, x, d, 1.0)
    assert abs(float(out["loss"]) - float(ref)) < 1e-6
    out_half = mm.mc_contrastive_loss(model, enc, x, d, 0.5, **WEIGHTS)
    assert abs(float(out_half["loss"]) - float(_reference_loss(model, enc, x, d, 0.5))) < 1e-6
    assert float(out["kl_shared_penalized"]) == float(out["kl_shared"])


def test_free_bits_one_raises_mean_shared_kl_per_dim():
    model, enc, x, d = _setup()
    base = mm.mc_contrastive_loss(model, enc, x, d, 1.0, **WEIGHTS)
    model.cfg.free_bits = 1.0
    fb = mm.mc_contrastive_loss(model, enc, x, d, 1.0, **WEIGHTS)
    n_dims = model.cfg.latent_dim_shared
    assert float(base["kl_shared"]) / n_dims < 1.0                      # random init is below the floor
    assert float(fb["kl_shared_penalized"]) / n_dims > float(base["kl_shared_penalized"]) / n_dims
    assert abs(float(fb["kl_shared_penalized"]) / n_dims - 1.0) < 1e-6  # every dim sits at the floor
    assert float(fb["kl_shared"]) == float(base["kl_shared"])           # the logged true KL is untouched
    assert float(fb["kl_bg"]) == float(base["kl_bg"])                   # background is never floored
    assert float(fb["loss"]) > float(base["loss"])


def test_free_bits_gives_no_kl_gradient_below_floor():
    model, enc, x, d = _setup()
    model.cfg.free_bits = 1e3                                           # floor far above any dim's KL
    loss = mm.mc_contrastive_loss(model, enc, x, d, 1.0, **{k: 0.0 for k in WEIGHTS})["loss"]
    loss.backward(retain_graph=True)
    g_floored = model.shared_head.mu.weight.grad.clone()
    model.zero_grad()
    # same objective with the shared/drug KL switched off entirely must give the same shared-head gradient
    model.cfg.free_bits = 0.0
    model.cfg.treatment_kl_beta = 0.0
    loss = mm.mc_contrastive_loss(model, enc, x, d, 1.0, **{k: 0.0 for k in WEIGHTS})["loss"]
    loss.backward()
    assert torch.allclose(g_floored, model.shared_head.mu.weight.grad, atol=1e-7)


def test_treatment_kl_beta_scales_shared_and_drug_only():
    model, enc, x, d = _setup()
    base = mm.mc_contrastive_loss(model, enc, x, d, 1.0, **WEIGHTS)
    model.cfg.treatment_kl_beta = 3.0
    scaled = mm.mc_contrastive_loss(model, enc, x, d, 1.0, **WEIGHTS)
    expect = float(base["loss"]) + 2.0 * (float(base["kl_shared"]) + float(base["kl_drug"]))
    assert abs(float(scaled["loss"]) - expect) < 1e-4


def test_kl_warmup_is_linear_then_flat():
    cfg = mm.MCContrastiveVIConfig(n_genes=4, n_drugs=2, kl_warmup_epochs=10)
    assert [mm.kl_weight_at_epoch(cfg, e) for e in (0, 5, 10, 25)] == [0.0, 0.5, 1.0, 1.0]


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print("ok", name)
