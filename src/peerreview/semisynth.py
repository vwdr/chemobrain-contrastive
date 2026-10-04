"""Semi-synthetic positive control: binomial thinning of raw counts (Gerard 2020, BMC Bioinformatics).

For a gene with target |log2 fold change| delta, every count is replaced by Binomial(count, p) with
p = 2**(-delta). An *up* program thins the gene in every cell of the relevant study except the
responders; a *down* program thins the gene in the responders. Either way responders end up with
expected expression 2**delta (up) or 2**(-delta) (down) relative to all other cells of the study.
"""
from __future__ import annotations

import numpy as np


def thin_counts(counts: np.ndarray, gene_cols, cell_mask: np.ndarray, delta: float, rng) -> np.ndarray:
    """Return a copy of ``counts`` with ``gene_cols`` thinned (keep prob. 2**-delta) in ``cell_mask`` rows."""
    out = counts.copy()
    p = 2.0 ** (-float(delta))
    rows = np.where(cell_mask)[0]
    cols = np.asarray(gene_cols, dtype=int)
    if len(rows) and len(cols):
        sub = out[np.ix_(rows, cols)]
        out[np.ix_(rows, cols)] = rng.binomial(np.round(sub).astype(np.int64), p).astype(out.dtype)
    return out


def inject(counts: np.ndarray, programs, rng) -> np.ndarray:
    """Apply a list of programs to raw integer counts.

    Each program is a dict with keys: ``up_cols``, ``down_cols`` (column indices into ``counts``),
    ``study_mask`` (cells of the relevant study or studies), ``responders`` (boolean mask, subset of
    ``study_mask``) and ``delta``. Up genes are thinned in ``study_mask & ~responders``; down genes in
    ``responders``. Columns of different programs must be disjoint.
    """
    out = counts.copy()
    for prog in programs:
        if prog['delta'] <= 0:
            continue
        resp = prog['responders'] & prog['study_mask']
        out = thin_counts(out, prog['up_cols'], prog['study_mask'] & ~resp, prog['delta'], rng)
        out = thin_counts(out, prog['down_cols'], resp, prog['delta'], rng)
    return out


def log_normalize(counts_model: np.ndarray, library_full: np.ndarray, target_sum: float = 1e4) -> np.ndarray:
    """log1p(counts / library * 1e4) with the library computed over all genes (canonical normalization)."""
    return np.log1p(counts_model.astype(np.float64) / library_full[:, None].astype(np.float64) * target_sum).astype(np.float32)
