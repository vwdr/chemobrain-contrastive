"""Model-agnostic integrated gradients with the adaptive Gauss-Legendre rule of scripts/09_interpretation.py.

``fun`` maps an (n, G) input tensor to an (n,) output (e.g. the squared norm of a posterior mean). The quadrature
order is increased over 256 -> 1024 -> 4096 until the median absolute completeness error divided by the mean
absolute output difference is < 0.01 (09_interpretation.py rule).
"""
from __future__ import annotations

import numpy as np
import torch


def integrated_gradients(fun, inputs: np.ndarray, baseline: np.ndarray):
    inp = torch.tensor(np.asarray(inputs, dtype=np.float32)); ref = torch.tensor(np.asarray(baseline, dtype=np.float32))[None, :]
    delta = inp - ref
    with torch.no_grad():
        diff = (fun(inp) - fun(ref)).numpy()
    for order in (256, 1024, 4096):
        nodes, weights = np.polynomial.legendre.leggauss(order); acc = torch.zeros_like(inp)
        for node, weight in zip(nodes, weights):
            xx = (ref + delta * float((node + 1) / 2)).detach().requires_grad_(True)
            grad = torch.autograd.grad(fun(xx).sum(), xx)[0]; acc += grad * float(weight / 2)
        ig = (delta * acc).detach().numpy(); err = np.abs(ig.sum(1) - diff)
        scaled = float(np.median(err) / (np.mean(np.abs(diff)) + 1e-8))
        if scaled < .01:
            break
    return ig, {'quadrature_points': order, 'median_scaled_error': scaled, 'median_absolute_error': float(np.median(err)),
                'mean_abs_output_difference': float(np.mean(np.abs(diff)))}


def recovery(mean_abs_ig: np.ndarray, genes, target_sets: dict, top: int = 50):
    """Top-``top`` hits against the union of ``target_sets`` (name -> gene set) and per-set breakdown."""
    genes = np.asarray(genes); topg = set(genes[np.argsort(mean_abs_ig)[-top:]].tolist())
    union = set().union(*target_sets.values())
    out = {'top50_hits_union': len(topg & union), 'union_size': len(union), 'precision_union': len(topg & union) / top,
           'recall_union': len(topg & union) / len(union), 'random_expectation_hits': top * len(union) / len(genes)}
    for k, v in target_sets.items():
        out[f'top50_hits_{k}'] = len(topg & set(v))
    return out, sorted(topg)
