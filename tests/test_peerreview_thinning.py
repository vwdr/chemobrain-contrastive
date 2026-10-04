"""Unit test for the binomial-thinning injection used in the semi-synthetic positive control."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.peerreview.semisynth import inject  # noqa: E402


def _realized_log2fc(x, resp, others, cols):
    return np.log2(x[np.ix_(resp, cols)].mean(0) / x[np.ix_(others, cols)].mean(0))


def test_thinning_realizes_delta_and_leaves_others_untouched():
    rng = np.random.default_rng(1)
    n_cells, n_genes = 100000, 30
    lam = rng.uniform(5, 50, size=n_genes)
    counts = rng.poisson(lam, size=(n_cells, n_genes)).astype(np.float32)
    study = np.zeros(n_cells, dtype=bool); study[: n_cells // 2] = True          # target study = first half
    treated = np.zeros(n_cells, dtype=bool); treated[n_cells // 4: n_cells // 2] = True
    responders = treated & (rng.random(n_cells) < 0.3)                         # 30% responder fraction
    up, down, untouched = [0, 1, 2, 3, 4], [5, 6, 7, 8, 9], list(range(10, n_genes))
    for delta in (0.25, 0.5, 1.0, 2.0):
        out = inject(counts, [dict(up_cols=up, down_cols=down, study_mask=study, responders=responders, delta=delta)],
                     np.random.default_rng(2))
        others = np.where(study & ~responders)[0]; resp = np.where(responders)[0]
        fc_up = _realized_log2fc(out, resp, others, up)
        fc_down = _realized_log2fc(out, resp, others, down)
        assert np.all(np.abs(fc_up - delta) < 0.05), (delta, fc_up)
        assert np.all(np.abs(fc_down + delta) < 0.05), (delta, fc_down)
        # non-target genes untouched everywhere
        assert np.array_equal(out[:, untouched], counts[:, untouched])
        # non-target study untouched for every gene
        assert np.array_equal(out[~study], counts[~study])
        # thinning never increases a count and keeps integers
        assert np.all(out <= counts) and np.allclose(out, np.round(out))
        # up genes untouched in responders; down genes untouched in non-responders
        assert np.array_equal(out[np.ix_(resp, up)], counts[np.ix_(resp, up)])
        assert np.array_equal(out[np.ix_(others, down)], counts[np.ix_(others, down)])


def test_zero_delta_is_identity():
    counts = np.random.default_rng(3).poisson(10, size=(100, 5)).astype(np.float32)
    mask = np.ones(100, dtype=bool)
    out = inject(counts, [dict(up_cols=[0], down_cols=[1], study_mask=mask, responders=mask, delta=0.0)],
                 np.random.default_rng(4))
    assert np.array_equal(out, counts)


if __name__ == '__main__':
    test_thinning_realizes_delta_and_leaves_others_untouched()
    test_zero_delta_is_identity()
    print('PASS thinning unit tests')
