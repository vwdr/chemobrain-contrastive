"""Regression check for the corrected GEO replicate-pair rescue reference."""
from pathlib import Path
import itertools
import numpy as np
import pandas as pd

R = Path(__file__).resolve().parents[1]
T = R / "analysis" / "corrected_20260920"

x = pd.read_csv(T / "rescue_library_projection.csv")
x = (
    x[x.study.eq("GSE216146") & x.arm.isin(["cisplatin", "cisplatin_GENUS"])]
    .groupby(["sample_id", "arm"], as_index=False)
    .mean(numeric_only=True)
)
lookup = dict(zip(x.sample_id, x.mean_distance))
pairs = [
    ("D20-6409", "D20-6410"),
    ("D21-2750", "D21-2752"),
    ("D21-2751", "D21-2753"),
]
diffs = np.asarray([lookup[r] - lookup[t] for t, r in pairs])
observed = float(diffs.mean())
reference = np.asarray(
    [
        np.mean(diffs * np.asarray(signs))
        for signs in itertools.product([-1, 1], repeat=3)
    ]
)
p = float(np.mean(np.abs(reference) >= abs(observed) - 1e-12))

assert len(reference) == 8
assert abs(observed - 0.00032405988888889) < 1e-12
assert p == 1.0
assert 2 / len(reference) == 0.25

print("PASS rescue replicate-pair sign-flip reference: 8 assignments, P=1.0, two-sided floor=0.25")
