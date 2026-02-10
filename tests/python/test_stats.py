"""Unit tests for the correlation helpers."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from cancerpg.stats import bh, partial_spearman, rowwise_spearman  # noqa: E402

RNG = np.random.default_rng(1)


def test_rowwise_spearman_matches_scipy_with_missing_values():
    a = pd.DataFrame(RNG.normal(size=(5, 40)))
    b = a + RNG.normal(scale=0.5, size=a.shape)
    b.iloc[0, :5] = np.nan
    r = rowwise_spearman(a, b)
    ok = b.iloc[0].notna()
    expected = stats.spearmanr(a.iloc[0][ok], b.iloc[0][ok]).statistic
    assert np.isclose(r.rho.iloc[0], expected)
    assert r.n.iloc[0] == 35


def test_bh_monotone_and_nan_safe():
    p = np.array([0.01, 0.04, np.nan, 0.03, 0.5])
    q = bh(p)
    assert np.isnan(q[2])
    assert np.all(q[~np.isnan(q)] >= p[~np.isnan(p)])
