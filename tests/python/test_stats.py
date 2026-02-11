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


def test_partial_spearman_removes_shared_confounder():
    n = 300
    z = RNG.normal(size=n)                       # e.g. stromal content
    x = z + RNG.normal(scale=0.3, size=n)        # mRNA driven by composition
    y = z + RNG.normal(scale=0.3, size=n)        # protein driven by composition
    a = pd.DataFrame([x], columns=range(n))
    b = pd.DataFrame([y], columns=range(n))
    raw = rowwise_spearman(a, b).rho.iloc[0]
    part = partial_spearman(a, b, pd.DataFrame({"z": z}, index=range(n))).iloc[0]
    assert raw > 0.8
    assert abs(part) < 0.15


def test_bh_monotone_and_nan_safe():
    p = np.array([0.01, 0.04, np.nan, 0.03, 0.5])
    q = bh(p)
    assert np.isnan(q[2])
    assert np.all(q[~np.isnan(q)] >= p[~np.isnan(p)])
