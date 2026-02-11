"""Small statistical helpers (vectorised where it matters)."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


def bh(p: pd.Series | np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    out = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    if ok.any():
        out[ok] = multipletests(p[ok], method="fdr_bh")[1]
    return out


def rowwise_spearman(a: pd.DataFrame, b: pd.DataFrame, min_n: int = 10) -> pd.DataFrame:
    """Spearman correlation of each row of `a` with the same row of `b` across
    columns (pairwise-complete). Rows and columns must be aligned."""
    assert a.shape == b.shape and (a.index == b.index).all() and (a.columns == b.columns).all()
    A, B = a.to_numpy(float), b.to_numpy(float)
    rho = np.full(len(A), np.nan)
    pval = np.full(len(A), np.nan)
    n = np.zeros(len(A), dtype=int)
    for i in range(len(A)):
        ok = ~(np.isnan(A[i]) | np.isnan(B[i]))
        n[i] = ok.sum()
        if n[i] >= min_n and np.ptp(A[i, ok]) > 0 and np.ptp(B[i, ok]) > 0:
            r = stats.spearmanr(A[i, ok], B[i, ok])
            rho[i], pval[i] = r.statistic, r.pvalue
    return pd.DataFrame({"n": n, "rho": rho, "p": pval}, index=a.index)


def rank_residuals(y: np.ndarray, covars: np.ndarray) -> np.ndarray:
    """Residuals of rank(y) after OLS on ranked covariates (NaNs kept as NaN)."""
    ok = ~np.isnan(y)
    out = np.full_like(y, np.nan, dtype=float)
    ry = stats.rankdata(y[ok])
    X = np.column_stack([np.ones(ok.sum()), np.apply_along_axis(stats.rankdata, 0, covars[ok])])
    beta, *_ = np.linalg.lstsq(X, ry, rcond=None)
    out[ok] = ry - X @ beta
    return out


def partial_spearman(a: pd.DataFrame, b: pd.DataFrame, covars: pd.DataFrame, min_n: int = 20) -> pd.Series:
    """Partial Spearman correlation per row, removing the (ranked) covariates."""
    C = covars.loc[a.columns].to_numpy(float)
    A, B = a.to_numpy(float), b.to_numpy(float)
    out = np.full(len(A), np.nan)
    for i in range(len(A)):
        ok = ~(np.isnan(A[i]) | np.isnan(B[i]))
        if ok.sum() < min_n:
            continue
        ra = rank_residuals(A[i, ok], C[ok])
        rb = rank_residuals(B[i, ok], C[ok])
        out[i] = np.corrcoef(ra, rb)[0, 1]
    return pd.Series(out, index=a.index)


def zscore_rows(df: pd.DataFrame) -> pd.DataFrame:
    return df.sub(df.mean(axis=1), axis=0).div(df.std(axis=1, ddof=1), axis=0)


def geneset_rank_test(values: pd.Series, sets: dict[str, set[str]], min_size: int = 10,
                      max_size: int = 500) -> pd.DataFrame:
    """Two-sided Mann-Whitney test of a per-gene statistic inside vs outside each set."""
    v = values.dropna()
    universe = set(v.index)
    rows = []
    for name, genes in sets.items():
        g = genes & universe
        if not (min_size <= len(g) <= max_size):
            continue
        inside = v[list(g)]
        outside = v.drop(list(g))
        u = stats.mannwhitneyu(inside, outside, alternative="two-sided")
        rows.append({"gene_set": name, "size": len(g), "median_in": inside.median(),
                     "median_out": outside.median(), "auc": u.statistic / (len(inside) * len(outside)),
                     "p": u.pvalue})
    out = pd.DataFrame(rows)
    out["fdr"] = bh(out.p)
    return out.sort_values("median_in", ascending=False).reset_index(drop=True)
