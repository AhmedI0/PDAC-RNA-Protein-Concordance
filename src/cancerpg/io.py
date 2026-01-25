from __future__ import annotations

import pandas as pd

from .config import load_config, paths


def _p():
    return paths(load_config())


def analysis_samples() -> pd.DataFrame:
    return pd.read_csv(_p().processed / "analysis_samples.tsv", sep="\t")


def cases() -> pd.DataFrame:
    return pd.read_csv(_p().interim / "cases.tsv", sep="\t")


def matrix(name: str) -> pd.DataFrame:
    return pd.read_csv(_p().processed / name, sep="\t", index_col=0)


def genesets(collections: list[str] | None = None) -> dict[str, set[str]]:
    gs = pd.read_csv(_p().processed / "genesets.tsv.gz", sep="\t")
    if collections:
        gs = gs[gs.collection.isin(collections)]
    return {k: set(v) for k, v in gs.groupby("gs_name")["gene_symbol"]}


def write_table(df: pd.DataFrame, name: str, index: bool = False) -> None:
    df.to_csv(_p().tables / name, sep="\t", index=index, float_format="%.6g")
