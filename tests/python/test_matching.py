"""Unit tests for the RNA <-> protein sample matching rules."""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from cancerpg import matching  # noqa: E402


def _rna(rows):
    return pd.DataFrame(rows, columns=["file_id", "case_submitter_id", "tissue", "sample_submitter_id",
                                       "aliquot_submitter_id"])


def _prot(rows):
    df = pd.DataFrame(rows, columns=["aliquot_submitter_id", "case_submitter_id", "tissue", "sample_submitter_ids"])
    df["matrix_column"] = "uuid:" + df.aliquot_submitter_id
    return df


def _links(prot):
    return (prot[["aliquot_submitter_id", "sample_submitter_ids"]]
            .assign(sample_submitter_id=lambda d: d.sample_submitter_ids.str.split(";"))
            .explode("sample_submitter_id")[["aliquot_submitter_id", "sample_submitter_id"]])


def test_split_samples_handles_pools_and_spaces():
    assert matching.split_samples("C3N-01166-01, C3N-01166-02") == ["C3N-01166-01", "C3N-01166-02"]
    assert matching.split_samples("C3L-00017-01") == ["C3L-00017-01"]


def test_match_levels_and_preference():
    rna = _rna([
        ("f1", "C3L-00001", "Tumor", "C3L-00001-01", "CPT01"),
        ("f2", "C3L-00002", "Tumor", "C3L-00002-02", "CPT02"),
        ("f3", "C3L-00003", "Tumor", "C3L-00003-03", "CPT03"),
        ("f4", "C3L-00004", "Tumor", "C3L-00004-01", "CPT04a"),   # two RNA files, one same-sample
        ("f5", "C3L-00004", "Tumor", "C3L-00004-02", "CPT04b"),
        ("f6", "C3L-00005", "Normal", "C3L-00005-06", "CPT05"),
    ])
    prot = _prot([
        ("P1", "C3L-00001", "Tumor", "C3L-00001-01"),
        ("P2", "C3L-00002", "Tumor", "C3L-00002-01;C3L-00002-02"),
        ("P3", "C3L-00003", "Tumor", "C3L-00003-01"),
        ("P4", "C3L-00004", "Tumor", "C3L-00004-02"),
        ("P6", "C3L-00006", "Tumor", "C3L-00006-01"),
    ])
    m = matching.match_rna_protein(rna, prot, _links(prot)).set_index("sample_key")
    assert m.loc["C3L-00001_T", "match_level"] == "same_sample"
    assert m.loc["C3L-00002_T", "match_level"] == "pooled_sample"
    assert m.loc["C3L-00003_T", "match_level"] == "same_case"
    assert m.loc["C3L-00004_T", "rna_file_id"] == "f5"           # the same-sample file wins
    assert m.loc["C3L-00006_T", "match_level"] == "protein_only"
    assert m.loc["C3L-00005_N", "match_level"] == "rna_only"


def test_no_cross_tissue_matching():
    rna = _rna([("f1", "C3L-00001", "Normal", "C3L-00001-06", "CPT01")])
    prot = _prot([("P1", "C3L-00001", "Tumor", "C3L-00001-01")])
    m = matching.match_rna_protein(rna, prot, _links(prot)).set_index("sample_key")
    assert m.loc["C3L-00001_T", "match_level"] == "protein_only"
    assert m.loc["C3L-00001_N", "match_level"] == "rna_only"
