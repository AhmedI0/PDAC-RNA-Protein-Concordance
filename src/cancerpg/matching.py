"""Link GDC RNA-seq aliquots to PDC proteome aliquots of the same tumour/normal sample.

CPTAC-3 tissue was cryopulverised and split for sequencing and mass spectrometry,
so the natural unit is the *sample* (e.g. C3L-01032-01), not the case. The PDC
proteome aliquot of a case is sometimes a pool of two or three tumour pieces
("C3N-01166-01, C3N-01166-02"); the RNA-seq aliquot comes from one of them.

Match levels, from strongest to weakest:
  same_sample      the protein aliquot comes from exactly the RNA sample
  pooled_sample    the protein aliquot is a pool that includes the RNA sample
  same_case        same case and tissue type, but different tissue pieces

Nothing is matched across tissue types (tumour RNA never pairs with normal protein).
"""
from __future__ import annotations

import re

import pandas as pd

TISSUE = {"Primary Tumor": "Tumor", "Solid Tissue Normal": "Normal"}
LEVEL_RANK = {"same_sample": 0, "pooled_sample": 1, "same_case": 2}
_CASE_RE = re.compile(r"^C3[LN]-\d{5}$")


def split_samples(value: str) -> list[str]:
    return [s.strip() for s in str(value).split(",") if s.strip()]


# --------------------------------------------------------------------------- PDC
def pdc_aliquots(bio: pd.DataFrame, matrix_columns: list[str], design: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (aliquots, aliquot_sample_links) for the PDC proteome study.

    `matrix_columns` are the CDAP column labels "<uuid>:<aliquot_submitter_id>".
    Reference pools, QC cell lines and other non-patient channels are flagged.
    """
    in_matrix = {c.split(":", 1)[1]: c for c in matrix_columns}
    plex = design.drop_duplicates("aliquot_submitter_id").set_index("aliquot_submitter_id")
    a = bio.copy()
    a["sample_submitter_ids"] = a["sample_submitter_id"].map(lambda s: ";".join(split_samples(s)))
    a["n_samples_pooled"] = a["sample_submitter_id"].map(lambda s: len(split_samples(s)))
    a["tissue"] = a["sample_type"].map(TISSUE)
    a["is_patient_sample"] = a["case_submitter_id"].astype(str).str.match(_CASE_RE) & a["tissue"].notna()
    a["in_cdap_matrix"] = a["aliquot_submitter_id"].isin(in_matrix)
    a["matrix_column"] = a["aliquot_submitter_id"].map(in_matrix)
    a["plex"] = a["aliquot_submitter_id"].map(plex["plex"])
    a["tmt_channel"] = a["aliquot_submitter_id"].map(plex["channel"])
    links = (
        a[["aliquot_submitter_id", "sample_submitter_id"]]
        .assign(sample_submitter_id=lambda d: d["sample_submitter_id"].map(split_samples))
        .explode("sample_submitter_id")
        .reset_index(drop=True)
    )
    cols = ["aliquot_submitter_id", "aliquot_id", "case_submitter_id", "sample_submitter_ids",
            "n_samples_pooled", "sample_type", "tissue", "disease_type", "is_patient_sample",
            "in_cdap_matrix", "matrix_column", "plex", "tmt_channel"]
    return a[cols].reset_index(drop=True), links


# --------------------------------------------------------------------------- GDC
def gdc_file_table(files: list[dict], kind: str) -> pd.DataFrame:
    """Flatten GDC /files hits into one row per (file, sample, aliquot)."""
    rows = []
    for f in files:
        for c in f.get("cases", []):
            for s in c.get("samples", []):
                for p in s.get("portions", []) or [{}]:
                    for an in p.get("analytes", []) or [{}]:
                        for al in an.get("aliquots", []) or [{}]:
                            rows.append({
                                "file_id": f["file_id"], "kind": kind, "file_name": f["file_name"],
                                "md5sum": f["md5sum"], "file_size": f["file_size"],
                                "workflow": f.get("analysis", {}).get("workflow_type"),
                                "case_id": c["case_id"], "case_submitter_id": c["submitter_id"],
                                "sample_id": s["sample_id"], "sample_submitter_id": s["submitter_id"],
                                "sample_type": s["sample_type"], "tissue": TISSUE.get(s["sample_type"], s.get("tissue_type")),
                                "aliquot_id": al.get("aliquot_id"), "aliquot_submitter_id": al.get("submitter_id"),
                            })
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- matching
def match_rna_protein(rna: pd.DataFrame, prot: pd.DataFrame, links: pd.DataFrame) -> pd.DataFrame:
    """One row per (case, tissue) with the best RNA file for the proteome aliquot.

    rna   : gdc_file_table(kind="rna") restricted to Tumor/Normal samples
    prot  : pdc_aliquots()[0] restricted to patient samples in the CDAP matrix
    links : pdc_aliquots()[1]
    """
    pooled = links.groupby("aliquot_submitter_id")["sample_submitter_id"].apply(set).to_dict()
    out = []
    for _, p in prot.iterrows():
        cand = rna[(rna.case_submitter_id == p.case_submitter_id) & (rna.tissue == p.tissue)]
        best, level = None, None
        for _, r in cand.iterrows():
            samples = pooled.get(p.aliquot_submitter_id, set())
            if r.sample_submitter_id in samples:
                lv = "same_sample" if len(samples) == 1 else "pooled_sample"
            else:
                lv = "same_case"
            # prefer stronger match; tie-break on aliquot id for determinism
            key = (LEVEL_RANK[lv], str(r.aliquot_submitter_id))
            if best is None or key < (LEVEL_RANK[level], str(best.aliquot_submitter_id)):
                best, level = r, lv
        out.append({
            "case_submitter_id": p.case_submitter_id,
            "tissue": p.tissue,
            "protein_aliquot": p.aliquot_submitter_id,
            "protein_samples": p.sample_submitter_ids,
            "protein_column": p.matrix_column,
            "rna_file_id": None if best is None else best.file_id,
            "rna_aliquot": None if best is None else best.aliquot_submitter_id,
            "rna_sample": None if best is None else best.sample_submitter_id,
            "n_rna_candidates": len(cand),
            "match_level": level if best is not None else "protein_only",
        })
    df = pd.DataFrame(out)
    # RNA files whose (case, tissue) has no proteome aliquot
    have = set(zip(df.case_submitter_id, df.tissue))
    extra = (rna[~rna.apply(lambda r: (r.case_submitter_id, r.tissue) in have, axis=1)]
             .sort_values(["case_submitter_id", "tissue", "aliquot_submitter_id"])
             .drop_duplicates(["case_submitter_id", "tissue"]))
    if len(extra):
        df = pd.concat([df, pd.DataFrame({
            "case_submitter_id": extra.case_submitter_id, "tissue": extra.tissue,
            "rna_file_id": extra.file_id, "rna_aliquot": extra.aliquot_submitter_id,
            "rna_sample": extra.sample_submitter_id, "match_level": "rna_only",
        })], ignore_index=True)
    df["sample_key"] = df.case_submitter_id + "_" + df.tissue.str[0]
    return df.sort_values(["tissue", "case_submitter_id"], ascending=[False, True]).reset_index(drop=True)
