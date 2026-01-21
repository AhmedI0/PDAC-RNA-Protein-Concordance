"""Step 03 - Build the SQLite database: clinical tables, GDC files, PDC aliquots
and the RNA <-> protein sample matches with inclusion flags.

Inclusion (per case and tissue):
  * case diagnosed as PDAC in GDC (primary_diagnosis "Infiltrating duct carcinoma, NOS")
    and labelled "Pancreatic Ductal Adenocarcinoma" in the PDC proteome study
  * a patient proteome aliquot present in the CDAP matrix
  * an RNA-seq file from the same case and tissue type (match level recorded)

    python scripts/python/03_build_database.py
"""
from __future__ import annotations

import json
import sqlite3
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from cancerpg import clinical, matching  # noqa: E402
from cancerpg.config import load_config, paths  # noqa: E402

PDAC_GDC = "Infiltrating duct carcinoma, NOS"
PDAC_PDC = "Pancreatic Ductal Adenocarcinoma"


def main() -> None:
    cfg = load_config()
    P = paths(cfg)
    gdir, pdir = P.raw / "gdc", P.raw / "pdc"

    cases_json = json.loads((gdir / "cases.json").read_text())
    cases = clinical.case_table(cases_json)
    fups = clinical.follow_up_table(cases_json)
    trts = clinical.treatment_table(cases_json)

    files = pd.concat(
        [matching.gdc_file_table(json.loads((gdir / f"files_{k}.json").read_text()), k) for k in ("rna", "maf", "cn")],
        ignore_index=True,
    )

    bio = pd.read_csv(pdir / "biospecimens.tsv", sep="\t")
    design = pd.read_csv(pdir / "tmt_design.tsv", sep="\t")
    mat_cols = pd.read_csv(pdir / "proteome_log2ratio.tsv.gz", sep="\t", index_col=0, nrows=0).columns.tolist()
    aliq, links = matching.pdc_aliquots(bio, mat_cols, design)

    rna = files[(files.kind == "rna") & files.tissue.isin(["Tumor", "Normal"])]
    prot = aliq[aliq.is_patient_sample & aliq.in_cdap_matrix]
    m = matching.match_rna_protein(rna, prot, links)

    pdac_cases = set(cases.loc[cases.primary_diagnosis == PDAC_GDC, "case_submitter_id"])
    pdc_disease = prot.set_index("aliquot_submitter_id")["disease_type"]

    def reason(r) -> str | None:
        if r.case_submitter_id not in pdac_cases:
            return "not PDAC in GDC"
        if pd.notna(r.protein_aliquot) and pdc_disease.get(r.protein_aliquot) != PDAC_PDC:
            return "not PDAC in PDC"
        if r.match_level == "rna_only":
            return "no proteome aliquot"
        if r.match_level == "protein_only":
            return "no RNA-seq"
        return None

    m["exclusion_reason"] = m.apply(reason, axis=1)
    m["included"] = m["exclusion_reason"].isna().astype(int)

    # ------------------------------------------------------------------ write database
    P.db.parent.mkdir(parents=True, exist_ok=True)
    if P.db.exists():
        P.db.unlink()
    con = sqlite3.connect(P.db)
    con.executescript((P.root / "db" / "schema.sql").read_text())
    known = set(cases.case_submitter_id)
    cases.assign(any_neoadjuvant=cases.any_neoadjuvant.astype(int)).to_sql("cases", con, if_exists="append", index=False)
    fups.to_sql("follow_ups", con, if_exists="append", index=False)
    trts.to_sql("treatments", con, if_exists="append", index=False)
    (files[files.case_submitter_id.isin(known)]
     .to_sql("gdc_files", con, if_exists="append", index=False))
    aliq.assign(is_patient_sample=aliq.is_patient_sample.astype(int),
                in_cdap_matrix=aliq.in_cdap_matrix.astype(int)).to_sql("pdc_aliquots", con, if_exists="append", index=False)
    links.to_sql("pdc_aliquot_samples", con, if_exists="append", index=False)
    m[m.case_submitter_id.isin(known)].to_sql("sample_matches", con, if_exists="append", index=False)
    con.commit()

    # ------------------------------------------------------------------ interim tables
    P.interim.mkdir(parents=True, exist_ok=True)
    cases.to_csv(P.interim / "cases.tsv", sep="\t", index=False)
    m.to_csv(P.interim / "sample_matches.tsv", sep="\t", index=False)
    files.to_csv(P.interim / "gdc_files.tsv", sep="\t", index=False)
    aliq.to_csv(P.interim / "pdc_aliquots.tsv", sep="\t", index=False)

    print(f"cases: {len(cases)} ({len(pdac_cases)} PDAC); follow-ups {len(fups)}; treatments {len(trts)}")
    print(f"GDC file rows: {files.groupby('kind').file_id.nunique().to_dict()}")
    print(f"PDC aliquots: {len(aliq)}; patient samples in CDAP matrix: {len(prot)}")
    print("\nmatch levels (all):")
    print(m.groupby(["tissue", "match_level"]).size().to_string())
    print("\nincluded:")
    print(m[m.included == 1].groupby(["tissue", "match_level"]).size().to_string())
    print("\nexclusions:")
    print(m[m.included == 0].groupby(["tissue", "exclusion_reason"]).size().to_string())
    for q in ["SELECT COUNT(*) FROM v_analysis_samples",
              "SELECT SUM(rna_tumor AND protein_tumor AND has_os) FROM v_case_availability"]:
        print(q, "->", con.execute(q).fetchone()[0])
    con.close()


if __name__ == "__main__":
    main()
