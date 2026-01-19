from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from cancerpg import gdc  # noqa: E402
from cancerpg.config import load_config, paths  # noqa: E402

CASE_EXPAND = [
    "demographic", "diagnoses", "diagnoses.treatments", "diagnoses.pathology_details",
    "follow_ups", "exposures", "project",
    "samples", "samples.portions", "samples.portions.analytes", "samples.portions.analytes.aliquots",
]

FILE_FIELDS = [
    "file_id", "file_name", "md5sum", "file_size", "data_type", "data_format", "access",
    "experimental_strategy", "analysis.workflow_type",
    "cases.case_id", "cases.submitter_id",
    "cases.samples.sample_id", "cases.samples.submitter_id", "cases.samples.sample_type",
    "cases.samples.tissue_type",
    "cases.samples.portions.analytes.aliquots.aliquot_id",
    "cases.samples.portions.analytes.aliquots.submitter_id",
]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--metadata-only", action="store_true")
    args = ap.parse_args()

    cfg = load_config()
    P = paths(cfg)
    out = P.raw / "gdc"
    out.mkdir(parents=True, exist_ok=True)
    g = cfg["gdc"]
    cohort = cfg["cohort"]

    (out / "gdc_status.json").write_text(json.dumps(gdc.status(), indent=1))

    case_filter = gdc.f_and(
        gdc.f_in("project.project_id", cohort["gdc_project"]),
        gdc.f_in("primary_site", cohort["gdc_primary_site"]),
    )
    cases = gdc.query_all("cases", case_filter, expand=CASE_EXPAND, page_size=200)
    (out / "cases.json").write_text(json.dumps(cases, indent=1))
    print(f"cases: {len(cases)}")

    base = [
        gdc.f_in("cases.project.project_id", cohort["gdc_project"]),
        gdc.f_in("cases.primary_site", cohort["gdc_primary_site"]),
        gdc.f_in("access", "open"),
    ]
    kinds = {
        "rna": [gdc.f_in("data_type", "Gene Expression Quantification"),
                gdc.f_in("analysis.workflow_type", g["rna_workflow"])],
        "maf": [gdc.f_in("data_type", "Masked Somatic Mutation"),
                gdc.f_in("analysis.workflow_type", g["maf_workflow"])],
        "cn": [gdc.f_in("data_type", "Gene Level Copy Number"),
               gdc.f_in("analysis.workflow_type", g["cn_workflow"])],
    }
    for kind, extra in kinds.items():
        files = gdc.query_all("files", gdc.f_and(*base, *extra), fields=FILE_FIELDS)
        (out / f"files_{kind}.json").write_text(json.dumps(files, indent=1))
        size_gb = sum(f["file_size"] for f in files) / 1e9
        print(f"{kind}: {len(files)} files, {size_gb:.2f} GB")
        if not args.metadata_only:
            gdc.download_many(files, out / kind, threads=g["download_threads"])


if __name__ == "__main__":
    main()
