"""Run the whole analysis in order (Windows, Linux or macOS).

    python run_pipeline.py                 # every step
    python run_pipeline.py --from 05       # resume from step 05
    python run_pipeline.py --only 09 10    # selected steps
    python run_pipeline.py --list

Rscript is taken from $RSCRIPT, else from PATH. Each step writes its log to
logs/<step>.log. Steps 01-02 download ~1.7 GB from GDC/PDC (skipped when the
files are already present and their md5 sums match).
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
PY = sys.executable
RSCRIPT = os.environ.get("RSCRIPT") or shutil.which("Rscript") or "Rscript"

STEPS = [
    ("00a", "python", "scripts/python/00_fetch_external.py", "pinned PDAC subtype resources (pdacR)"),
    ("00b", "r", "scripts/R/00_export_genesets.R", "MSigDB gene sets, subtype signatures, PurIST model"),
    ("01", "python", "scripts/python/01_fetch_gdc.py", "GDC clinical, RNA-seq, mutations, copy number"),
    ("02", "python", "scripts/python/02_fetch_pdc.py", "PDC proteome (CDAP log2 ratios) and TMT design"),
    ("03", "python", "scripts/python/03_build_database.py", "SQLite database and RNA-protein sample matching"),
    ("04", "python", "scripts/python/04_assemble_matrices.py", "expression, mutation and copy-number matrices"),
    ("05", "r", "scripts/R/05_qc_purity.R", "QC, tissue identity, purity (Fig 2)"),
    ("06", "python", "scripts/python/06_cohort_overview.py", "Table 1 and Fig 1"),
    ("07", "python", "scripts/python/07_rna_protein_concordance.py", "mRNA-protein concordance (Fig 3)"),
    ("08", "python", "scripts/python/08_genomic_cis_effects.py", "copy-number dosage and driver cis effects (Fig 7)"),
    ("09", "r", "scripts/R/09_differential_expression.R", "tumour vs normal, DESeq2 / limma (Fig 4)"),
    ("10", "r", "scripts/R/10_pathways.R", "fgsea, GSVA, pathway-level concordance (Fig 5)"),
    ("11", "r", "scripts/R/11_subtypes.R", "PurIST, Moffitt, consensus clustering (Fig 6)"),
    ("12", "r", "scripts/R/12_survival.R", "survival by subtype, gene and pathway (Fig S1)"),
    ("13", "r", "scripts/R/13_prognostic_models.R", "cross-validated prognostic models (Fig 8)"),
    ("14", "python", "scripts/python/14_export_results_to_db.py", "load results into the database"),
]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--from", dest="start")
    ap.add_argument("--only", nargs="*")
    ap.add_argument("--list", action="store_true")
    a = ap.parse_args()
    if a.list:
        for sid, lang, script, desc in STEPS:
            print(f"{sid:4s} {lang:6s} {script:48s} {desc}")
        return
    steps = STEPS
    if a.only:
        steps = [s for s in STEPS if s[0] in a.only]
    elif a.start:
        steps = STEPS[[s[0] for s in STEPS].index(a.start):]
    (ROOT / "logs").mkdir(exist_ok=True)
    for sid, lang, script, desc in steps:
        cmd = [PY, script] if lang == "python" else [RSCRIPT, script]
        log = ROOT / "logs" / f"{Path(script).stem}.log"
        print(f"[{sid}] {desc} ...", flush=True)
        t0 = time.time()
        with open(log, "w", encoding="utf-8") as fh:
            rc = subprocess.call(cmd, cwd=ROOT, stdout=fh, stderr=subprocess.STDOUT)
        print(f"      {'ok' if rc == 0 else 'FAILED'} in {time.time() - t0:.0f} s  (log: {log.relative_to(ROOT)})", flush=True)
        if rc != 0:
            sys.exit(rc)


if __name__ == "__main__":
    main()
