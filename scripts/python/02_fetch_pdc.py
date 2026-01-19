from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from cancerpg import pdc  # noqa: E402
from cancerpg.config import load_config, paths  # noqa: E402


def main() -> None:
    cfg = load_config()
    P = paths(cfg)
    out = P.raw / "pdc"
    out.mkdir(parents=True, exist_ok=True)
    sid = cfg["cohort"]["pdc_proteome_study"]

    info = pdc.study_info(sid)
    (out / "study.json").write_text(json.dumps(info, indent=1))
    print(f"{sid}: {info['study_name']} | {info['experiment_type']} | "
          f"{info['cases_count']} cases, {info['aliquots_count']} aliquots")

    bio = pdc.biospecimens(sid)
    bio.to_csv(out / "biospecimens.tsv", sep="\t", index=False)
    print(f"biospecimens: {len(bio)} aliquots; sample types: {bio.sample_type.value_counts().to_dict()}")

    design = pdc.experimental_design(sid)
    design.to_csv(out / "tmt_design.tsv", sep="\t", index=False)
    print(f"TMT design: {design.plex.nunique()} plexes, {len(design)} channel assignments")

    mat = pdc.quant_matrix(sid, cfg["cohort"]["pdc_quant_type"])
    mat.to_csv(out / "proteome_log2ratio.tsv.gz", sep="\t", float_format="%.6g", compression="gzip")
    print(f"proteome matrix: {mat.shape[0]} genes x {mat.shape[1]} aliquots, "
          f"{mat.notna().to_numpy().mean():.1%} values present")


if __name__ == "__main__":
    main()
