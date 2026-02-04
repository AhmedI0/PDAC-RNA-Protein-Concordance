from __future__ import annotations

import hashlib
import sys
from pathlib import Path

import requests

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from cancerpg.config import paths  # noqa: E402

PDACR_SHA = "5b044b2f3258bb2ac9d7d4300b31b6b65d9520af"
FILES = ["gene_lists.rds", "Moffitt_classifier_2019.RData"]


def main() -> None:
    out = paths().raw / "external" / f"pdacR_{PDACR_SHA}"
    out.mkdir(parents=True, exist_ok=True)
    for f in FILES:
        url = f"https://raw.githubusercontent.com/rmoffitt/pdacR/{PDACR_SHA}/data/{f}"
        r = requests.get(url, timeout=120)
        r.raise_for_status()
        (out / f).write_bytes(r.content)
        print(f"{f}: {len(r.content)} bytes, sha256 {hashlib.sha256(r.content).hexdigest()[:16]}")


if __name__ == "__main__":
    main()
