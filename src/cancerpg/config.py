from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]


@dataclass(frozen=True)
class Paths:
    root: Path
    raw: Path
    interim: Path
    processed: Path
    db: Path
    tables: Path
    figures: Path

    def ensure(self) -> "Paths":
        for p in (self.raw, self.interim, self.processed, self.tables, self.figures, self.db.parent):
            p.mkdir(parents=True, exist_ok=True)
        return self


@lru_cache(maxsize=1)
def load_config(path: str | Path | None = None) -> dict:
    path = Path(path) if path else ROOT / "config" / "config.yaml"
    with open(path, encoding="utf-8") as fh:
        return yaml.safe_load(fh)


def paths(cfg: dict | None = None) -> Paths:
    cfg = cfg or load_config()
    p = cfg["paths"]
    return Paths(
        root=ROOT,
        raw=ROOT / p["raw"],
        interim=ROOT / p["interim"],
        processed=ROOT / p["processed"],
        db=ROOT / p["db"],
        tables=ROOT / p["tables"],
        figures=ROOT / p["figures"],
    ).ensure()
