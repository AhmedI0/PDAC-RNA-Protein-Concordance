from __future__ import annotations

import hashlib
import json
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Iterable

import requests

API = "https://api.gdc.cancer.gov/"
_SESSION = requests.Session()
_SESSION.headers.update({"User-Agent": "cancerpg/0.1 (research data client)"})


# --------------------------------------------------------------------------- filters
def f_in(field: str, values: Any) -> dict:
    return {"op": "in", "content": {"field": field, "value": values if isinstance(values, list) else [values]}}


def f_and(*clauses: dict) -> dict:
    return {"op": "and", "content": list(clauses)}


# --------------------------------------------------------------------------- queries
def _post(endpoint: str, payload: dict, retries: int = 5) -> dict:
    for attempt in range(retries):
        try:
            r = _SESSION.post(API + endpoint, json=payload, timeout=300)
            r.raise_for_status()
            return r.json()
        except (requests.RequestException, json.JSONDecodeError):
            if attempt == retries - 1:
                raise
            time.sleep(2 ** attempt)
    raise RuntimeError("unreachable")


def query_all(endpoint: str, filters: dict, fields: Iterable[str] = (), expand: Iterable[str] = (),
              page_size: int = 500) -> list[dict]:
    """Return every hit of a /cases or /files query (follows pagination)."""
    hits: list[dict] = []
    start = 0
    while True:
        payload = {"filters": filters, "size": page_size, "from": start, "format": "JSON"}
        if fields:
            payload["fields"] = ",".join(fields)
        if expand:
            payload["expand"] = ",".join(expand)
        data = _post(endpoint, payload)["data"]
        hits.extend(data["hits"])
        total = data["pagination"]["total"]
        start += page_size
        if start >= total:
            return hits


def status() -> dict:
    r = _SESSION.get(API + "status", timeout=60)
    r.raise_for_status()
    return r.json()


# --------------------------------------------------------------------------- downloads
def md5(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def download_file(file_id: str, dest: Path, expected_md5: str | None = None, retries: int = 6) -> Path:
    """Download one open-access file to `dest`; skip if present with matching md5."""
    if dest.exists() and (expected_md5 is None or md5(dest) == expected_md5):
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    for attempt in range(retries):
        try:
            with _SESSION.get(API + "data/" + file_id, stream=True, timeout=600) as r:
                r.raise_for_status()
                with open(tmp, "wb") as fh:
                    for block in r.iter_content(1 << 20):
                        fh.write(block)
            if expected_md5 and md5(tmp) != expected_md5:
                raise IOError(f"md5 mismatch for {file_id}")
            tmp.replace(dest)
            return dest
        except (requests.RequestException, IOError):
            if attempt == retries - 1:
                raise
            time.sleep(2 ** attempt)
    raise RuntimeError("unreachable")


def download_many(records: list[dict], out_dir: Path, threads: int = 6) -> list[Path]:
    """Download records with keys file_id, file_name, md5sum into out_dir/<file_id>/<file_name>."""
    out: list[Path] = []
    with ThreadPoolExecutor(max_workers=threads) as pool:
        futs = {
            pool.submit(download_file, r["file_id"], out_dir / r["file_id"] / r["file_name"], r.get("md5sum")): r
            for r in records
        }
        for i, fut in enumerate(as_completed(futs), 1):
            out.append(fut.result())
            if i % 50 == 0 or i == len(futs):
                print(f"  downloaded {i}/{len(futs)}", flush=True)
    return out
