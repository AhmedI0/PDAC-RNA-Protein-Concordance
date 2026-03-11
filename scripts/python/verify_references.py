from __future__ import annotations

import html
import re
import sys
import time
import unicodedata
import urllib.parse
from pathlib import Path

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[2]
HEADERS = {"User-Agent": "cancer-proteogenomics-reference-check/0.1"}


def norm(s: str) -> str:
    s = unicodedata.normalize("NFKD", html.unescape(s)).encode("ascii", "ignore").decode()
    return re.sub(r"[^a-z0-9]+", " ", s.lower()).strip()


def fetch(doi: str, retries: int = 5) -> dict:
    url = "https://api.crossref.org/works/" + urllib.parse.quote(doi, safe="")
    for attempt in range(retries):
        r = requests.get(url, headers=HEADERS, timeout=60)
        if r.status_code == 429:
            time.sleep(10 * (attempt + 1))
            continue
        r.raise_for_status()
        return r.json()["message"]
    raise RuntimeError(f"rate limited: {doi}")


def bibtex(key: str, m: dict) -> str:
    authors = " and ".join(
        f"{a.get('family', '')}, {a.get('given', '')}".strip(", ") if "family" in a else a.get("name", "")
        for a in m.get("author", [])
    )
    year = (m.get("issued", {}).get("date-parts") or [[None]])[0][0]
    kind = "article" if m.get("type") == "journal-article" else "misc"
    fields = {
        "author": authors,
        "title": html.unescape(m["title"][0]),
        "journal": html.unescape((m.get("container-title") or [""])[0]) or m.get("publisher", ""),
        "year": year,
        "volume": m.get("volume"),
        "number": m.get("issue"),
        "pages": m.get("page"),
        "doi": m["DOI"],
    }
    body = ",\n".join(f"  {k} = {{{v}}}" for k, v in fields.items() if v)
    return f"@{kind}{{{key},\n{body}\n}}\n"


def main() -> None:
    refs = pd.read_csv(ROOT / "manuscript" / "references.tsv", sep="\t")
    entries, report = [], []
    for _, r in refs.iterrows():
        m = fetch(r.doi)
        ok = norm(m["title"][0]).startswith(norm(r.title_starts))
        year = (m.get("issued", {}).get("date-parts") or [[None]])[0][0]
        report.append({"key": r.key, "doi": r.doi, "title_ok": ok, "year": year,
                       "journal": html.unescape((m.get("container-title") or [""])[0]),
                       "first_author": (m.get("author") or [{}])[0].get("family", "")})
        if ok:
            entries.append(bibtex(r.key, m))
        print(f"{'OK ' if ok else 'BAD'} {r.key:18s} {r.doi}")
        time.sleep(1.2)
    (ROOT / "manuscript" / "references.bib").write_text("\n".join(entries), encoding="utf-8")
    rep = pd.DataFrame(report)
    rep.to_csv(ROOT / "manuscript" / "references_check.tsv", sep="\t", index=False)
    bad = rep[~rep.title_ok]
    print(f"\n{len(rep) - len(bad)}/{len(rep)} DOIs verified")
    sys.exit(1 if len(bad) else 0)


if __name__ == "__main__":
    main()
