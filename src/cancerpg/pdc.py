from __future__ import annotations

import time

import pandas as pd
import requests

API = "https://pdc.cancer.gov/graphql"
_SESSION = requests.Session()
_SESSION.headers.update({"User-Agent": "cancerpg/0.1 (research data client)", "Content-Type": "application/json"})


def gql(query: str, retries: int = 5) -> dict:
    for attempt in range(retries):
        try:
            r = _SESSION.post(API, json={"query": query}, timeout=900)
            r.raise_for_status()
            out = r.json()
            if "errors" in out:
                raise RuntimeError(out["errors"])
            return out["data"]
        except requests.RequestException:
            if attempt == retries - 1:
                raise
            time.sleep(2 ** attempt)
    raise RuntimeError("unreachable")


def study_info(pdc_study_id: str) -> dict:
    q = f"""{{ study(pdc_study_id: "{pdc_study_id}") {{
        pdc_study_id study_submitter_id study_name disease_type primary_site
        analytical_fraction experiment_type acquisition_type cases_count aliquots_count
        study_shortname program_name project_name }} }}"""
    return gql(q)["study"][0]


def biospecimens(pdc_study_id: str) -> pd.DataFrame:
    q = f"""{{ biospecimenPerStudy(pdc_study_id: "{pdc_study_id}") {{
        aliquot_id aliquot_submitter_id sample_id sample_submitter_id case_id case_submitter_id
        sample_type disease_type primary_site pool taxon aliquot_status }} }}"""
    return pd.DataFrame(gql(q)["biospecimenPerStudy"])


def experimental_design(pdc_study_id: str) -> pd.DataFrame:
    """One row per TMT plex with the aliquot in each channel."""
    channels = ["tmt_126", "tmt_127n", "tmt_127c", "tmt_128n", "tmt_128c", "tmt_129n",
                "tmt_129c", "tmt_130n", "tmt_130c", "tmt_131", "tmt_131c"]
    sub = " ".join(f"{c} {{ aliquot_submitter_id }}" for c in channels)
    q = f"""{{ studyExperimentalDesign(pdc_study_id: "{pdc_study_id}") {{
        study_run_metadata_submitter_id plex_dataset_name experiment_number {sub} }} }}"""
    rows = []
    for plex in gql(q)["studyExperimentalDesign"]:
        for c in channels:
            for a in plex.get(c) or []:
                rows.append({
                    "plex": plex["plex_dataset_name"] or plex["study_run_metadata_submitter_id"],
                    "experiment_number": plex["experiment_number"],
                    "channel": c,
                    "aliquot_submitter_id": a["aliquot_submitter_id"],
                })
    return pd.DataFrame(rows)


def quant_matrix(pdc_study_id: str, data_type: str = "log2_ratio") -> pd.DataFrame:
    """Gene x aliquot matrix from the PDC common data analysis pipeline (CDAP).

    Column labels returned by PDC are "<aliquot_uuid>:<aliquot_submitter_id>".
    """
    q = f'{{ quantDataMatrix(pdc_study_id: "{pdc_study_id}", data_type: "{data_type}") }}'
    m = gql(q)["quantDataMatrix"]
    header, body = m[0], m[1:]
    df = pd.DataFrame(body, columns=header).set_index(header[0])
    df.index.name = "gene"
    return df.apply(pd.to_numeric, errors="coerce")
