"""Flatten GDC case records into analysis-ready clinical tables."""
from __future__ import annotations

import numpy as np
import pandas as pd

_STAGE_GROUP = {
    "Stage IA": "I", "Stage IB": "I", "Stage I": "I",
    "Stage IIA": "II", "Stage IIB": "II", "Stage II": "II",
    "Stage III": "III", "Stage IV": "IV",
}


def _num(x) -> float:
    try:
        return float(x)
    except (TypeError, ValueError):
        return np.nan


def _first(items: list[dict] | None) -> dict:
    return items[0] if items else {}


def case_table(cases: list[dict]) -> pd.DataFrame:
    rows = []
    for c in cases:
        demo = c.get("demographic") or {}
        dx = _first(c.get("diagnoses"))
        path = _first(dx.get("pathology_details"))
        fups = c.get("follow_ups") or []
        fup_days = [_num(f.get("days_to_follow_up")) for f in fups]
        last_fu = np.nanmax([_num(dx.get("days_to_last_follow_up")), *fup_days])
        dead = demo.get("vital_status") == "Dead"
        days_death = _num(demo.get("days_to_death"))
        stage = dx.get("ajcc_pathologic_stage")
        rec_days = [_num(dx.get("days_to_recurrence"))] + [_num(f.get("days_to_recurrence")) for f in fups]
        rec_days = [d for d in rec_days if not np.isnan(d)]
        rec_flag = (str(dx.get("progression_or_recurrence", "")).lower() == "yes"
                    or any(str(f.get("progression_or_recurrence", "")).lower() == "yes" for f in fups))
        rows.append({
            "case_id": c["case_id"],
            "case_submitter_id": c["submitter_id"],
            "disease_type": c.get("disease_type"),
            "primary_diagnosis": dx.get("primary_diagnosis"),
            "morphology": dx.get("morphology"),
            "index_date": c.get("index_date"),
            "age_years": _num(dx.get("age_at_diagnosis")) / 365.25,
            "sex": demo.get("sex_at_birth") or demo.get("gender"),
            "race": demo.get("race"),
            "ethnicity": demo.get("ethnicity"),
            "country_of_birth": demo.get("country_of_birth"),
            "ajcc_stage": stage,
            "stage_group": _STAGE_GROUP.get(stage),
            "ajcc_t": dx.get("ajcc_pathologic_t"),
            "ajcc_n": dx.get("ajcc_pathologic_n"),
            "ajcc_m": dx.get("ajcc_pathologic_m") or dx.get("ajcc_clinical_m"),
            "ajcc_edition": dx.get("ajcc_staging_system_edition"),
            "tumor_grade": dx.get("tumor_grade"),
            "residual_disease": dx.get("residual_disease"),
            "tumor_site": dx.get("site_of_resection_or_biopsy"),
            "tumor_size_cm": _num(path.get("greatest_tumor_dimension")),
            "lymph_nodes_positive": _num(path.get("lymph_nodes_positive")),
            "lymph_nodes_tested": _num(path.get("lymph_nodes_tested")),
            "perineural_invasion": path.get("perineural_invasion_present"),
            "vascular_invasion": path.get("vascular_invasion_present"),
            "vital_status": demo.get("vital_status"),
            "cause_of_death": demo.get("cause_of_death"),
            "days_to_death": days_death,
            "days_to_last_follow_up": last_fu,
            "os_days": days_death if dead else last_fu,
            "os_event": int(dead) if demo.get("vital_status") in ("Dead", "Alive") else np.nan,
            "recurrence": int(rec_flag),
            "days_to_recurrence": min(rec_days) if rec_days else np.nan,
            "n_follow_ups": len(fups),
            "any_neoadjuvant": any(
                str(t.get("treatment_intent_type", "")).lower() == "neoadjuvant" for t in dx.get("treatments") or []
            ),
        })
    df = pd.DataFrame(rows)
    df.loc[df.os_days <= 0, "os_days"] = np.nan
    df.loc[df.os_days.isna(), "os_event"] = np.nan
    return df.sort_values("case_submitter_id").reset_index(drop=True)


def follow_up_table(cases: list[dict]) -> pd.DataFrame:
    keep = ["days_to_follow_up", "timepoint_category", "disease_response", "progression_or_recurrence",
            "days_to_recurrence", "ecog_performance_status", "karnofsky_performance_status"]
    rows = [
        {"case_submitter_id": c["submitter_id"], **{k: f.get(k) for k in keep}}
        for c in cases for f in (c.get("follow_ups") or [])
    ]
    return pd.DataFrame(rows, columns=["case_submitter_id", *keep])


def treatment_table(cases: list[dict]) -> pd.DataFrame:
    keep = ["treatment_type", "treatment_intent_type", "treatment_or_therapy", "treatment_outcome",
            "therapeutic_agents", "days_to_treatment_start"]
    rows = [
        {"case_submitter_id": c["submitter_id"], **{k: t.get(k) for k in keep}}
        for c in cases for dx in (c.get("diagnoses") or []) for t in (dx.get("treatments") or [])
    ]
    return pd.DataFrame(rows, columns=["case_submitter_id", *keep])
