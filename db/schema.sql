-- PDAC-RNA-Protein-Concordance SQLite database
-- Built by scripts/python/03_build_database.py; analysis results are added by
-- scripts/python/16_export_results_to_db.py.

PRAGMA foreign_keys = ON;

-- ----------------------------------------------------------------- clinical
CREATE TABLE cases (
    case_submitter_id      TEXT PRIMARY KEY,
    case_id                TEXT UNIQUE NOT NULL,      -- GDC case UUID
    disease_type           TEXT,
    primary_diagnosis      TEXT,
    morphology             TEXT,
    index_date             TEXT,
    age_years              REAL,
    sex                    TEXT,
    race                   TEXT,
    ethnicity              TEXT,
    country_of_birth       TEXT,
    ajcc_stage             TEXT,
    stage_group            TEXT,
    ajcc_t                 TEXT,
    ajcc_n                 TEXT,
    ajcc_m                 TEXT,
    ajcc_edition           TEXT,
    tumor_grade            TEXT,
    residual_disease       TEXT,
    tumor_site             TEXT,
    tumor_size_cm          REAL,
    lymph_nodes_positive   REAL,
    lymph_nodes_tested     REAL,
    perineural_invasion    TEXT,
    vascular_invasion      TEXT,
    vital_status           TEXT,
    cause_of_death         TEXT,
    days_to_death          REAL,
    days_to_last_follow_up REAL,
    os_days                REAL,
    os_event               INTEGER,
    recurrence             INTEGER,
    days_to_recurrence     REAL,
    n_follow_ups           INTEGER,
    any_neoadjuvant        INTEGER
);

CREATE TABLE follow_ups (
    case_submitter_id      TEXT NOT NULL REFERENCES cases(case_submitter_id),
    days_to_follow_up      REAL,
    timepoint_category     TEXT,
    disease_response       TEXT,
    progression_or_recurrence TEXT,
    days_to_recurrence     REAL,
    ecog_performance_status TEXT,
    karnofsky_performance_status TEXT
);

CREATE TABLE treatments (
    case_submitter_id      TEXT NOT NULL REFERENCES cases(case_submitter_id),
    treatment_type         TEXT,
    treatment_intent_type  TEXT,
    treatment_or_therapy   TEXT,
    treatment_outcome      TEXT,
    therapeutic_agents     TEXT,
    days_to_treatment_start REAL
);

-- ----------------------------------------------------------------- biospecimens and files
CREATE TABLE gdc_files (
    file_id                TEXT NOT NULL,             -- GDC file UUID
    kind                   TEXT NOT NULL CHECK (kind IN ('rna', 'maf', 'cn')),
    file_name              TEXT NOT NULL,
    md5sum                 TEXT NOT NULL,
    file_size              INTEGER,
    workflow               TEXT,
    case_submitter_id      TEXT NOT NULL REFERENCES cases(case_submitter_id),
    sample_submitter_id    TEXT NOT NULL,
    sample_type            TEXT,
    tissue                 TEXT,
    aliquot_submitter_id   TEXT,
    aliquot_id             TEXT,
    PRIMARY KEY (file_id, sample_submitter_id, aliquot_submitter_id)
);

CREATE TABLE pdc_aliquots (
    aliquot_submitter_id   TEXT PRIMARY KEY,
    aliquot_id             TEXT,                      -- PDC aliquot UUID
    case_submitter_id      TEXT,
    sample_submitter_ids   TEXT,                      -- ';'-separated when pooled
    n_samples_pooled       INTEGER,
    sample_type            TEXT,
    tissue                 TEXT,
    disease_type           TEXT,
    is_patient_sample      INTEGER,
    in_cdap_matrix         INTEGER,
    matrix_column          TEXT,
    plex                   TEXT,
    tmt_channel            TEXT
);

CREATE TABLE pdc_aliquot_samples (
    aliquot_submitter_id   TEXT NOT NULL REFERENCES pdc_aliquots(aliquot_submitter_id),
    sample_submitter_id    TEXT NOT NULL,
    PRIMARY KEY (aliquot_submitter_id, sample_submitter_id)
);

-- One row per (case, tissue): the RNA file chosen for the proteome aliquot.
CREATE TABLE sample_matches (
    sample_key             TEXT PRIMARY KEY,          -- <case>_T or <case>_N
    case_submitter_id      TEXT NOT NULL REFERENCES cases(case_submitter_id),
    tissue                 TEXT NOT NULL CHECK (tissue IN ('Tumor', 'Normal')),
    protein_aliquot        TEXT REFERENCES pdc_aliquots(aliquot_submitter_id),
    protein_samples        TEXT,
    protein_column         TEXT,
    rna_file_id            TEXT,
    rna_aliquot            TEXT,
    rna_sample             TEXT,
    n_rna_candidates       INTEGER,
    match_level            TEXT NOT NULL CHECK (match_level IN
                           ('same_sample', 'pooled_sample', 'same_case', 'protein_only', 'rna_only')),
    included               INTEGER NOT NULL,
    exclusion_reason       TEXT
);

-- ----------------------------------------------------------------- convenience views
CREATE VIEW v_analysis_samples AS
SELECT m.*, c.primary_diagnosis, c.stage_group, c.os_days, c.os_event
FROM sample_matches m JOIN cases c USING (case_submitter_id)
WHERE m.included = 1;

CREATE VIEW v_case_availability AS
SELECT c.case_submitter_id,
       c.primary_diagnosis,
       MAX(CASE WHEN m.tissue = 'Tumor'  AND m.rna_file_id IS NOT NULL THEN 1 ELSE 0 END) AS rna_tumor,
       MAX(CASE WHEN m.tissue = 'Tumor'  AND m.protein_aliquot IS NOT NULL THEN 1 ELSE 0 END) AS protein_tumor,
       MAX(CASE WHEN m.tissue = 'Normal' AND m.rna_file_id IS NOT NULL THEN 1 ELSE 0 END) AS rna_normal,
       MAX(CASE WHEN m.tissue = 'Normal' AND m.protein_aliquot IS NOT NULL THEN 1 ELSE 0 END) AS protein_normal,
       (SELECT COUNT(*) FROM gdc_files f WHERE f.case_submitter_id = c.case_submitter_id AND f.kind = 'maf' AND f.tissue = 'Tumor') > 0 AS has_maf,
       (SELECT COUNT(*) FROM gdc_files f WHERE f.case_submitter_id = c.case_submitter_id AND f.kind = 'cn'  AND f.tissue = 'Tumor') > 0 AS has_cn,
       c.os_days IS NOT NULL AS has_os
FROM cases c LEFT JOIN sample_matches m USING (case_submitter_id)
GROUP BY c.case_submitter_id;
