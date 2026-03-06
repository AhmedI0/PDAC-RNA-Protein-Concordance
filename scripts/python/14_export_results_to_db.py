from __future__ import annotations

import sqlite3
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from cancerpg.config import load_config, paths  # noqa: E402

DESCRIPTIONS = {
    "table1_cohort": ("06", "Clinical characteristics of the matched tumour cohort (Table 1)"),
    "cohort_data_availability": ("06", "PDAC cases with each data type after QC"),
    "qc_rna_samples": ("05", "RNA library size, detected genes and size factors"),
    "qc_protein_samples": ("05", "Proteins quantified per sample, TMT plex"),
    "qc_protein_pc_variance": ("05", "Protein PCs: variance explained by tissue and by TMT plex"),
    "qc_tissue_identity": ("05", "Intestinal vs acinar marker scores; duodenal-mucosa flag"),
    "purity_estimates": ("05", "ESTIMATE stromal/immune scores and KRAS-VAF cellularity"),
    "genesets_pdac_subtypes": ("00", "PDAC-related gene sets found in MSigDB C2:CGP"),
    "concordance_gene": ("07", "Gene-wise mRNA-protein Spearman correlations (tumour, normal, pooled, partial)"),
    "concordance_sample": ("07", "Per-tumour agreement of gene-centred mRNA and protein profiles"),
    "concordance_summary": ("07", "Summary statistics of mRNA-protein concordance"),
    "concordance_determinants": ("07", "Correlation of gene-wise rho with expression level and variability"),
    "concordance_hallmark": ("07", "Hallmark gene sets: median gene-wise rho, rank test"),
    "concordance_hallmark_partial": ("07", "Hallmark gene sets, composition-adjusted rho"),
    "concordance_reactome": ("07", "Reactome gene sets: median gene-wise rho, rank test"),
    "concordance_complexes": ("07", "Protein complexes (GO CC): median gene-wise rho"),
    "cis_cn_dosage": ("08", "Copy number vs own mRNA / protein, attenuation per gene"),
    "cis_cn_summary": ("08", "Summary of dosage transmission and buffering"),
    "cis_driver_effects": ("08", "mRNA and protein of KRAS/TP53/CDKN2A/SMAD4 by alteration class"),
    "driver_alteration_status": ("08", "Alteration class per tumour for the four driver genes"),
    "de_rna_tumour_vs_normal": ("09", "DESeq2 tumour vs normal (apeglm-shrunken log2FC)"),
    "de_protein_tumour_vs_normal": ("09", "limma tumour vs normal (plex covariate, patient block)"),
    "de_rna_protein_categories": ("09", "Per-gene concordance category of tumour-vs-normal changes"),
    "de_rna_protein_categories_matched": ("09", "Same, both layers restricted to identical matched samples"),
    "de_summary": ("09", "Differential expression summary"),
    "de_category_counts": ("09", "Genes per concordance category"),
    "gsea_tumour_vs_normal": ("10", "fgsea results (Hallmark, Reactome) per layer"),
    "gsea_rna_vs_protein": ("10", "NES of each pathway in mRNA and protein, agreement class"),
    "gsea_nes_agreement": ("10", "Spearman of NES between layers (all / without complex-dominated sets)"),
    "pathway_level_concordance": ("10", "GSVA Hallmark score correlation mRNA vs protein"),
    "subtype_signature_coverage": ("11", "Moffitt genes measured in each layer"),
    "subtype_calls": ("11", "PurIST, Moffitt scores and consensus calls, unsupervised clusters per tumour"),
    "subtype_agreement": ("11", "Cohen's kappa between subtype calls"),
    "subtype_score_auc": ("11", "AUC of each layer's Moffitt score against PurIST"),
    "consensus_clustering_pac": ("11", "Proportion of ambiguous clustering by k and layer"),
    "consensus_clustering_ari": ("11", "Adjusted Rand index between layer clusterings by k"),
    "cluster_hallmark_protein": ("11", "Hallmark GSVA scores by protein cluster"),
    "cluster_hallmark_rna": ("11", "Hallmark GSVA scores by mRNA cluster"),
    "survival_cox_models": ("12", "Adjusted Cox models: clinical, subtypes, per-tumour concordance"),
    "survival_ph_checks": ("12", "Proportional-hazards tests (Schoenfeld residuals)"),
    "survival_logrank": ("12", "Log-rank tests by subtype / cluster"),
    "survival_gene_cox": ("12", "Stage-adjusted Cox per gene, mRNA and protein side by side"),
    "survival_gene_cox_long": ("12", "Stage-adjusted Cox per gene and layer (long format)"),
    "survival_gene_cox_summary": ("12", "Counts of prognostic genes per layer"),
    "survival_hr_agreement_by_concordance": ("12", "Agreement of gene log HRs by RNA-protein concordance bin"),
    "survival_pathway_cox": ("12", "Adjusted Cox per Hallmark GSVA score and layer"),
    "prognostic_cv_cindex_by_repeat": ("13", "Cross-validated C-index per model and repeat"),
    "prognostic_cv_summary": ("13", "Mean C-index and gain over the clinical model"),
    "prognostic_signature_genes": ("13", "Stable signature genes (selection frequency) and final coefficients"),
    "prognostic_cv_risk_groups": ("13", "Out-of-fold risk group per tumour"),
    "prognostic_cv_risk_cox": ("13", "Adjusted Cox of the out-of-fold risk groups"),
}


def main() -> None:
    P = paths(load_config())
    con = sqlite3.connect(P.db)
    qc = pd.read_csv(P.processed / "analysis_samples.tsv", sep="\t")
    qc.to_sql("sample_qc", con, if_exists="replace", index=False)
    index = []
    for f in sorted(P.tables.glob("*.tsv")):
        name = f.stem
        df = pd.read_csv(f, sep="\t")
        df.to_sql(f"res_{name}", con, if_exists="replace", index=False)
        step, desc = DESCRIPTIONS.get(name, ("", ""))
        index.append({"table": f"res_{name}", "source_file": f"results/tables/{f.name}", "step": step,
                      "description": desc, "rows": len(df), "columns": df.shape[1]})
    pd.DataFrame(index).to_sql("results_index", con, if_exists="replace", index=False)
    con.commit()
    con.execute("VACUUM")
    con.close()
    missing = [i["table"] for i in index if not i["description"]]
    print(f"{len(index)} result tables loaded; database {P.db.stat().st_size / 1e6:.1f} MB")
    if missing:
        print("tables without description:", ", ".join(missing))


if __name__ == "__main__":
    main()
