"""Step 04 - Assemble analysis matrices from the downloaded files.

Samples are keyed "<case>_T" / "<case>_N" (one RNA file and/or one proteome
aliquot per case and tissue, chosen in step 03). Only PDAC cases are kept.

Outputs (data/processed/):
  samples.tsv                 sample annotation (layers present, match level, TMT plex, ...)
  genes.tsv                   GENCODE v36 gene annotation (+ chromosome position from the CN files)
  rna_counts.tsv.gz           STAR gene counts, strand-specific column chosen per library
  rna_tpm.tsv.gz              TPM (unstranded, as provided by GDC)
  rna_strandedness.tsv        per-library strand call
  protein_log2ratio.tsv.gz    CDAP log2 ratio to the pooled reference (as downloaded)
  protein_norm.tsv.gz         per-sample median-centred log2 ratios
  mutations.tsv.gz            protein-altering somatic variants (one tumour MAF per case)
  mutation_matrix.tsv         case x driver gene (0/1) + KRAS allele and VAF
  cn_log2.tsv.gz              gene-level log2(copy number / sample median copy number)
  cn_absolute.tsv.gz          gene-level absolute copy number (ASCAT)

    python scripts/python/04_assemble_matrices.py
"""
from __future__ import annotations

import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from cancerpg.config import load_config, paths  # noqa: E402

NONSYN = {
    "Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins",
    "In_Frame_Del", "In_Frame_Ins", "Splice_Site", "Translation_Start_Site", "Nonstop_Mutation",
}
MAF_COLS = ["Hugo_Symbol", "Chromosome", "Start_Position", "Variant_Classification", "Variant_Type",
            "Reference_Allele", "Tumor_Seq_Allele2", "HGVSp_Short", "t_depth", "t_alt_count",
            "Tumor_Sample_Barcode", "IMPACT", "callers"]


def read_star(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", comment="#")
    return df[df.gene_id.str.startswith("ENSG")].set_index("gene_id")


def strand_call(df: pd.DataFrame) -> tuple[str, float]:
    first, second = df["stranded_first"].sum(), df["stranded_second"].sum()
    ratio = (second + 1) / (first + 1)
    if ratio > 5:
        return "stranded_second", ratio
    if ratio < 0.2:
        return "stranded_first", ratio
    return "unstranded", ratio


def main() -> None:
    cfg = load_config()
    P = paths(cfg)
    gdir, pdir, out = P.raw / "gdc", P.raw / "pdc", P.processed
    out.mkdir(parents=True, exist_ok=True)

    cases = pd.read_csv(P.interim / "cases.tsv", sep="\t")
    m = pd.read_csv(P.interim / "sample_matches.tsv", sep="\t")
    files = pd.read_csv(P.interim / "gdc_files.tsv", sep="\t")
    aliq = pd.read_csv(P.interim / "pdc_aliquots.tsv", sep="\t")
    pdac = set(cases.loc[cases.primary_diagnosis == "Infiltrating duct carcinoma, NOS", "case_submitter_id"])
    pdac_pdc = set(aliq.loc[aliq.disease_type == "Pancreatic Ductal Adenocarcinoma", "aliquot_submitter_id"])

    m = m[m.case_submitter_id.isin(pdac) & (m.protein_aliquot.isna() | m.protein_aliquot.isin(pdac_pdc))].copy()

    # ------------------------------------------------------------------ RNA
    rna_rows = m[m.rna_file_id.notna()]
    fmeta = files[files.kind == "rna"].drop_duplicates("file_id").set_index("file_id")
    counts, tpm, strands, genes = {}, {}, [], None
    for _, r in rna_rows.iterrows():
        f = fmeta.loc[r.rna_file_id]
        df = read_star(gdir / "rna" / r.rna_file_id / f.file_name)
        col, ratio = strand_call(df)
        counts[r.sample_key] = df[col].astype(np.int64)
        tpm[r.sample_key] = df["tpm_unstranded"].astype(float)
        strands.append({"sample_key": r.sample_key, "count_column": col, "second_to_first_ratio": round(ratio, 2),
                        "assigned_reads": int(df[col].sum())})
        if genes is None:
            genes = df[["gene_name", "gene_type"]].copy()
    counts = pd.DataFrame(counts)
    tpm = pd.DataFrame(tpm)
    strands = pd.DataFrame(strands)

    # gene positions from one ASCAT gene-level CN file (same GENCODE v36 gene IDs)
    # one row per (file, sample): keep the tumour row of each tumour/normal WGS pair
    cn_meta = files[files.kind == "cn"].drop_duplicates("file_id")
    first_cn = cn_meta.iloc[0]
    pos = pd.read_csv(gdir / "cn" / first_cn.file_id / first_cn.file_name, sep="\t",
                      usecols=["gene_id", "chromosome", "start", "end"]).set_index("gene_id")
    genes = genes.join(pos, how="left")
    genes.index.name = "gene_id"
    genes["gene_id_stable"] = genes.index.str.replace(r"\.\d+$", "", regex=True)

    # ------------------------------------------------------------------ protein
    prot_raw = pd.read_csv(pdir / "proteome_log2ratio.tsv.gz", sep="\t", index_col=0)
    prot_rows = m[m.protein_column.notna()]
    prot = prot_raw[prot_rows.protein_column.tolist()].copy()
    prot.columns = prot_rows.sample_key.tolist()
    prot = prot.loc[prot.notna().any(axis=1)]
    prot_norm = prot - prot.median(axis=0, skipna=True)

    # ------------------------------------------------------------------ samples
    samples = m[["sample_key", "case_submitter_id", "tissue", "match_level", "included",
                 "rna_aliquot", "rna_sample", "protein_aliquot", "protein_samples"]].copy()
    samples["has_rna"] = samples.sample_key.isin(counts.columns).astype(int)
    samples["has_protein"] = samples.sample_key.isin(prot.columns).astype(int)
    samples["n_samples_pooled"] = samples.protein_samples.fillna("").map(lambda s: len([x for x in s.split(";") if x]))
    plex = aliq.set_index("aliquot_submitter_id")[["plex", "tmt_channel"]]
    samples = samples.join(plex, on="protein_aliquot")
    samples = samples.merge(strands[["sample_key", "count_column", "assigned_reads"]], on="sample_key", how="left")
    samples["protein_n_quantified"] = samples.sample_key.map(prot.notna().sum())

    # ------------------------------------------------------------------ mutations
    maf_meta = files[(files.kind == "maf") & (files.tissue == "Tumor")]
    tumor_rna_sample = m[m.tissue == "Tumor"].set_index("case_submitter_id")["rna_sample"].to_dict()
    muts = []
    for case, g in maf_meta.groupby("case_submitter_id"):
        if case not in pdac:
            continue
        # prefer the WXS tumour sample that was also used for RNA-seq, then a stable order
        g = g.assign(pref=(g.sample_submitter_id != tumor_rna_sample.get(case)).astype(int)).sort_values(["pref", "file_id"])
        f = g.iloc[0]
        with gzip.open(gdir / "maf" / f.file_id / f.file_name, "rt") as fh:
            df = pd.read_csv(fh, sep="\t", comment="#", usecols=MAF_COLS, low_memory=False)
        df = df[df.Variant_Classification.isin(NONSYN)].copy()
        df.insert(0, "case_submitter_id", case)
        df["maf_file_id"] = f.file_id
        df["vaf"] = df.t_alt_count / df.t_depth
        if len(df):
            muts.append(df)
    muts = pd.concat(muts, ignore_index=True)
    maf_cases = sorted(muts.case_submitter_id.unique().tolist()
                       + [c for c in maf_meta.case_submitter_id.unique() if c in pdac and c not in set(muts.case_submitter_id)])
    drivers = cfg["drivers"]
    mm = pd.DataFrame(0, index=pd.Index(sorted(set(maf_cases)), name="case_submitter_id"), columns=drivers)
    for gname in drivers:
        hit = muts.loc[muts.Hugo_Symbol == gname, "case_submitter_id"].unique()
        mm.loc[mm.index.intersection(hit), gname] = 1
    kras = (muts[muts.Hugo_Symbol == "KRAS"].sort_values("vaf", ascending=False)
            .drop_duplicates("case_submitter_id").set_index("case_submitter_id"))
    mm["KRAS_allele"] = kras["HGVSp_Short"].reindex(mm.index)
    mm["KRAS_vaf"] = kras["vaf"].reindex(mm.index)
    mm["tmb_nonsyn"] = muts.groupby("case_submitter_id").size().reindex(mm.index).fillna(0).astype(int)

    # ------------------------------------------------------------------ copy number
    cn_abs = {}
    for case, g in cn_meta.groupby("case_submitter_id"):
        if case not in pdac:
            continue
        g = g.assign(pref=(g.sample_submitter_id != tumor_rna_sample.get(case)).astype(int)).sort_values(["pref", "file_id"])
        f = g.iloc[0]
        df = pd.read_csv(gdir / "cn" / f.file_id / f.file_name, sep="\t", usecols=["gene_id", "copy_number"])
        cn_abs[case] = df.set_index("gene_id")["copy_number"]
    cn_abs = pd.DataFrame(cn_abs)
    med = cn_abs.median(axis=0)
    cn_log2 = np.log2((cn_abs + 0.5) / (med + 0.5))

    # ------------------------------------------------------------------ write
    samples.to_csv(out / "samples.tsv", sep="\t", index=False)
    genes.to_csv(out / "genes.tsv", sep="\t")
    counts.to_csv(out / "rna_counts.tsv.gz", sep="\t", compression="gzip")
    tpm.to_csv(out / "rna_tpm.tsv.gz", sep="\t", float_format="%.4f", compression="gzip")
    strands.to_csv(out / "rna_strandedness.tsv", sep="\t", index=False)
    prot.to_csv(out / "protein_log2ratio.tsv.gz", sep="\t", float_format="%.5f", compression="gzip")
    prot_norm.to_csv(out / "protein_norm.tsv.gz", sep="\t", float_format="%.5f", compression="gzip")
    muts.to_csv(out / "mutations.tsv.gz", sep="\t", index=False, compression="gzip")
    mm.to_csv(out / "mutation_matrix.tsv", sep="\t")
    cn_abs.to_csv(out / "cn_absolute.tsv.gz", sep="\t", compression="gzip")
    cn_log2.to_csv(out / "cn_log2.tsv.gz", sep="\t", float_format="%.4f", compression="gzip")

    print(f"RNA: {counts.shape[0]} genes x {counts.shape[1]} samples; strand calls {strands.count_column.value_counts().to_dict()}")
    print(f"protein: {prot.shape[0]} genes x {prot.shape[1]} samples; complete in all samples: {prot.notna().all(axis=1).sum()}")
    print(f"mutations: {len(muts)} protein-altering variants in {muts.case_submitter_id.nunique()} cases")
    print("driver frequency:", (mm[drivers].mean() * 100).round(1).to_dict())
    print(f"CN: {cn_abs.shape[0]} genes x {cn_abs.shape[1]} tumours")
    print(samples.groupby(["tissue", "has_rna", "has_protein"]).size().to_string())


if __name__ == "__main__":
    main()
