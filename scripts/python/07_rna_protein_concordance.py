"""Step 07 - How well does mRNA abundance track protein abundance?

Matched samples only (RNA and protein from the same tissue, QC pass).

* gene-wise Spearman correlation across tumours (and across normal-adjacent
  tissues), with a sample-permutation null
* pooled tumour + normal correlation, to show how tissue differences inflate it
* partial correlation after removing tissue composition (ESTIMATE stromal and
  immune scores, acinar marker score), so that co-variation caused by how much
  stroma or acinar tissue a piece contains is not counted as "concordance"
* determinants of concordance (expression level, RNA and protein variability)
* sample-wise agreement of gene-centred profiles and its relation to purity
* gene-set level: Hallmark, Reactome and GO cellular-component complexes

Outputs: results/tables/concordance_*.tsv and results/figures/Fig3_concordance.*

    python scripts/python/07_rna_protein_concordance.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from cancerpg import io  # noqa: E402
from cancerpg.config import load_config, paths  # noqa: E402
from cancerpg.plotting import GREY, LAYER, TISSUE, panel_label, save_fig, set_style  # noqa: E402
from cancerpg.stats import bh, geneset_rank_test, partial_spearman, rowwise_spearman, zscore_rows  # noqa: E402

N_PERM = 20
COMPLEXES = {
    "GOCC_CYTOSOLIC_RIBOSOME": "Cytosolic ribosome",
    "GOCC_MITOCHONDRIAL_RIBOSOME": "Mitochondrial ribosome",
    "GOCC_PROTEASOME_COMPLEX": "Proteasome",
    "GOCC_RESPIRATORY_CHAIN_COMPLEX": "Respiratory chain",
    "GOCC_SPLICEOSOMAL_COMPLEX": "Spliceosome",
    "GOCC_MHC_CLASS_I_PROTEIN_COMPLEX": "MHC class I",
    "GOCC_COLLAGEN_TRIMER": "Collagen trimer",
    "GOCC_IMMUNOGLOBULIN_COMPLEX": "Immunoglobulin complex",
}


def main() -> None:
    cfg = load_config()
    P = paths(cfg)
    set_style()
    rng = np.random.default_rng(cfg["survival"]["seed"])

    an = io.analysis_samples()
    T = an[(an.tissue == "Tumor") & an.matched].sample_key.tolist()
    N = an[(an.tissue == "Normal") & an.matched].sample_key.tolist()
    rna = io.matrix("rna_vst.tsv.gz")
    prot = io.matrix("protein_norm.tsv.gz")

    min_frac = cfg["filters"]["protein_min_fraction_quantified"]
    genes = sorted(set(rna.index) & set(prot.index))
    genes = [g for g in genes if prot.loc[g, T].notna().mean() >= min_frac]
    rT, pT = rna.loc[genes, T], prot.loc[genes, T]
    print(f"genes with RNA and protein (protein in >= {min_frac:.0%} of {len(T)} tumours): {len(genes)}")

    # ------------------------------------------------------------------ gene-wise
    gT = rowwise_spearman(rT, pT)
    gT["fdr"] = bh(gT.p)
    gN = rowwise_spearman(rna.loc[genes, N], prot.loc[genes, N], min_n=15)
    gA = rowwise_spearman(rna.loc[genes, T + N], prot.loc[genes, T + N])

    null = []
    for _ in range(N_PERM):
        perm = list(rng.permutation(T))
        pp = pT[perm].copy()
        pp.columns = T
        null.append(rowwise_spearman(rT, pp).rho.to_numpy())
    null = np.concatenate(null)

    # composition covariates: ESTIMATE stromal + immune, acinar marker score
    pur = pd.read_csv(P.tables / "purity_estimates.tsv", sep="\t").set_index("sample_key")
    cov = pd.concat([pur[["stromal_score", "immune_score"]], an.set_index("sample_key")["acinar_rna"]], axis=1).loc[T]
    gT["rho_partial"] = partial_spearman(rT, pT, cov)

    gene_tab = pd.DataFrame({
        "gene": genes,
        "n_tumours": gT.n.values, "rho_tumour": gT.rho.values, "p_tumour": gT.p.values, "fdr_tumour": gT.fdr.values,
        "rho_tumour_partial": gT.rho_partial.values,
        "n_normals": gN.n.values, "rho_normal": gN.rho.values,
        "rho_pooled_T_N": gA.rho.values,
        "rna_mean": rT.mean(axis=1).values, "rna_sd": rT.std(axis=1).values,
        "protein_sd": pT.std(axis=1).values, "protein_missing_frac": pT.isna().mean(axis=1).values,
    })
    io.write_table(gene_tab, "concordance_gene.tsv")

    # ------------------------------------------------------------------ determinants
    det = []
    for v in ["rna_mean", "rna_sd", "protein_sd", "protein_missing_frac"]:
        r = stats.spearmanr(gene_tab[v], gene_tab.rho_tumour, nan_policy="omit")
        det.append({"feature": v, "spearman_with_rho": r.statistic, "p": r.pvalue})
    io.write_table(pd.DataFrame(det), "concordance_determinants.tsv")

    # ------------------------------------------------------------------ sample-wise
    both = [g for g in genes if pT.loc[g].notna().all()]
    zr, zp = zscore_rows(rT.loc[both]), zscore_rows(pT.loc[both])
    sw = pd.DataFrame({
        "sample_key": T,
        "rho_samplewise": [stats.spearmanr(zr[s], zp[s]).statistic for s in T],
    }).merge(pur[["kras_cellularity", "estimate_score"]], left_on="sample_key", right_index=True, how="left")
    sw = sw.merge(an[["sample_key", "match_level"]], on="sample_key", how="left")
    io.write_table(sw, "concordance_sample.tsv")
    # does a pooled proteome aliquot (2-3 tumour pieces) weaken the coupling?
    same, pooled = (sw.loc[sw.match_level == lv, "rho_samplewise"] for lv in ("same_sample", "pooled_sample"))

    # ------------------------------------------------------------------ gene sets
    rho = gene_tab.set_index("gene").rho_tumour
    hall = geneset_rank_test(rho, io.genesets(["HALLMARK"]))
    reac = geneset_rank_test(rho, io.genesets(["REACTOME"]), min_size=15, max_size=300)
    cc_sets = io.genesets(["GO_CC"])
    comp = geneset_rank_test(rho, {k: cc_sets[k] for k in COMPLEXES if k in cc_sets}, min_size=5)
    comp["label"] = comp.gene_set.map(COMPLEXES)
    hall_partial = geneset_rank_test(gene_tab.set_index("gene").rho_tumour_partial, io.genesets(["HALLMARK"]))
    io.write_table(hall, "concordance_hallmark.tsv")
    io.write_table(reac, "concordance_reactome.tsv")
    io.write_table(comp, "concordance_complexes.tsv")
    io.write_table(hall_partial, "concordance_hallmark_partial.tsv")

    # ------------------------------------------------------------------ summary
    summ = {
        "genes": len(genes), "tumours": len(T), "normals": len(N),
        "median_rho_tumour": np.nanmedian(gT.rho), "median_rho_normal": np.nanmedian(gN.rho),
        "median_rho_pooled": np.nanmedian(gA.rho), "median_rho_tumour_partial": np.nanmedian(gT.rho_partial),
        "median_rho_null": np.nanmedian(null), "null_95pct": np.nanquantile(null, 0.95),
        "frac_positive_fdr05": float(((gT.fdr < 0.05) & (gT.rho > 0)).mean()),
        "frac_negative_fdr05": float(((gT.fdr < 0.05) & (gT.rho < 0)).mean()),
        "median_samplewise_rho": sw.rho_samplewise.median(),
        "samplewise_vs_cellularity_rho": stats.spearmanr(sw.rho_samplewise, sw.kras_cellularity, nan_policy="omit").statistic,
        "samplewise_vs_cellularity_p": stats.spearmanr(sw.rho_samplewise, sw.kras_cellularity, nan_policy="omit").pvalue,
        "median_samplewise_rho_same_sample": same.median(), "median_samplewise_rho_pooled_sample": pooled.median(),
        "samplewise_same_vs_pooled_p": stats.mannwhitneyu(same, pooled).pvalue,
    }
    io.write_table(pd.Series(summ).rename("value").rename_axis("metric").reset_index(), "concordance_summary.tsv")
    for k, v in summ.items():
        print(f"{k:32s} {v:.4g}" if isinstance(v, float) else f"{k:32s} {v}")

    # ------------------------------------------------------------------ figure 3
    # two independent rows: the long Hallmark labels of panel D must not squeeze panel A
    fig = plt.figure(figsize=(7.6, 5.8))
    top = fig.add_gridspec(1, 3, left=0.075, right=0.985, top=0.93, bottom=0.60, wspace=0.32)
    bot = fig.add_gridspec(1, 3, left=0.19, right=0.985, top=0.44, bottom=0.075, wspace=0.92)
    axes = np.array([[fig.add_subplot(top[0, i]) for i in range(3)], [fig.add_subplot(bot[0, i]) for i in range(3)]])
    bins = np.linspace(-1, 1, 61)

    ax = axes[0, 0]
    ax.hist(null, bins=bins, density=True, color="0.85", label="permuted samples")
    ax.hist(gT.rho.dropna(), bins=bins, density=True, histtype="step", lw=1.2, color=TISSUE["Tumor"], label="tumours")
    ax.axvline(np.nanmedian(gT.rho), color=TISSUE["Tumor"], lw=0.8, ls="--")
    ax.set_xlabel("Gene-wise rho (mRNA vs protein)")
    ax.set_ylabel("Density")
    ax.set_title(f"{len(genes):,} genes x {len(T)} tumours")
    ax.text(0.02, 0.98, f"median {np.nanmedian(gT.rho):.2f}\n{summ['frac_positive_fdr05']:.0%} positive, FDR<5%",
            transform=ax.transAxes, va="top", fontsize=6.5)
    ax.set_ylim(0, 5.6)
    ax.legend(loc="upper right", fontsize=5.8)
    panel_label(ax, "A", x=-0.22, y=1.08)

    ax = axes[0, 1]
    for lab, vals, col, ls in [("tumours", gT.rho, TISSUE["Tumor"], "-"), ("normals (n=%d)" % len(N), gN.rho, TISSUE["Normal"], "-"),
                               ("tumour + normal pooled", gA.rho, "0.2", ":"),
                               ("tumours, composition-adjusted", gT.rho_partial, "#E08214", "--")]:
        v = vals.dropna()
        kde = stats.gaussian_kde(v)
        xs = np.linspace(-1, 1, 300)
        ax.plot(xs, kde(xs), color=col, lw=1.1, ls=ls, label=f"{lab}: {np.median(v):.2f}")
    ax.set_xlabel("Gene-wise Spearman rho")
    ax.set_ylabel("Density")
    ax.set_title("Tissue context")
    ax.set_ylim(0, 2.7)
    ax.legend(fontsize=5.3, loc="upper left", title="median rho", title_fontsize=5.3, handlelength=1.8)
    panel_label(ax, "B", x=-0.22, y=1.08)

    ax = axes[0, 2]
    gene_tab["rna_sd_decile"] = pd.qcut(gene_tab.rna_sd, 10, labels=False) + 1
    data = [gene_tab.loc[gene_tab.rna_sd_decile == d, "rho_tumour"].dropna() for d in range(1, 11)]
    ax.boxplot(data, widths=0.6, showfliers=False, medianprops=dict(color=TISSUE["Tumor"]),
               boxprops=dict(lw=0.6), whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6))
    ax.axhline(0, color="0.6", lw=0.5)
    ax.set_xlabel("mRNA SD decile (across tumours)")
    ax.set_ylabel("Spearman rho")
    ax.set_title("mRNA variability")
    ax.tick_params(axis="x", labelsize=6)
    panel_label(ax, "C", x=-0.22, y=1.08)

    ax = axes[1, 0]
    top = pd.concat([hall.head(8), hall.tail(8)])
    lab = top.gene_set.str.replace("HALLMARK_", "").str.replace("_", " ").str.lower().str.capitalize()
    cols = [TISSUE["Tumor"] if f < 0.05 else GREY for f in top.fdr]
    ax.barh(range(len(top))[::-1], top.median_in, color=cols, height=0.7)
    ax.axvline(np.nanmedian(gT.rho), color="0.3", lw=0.6, ls="--")
    ax.set_yticks(range(len(top))[::-1])
    ax.set_yticklabels(lab, fontsize=5.5)
    ax.set_xlabel("Median gene-wise rho")
    ax.set_title("Hallmark sets")
    panel_label(ax, "D", x=-1.05, y=1.08)

    ax = axes[1, 1]
    comp = comp.sort_values("median_in")
    ax.barh(comp.label, comp.median_in, color=[LAYER["Protein"] if f < 0.05 else GREY for f in comp.fdr], height=0.6)
    ax.axvline(np.nanmedian(gT.rho), color="0.3", lw=0.6, ls="--")
    for i, (m, s) in enumerate(zip(comp.median_in, comp["size"])):
        ax.text(max(m, 0) + 0.02, i, f"n={s}", va="center", fontsize=5.5)
    ax.set_xlabel("Median gene-wise rho")
    ax.set_title("Complexes / assemblies")
    ax.tick_params(axis="y", labelsize=6)
    panel_label(ax, "E", x=-0.95, y=1.08)

    ax = axes[1, 2]
    s = sw.dropna(subset=["kras_cellularity"])
    ax.scatter(s.kras_cellularity, s.rho_samplewise, s=7, color="0.25", lw=0)
    b = np.polyfit(s.kras_cellularity, s.rho_samplewise, 1)
    xs = np.linspace(s.kras_cellularity.min(), s.kras_cellularity.max(), 10)
    ax.plot(xs, np.polyval(b, xs), color=TISSUE["Tumor"], lw=0.8)
    ax.text(0.03, 0.97, f"rho = {summ['samplewise_vs_cellularity_rho']:.2f}\np = {summ['samplewise_vs_cellularity_p']:.1g}\nn = {len(s)}",
            transform=ax.transAxes, va="top", fontsize=6.5)
    ax.set_xlabel("Cellularity (2 x KRAS VAF)")
    ax.set_ylabel("Sample-wise rho (gene-centred)")
    ax.set_title("Per-tumour agreement")
    panel_label(ax, "F", x=-0.45, y=1.08)

    save_fig(fig, P.figures / "Fig3_concordance")


if __name__ == "__main__":
    main()
