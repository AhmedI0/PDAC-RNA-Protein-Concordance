"""Step 08 - Do DNA alterations reach the protein?

1. Copy-number dosage: for every gene, Spearman correlation of its copy number
   with its own mRNA and with its own protein across matched tumours. The drop
   from the CN-mRNA to the CN-protein correlation ("attenuation") measures how
   much protein levels are buffered against gene dosage; members of stable
   multi-protein complexes are expected to be buffered most.
2. Driver genes: mRNA and protein of KRAS, TP53, CDKN2A and SMAD4 by
   alteration status (somatic mutation class from WXS; homozygous deletion
   = ASCAT absolute copy number 0).

Outputs: results/tables/cis_cn_*.tsv, cis_driver_*.tsv; results/figures/Fig7_genomic_cis.*

    python scripts/python/08_genomic_cis_effects.py
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
from cancerpg.plotting import GREY, LAYER, panel_label, save_fig, set_style  # noqa: E402
from cancerpg.stats import bh, rowwise_spearman  # noqa: E402

TRUNC = {"Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Nonstop_Mutation",
         "Translation_Start_Site"}


def complex_members(max_size: int = 400) -> set[str]:
    cc = io.genesets(["GO_CC"])
    out: set[str] = set()
    for k, v in cc.items():
        if "COMPLEX" in k and len(v) <= max_size and "PROTEIN_CONTAINING_COMPLEX" not in k:
            out |= v
    return out


def driver_status(muts: pd.DataFrame, cn_abs: pd.DataFrame, genes: pd.DataFrame, gene: str,
                  cases: list[str]) -> pd.Series:
    """WT / missense / truncating / homdel (homozygous deletion wins, then truncating)."""
    st = pd.Series("WT", index=cases, dtype=object)
    m = muts[muts.Hugo_Symbol == gene]
    for c, g in m.groupby("case_submitter_id"):
        if c in st.index:
            st[c] = "truncating" if g.Variant_Classification.isin(TRUNC).any() else "missense/inframe"
    gid = genes.loc[genes.gene_name == gene, "gene_id"]
    if len(gid) and gid.iloc[0] in cn_abs.index:
        cn = cn_abs.loc[gid.iloc[0]]
        for c in cases:
            if c in cn.index and cn[c] == 0:
                st[c] = "homdel"
            elif c not in cn.index and st[c] == "WT":
                pass
    return st


def main() -> None:
    cfg = load_config()
    P = paths(cfg)
    set_style()
    an = io.analysis_samples()
    T = an[(an.tissue == "Tumor") & an.matched]
    case_of = dict(zip(T.sample_key, T.case_submitter_id))
    rna = io.matrix("rna_vst.tsv.gz")
    prot = io.matrix("protein_norm.tsv.gz")
    genes = pd.read_csv(P.processed / "genes.tsv", sep="\t")
    cn_log2 = io.matrix("cn_log2.tsv.gz")
    cn_abs = io.matrix("cn_absolute.tsv.gz")
    muts = pd.read_csv(P.processed / "mutations.tsv.gz", sep="\t")
    mm = pd.read_csv(P.processed / "mutation_matrix.tsv", sep="\t", index_col=0)

    # ------------------------------------------------------------------ CN dosage
    conc = pd.read_csv(P.tables / "concordance_gene.tsv", sep="\t")
    sym2id = genes.drop_duplicates("gene_name").set_index("gene_name")["gene_id"]
    use = conc[conc.gene.isin(sym2id.index)].gene
    use = [g for g in use if sym2id[g] in cn_log2.index]
    samples = [s for s in T.sample_key if case_of[s] in cn_log2.columns]
    cases = [case_of[s] for s in samples]
    cn = cn_log2.loc[[sym2id[g] for g in use], cases]
    cn.index = use
    cn.columns = samples
    var_ok = (cn.abs() > 0.2).sum(axis=1) >= 10           # copy-number change in >= 10 tumours
    cn = cn[var_ok]
    g = cn.index.tolist()
    r_rna = rowwise_spearman(cn, rna.loc[g, samples], min_n=30)
    r_pro = rowwise_spearman(cn, prot.loc[g, samples], min_n=30)
    cmem = complex_members()
    dos = pd.DataFrame({
        "gene": g, "n": r_rna.n.values,
        "rho_cn_rna": r_rna.rho.values, "fdr_cn_rna": bh(r_rna.p.values),
        "rho_cn_protein": r_pro.rho.values, "fdr_cn_protein": bh(r_pro.p.values),
        "complex_member": [x in cmem for x in g],
    })
    dos["attenuation"] = dos.rho_cn_rna - dos.rho_cn_protein
    dos["class"] = np.select(
        [(dos.fdr_cn_rna < 0.05) & (dos.rho_cn_rna > 0) & (dos.fdr_cn_protein < 0.05) & (dos.rho_cn_protein > 0),
         (dos.fdr_cn_rna < 0.05) & (dos.rho_cn_rna > 0)],
        ["transmitted to protein", "buffered at protein"], default="no RNA dosage effect")
    io.write_table(dos, "cis_cn_dosage.tsv")
    sig = dos[dos["class"] != "no RNA dosage effect"]
    mw = stats.mannwhitneyu(sig.loc[sig.complex_member, "attenuation"].dropna(),
                            sig.loc[~sig.complex_member, "attenuation"].dropna())
    summ = pd.DataFrame([
        {"metric": "tumours with CN", "value": len(samples)},
        {"metric": "genes tested", "value": len(dos)},
        {"metric": "median rho CN-RNA", "value": dos.rho_cn_rna.median()},
        {"metric": "median rho CN-protein", "value": dos.rho_cn_protein.median()},
        {"metric": "genes with RNA dosage effect (FDR<5%)", "value": len(sig)},
        {"metric": "of which transmitted to protein", "value": int((sig["class"] == "transmitted to protein").sum())},
        {"metric": "median attenuation, complex members", "value": sig.loc[sig.complex_member, "attenuation"].median()},
        {"metric": "median attenuation, other genes", "value": sig.loc[~sig.complex_member, "attenuation"].median()},
        {"metric": "Mann-Whitney p (complex vs other)", "value": mw.pvalue},
    ])
    io.write_table(summ, "cis_cn_summary.tsv")
    print(summ.to_string(index=False))

    # ------------------------------------------------------------------ drivers
    all_cases = T.case_submitter_id.tolist()
    rows, status_tab = [], {}
    for gene in ["KRAS", "TP53", "CDKN2A", "SMAD4"]:
        st = driver_status(muts, cn_abs, genes, gene, all_cases)
        st[~st.index.isin(mm.index)] = np.nan
        status_tab[gene] = st
        for layer, mat in [("RNA", rna), ("Protein", prot)]:
            if gene not in mat.index:
                continue
            v = pd.Series(mat.loc[gene, T.sample_key].values, index=all_cases)
            for grp in ["missense/inframe", "truncating", "homdel"]:
                a, b = v[st == grp].dropna(), v[st == "WT"].dropna()
                if len(a) >= 3 and len(b) >= 3:
                    t = stats.mannwhitneyu(a, b)
                    rows.append({"gene": gene, "layer": layer, "group": grp, "n_group": len(a), "n_wt": len(b),
                                 "median_diff": a.median() - b.median(), "p": t.pvalue})
    drv = pd.DataFrame(rows)
    drv["fdr"] = bh(drv.p)
    io.write_table(drv, "cis_driver_effects.tsv")
    st_df = pd.DataFrame(status_tab)
    st_df.index.name = "case_submitter_id"
    io.write_table(st_df.reset_index(), "driver_alteration_status.tsv")
    print(drv.to_string(index=False))

    # ------------------------------------------------------------------ figure 7
    fig = plt.figure(figsize=(7.4, 5.4))
    gs = fig.add_gridspec(2, 4, hspace=0.6, wspace=0.85)

    ax = fig.add_subplot(gs[0, :2])
    hb = ax.hexbin(dos.rho_cn_rna, dos.rho_cn_protein, gridsize=45, cmap="Greys", mincnt=1, bins="log", lw=0)
    ax.plot([-0.4, 1], [-0.4, 1], color="#B2182B", lw=0.7, ls="--")
    ax.set_xlabel("rho (copy number vs mRNA)")
    ax.set_ylabel("rho (copy number vs protein)")
    ax.set_title(f"Gene dosage ({len(dos):,} genes, {len(samples)} tumours)")
    ax.text(0.03, 0.97, f"median CN-mRNA {dos.rho_cn_rna.median():.2f}\nmedian CN-protein {dos.rho_cn_protein.median():.2f}",
            transform=ax.transAxes, va="top", fontsize=6.5)
    cb = fig.colorbar(hb, ax=ax, shrink=0.6, pad=0.02)
    cb.ax.set_title("genes", fontsize=6)
    cb.ax.tick_params(labelsize=5.5)
    panel_label(ax, "A", x=-0.18)

    ax = fig.add_subplot(gs[0, 2:])
    d1 = sig.loc[sig.complex_member, "attenuation"].dropna()
    d0 = sig.loc[~sig.complex_member, "attenuation"].dropna()
    parts = ax.violinplot([d0, d1], showmedians=True, widths=0.7)
    for pc, col in zip(parts["bodies"], [GREY, LAYER["Protein"]]):
        pc.set_facecolor(col)
        pc.set_alpha(0.6)
    for k in ("cmedians", "cmins", "cmaxes", "cbars"):
        parts[k].set_color("0.25")
        parts[k].set_linewidth(0.7)
    ax.set_xticks([1, 2])
    ax.set_xticklabels([f"other genes\n(n={len(d0)})", f"complex subunits\n(n={len(d1)})"])
    ax.set_ylabel("Attenuation\n(rho CN-mRNA - rho CN-protein)", fontsize=7)
    ax.set_title("Protein buffering of dosage")
    ax.text(0.5, 0.97, f"Mann-Whitney\np = {mw.pvalue:.1g}", transform=ax.transAxes, ha="center", va="top", fontsize=6.5)
    panel_label(ax, "B", x=-0.3)

    order = ["WT", "missense/inframe", "truncating", "homdel"]
    short = {"WT": "WT", "missense/inframe": "mis", "truncating": "trunc", "homdel": "HD"}
    for i, gene in enumerate(["TP53", "SMAD4", "CDKN2A", "KRAS"]):
        ax = fig.add_subplot(gs[1, i])
        st = status_tab[gene]
        groups = [o for o in order if (st == o).sum() >= 3]
        for j, (layer, mat, col) in enumerate([("RNA", rna, LAYER["RNA"]), ("Protein", prot, LAYER["Protein"])]):
            if gene not in mat.index:
                continue
            v = pd.Series(mat.loc[gene, T.sample_key].values, index=all_cases)
            z = (v - v[st == "WT"].mean()) / v[st == "WT"].std()
            data = [z[st == o].dropna() for o in groups]
            pos = np.arange(len(groups)) + (j - 0.5) * 0.36
            bp = ax.boxplot(data, positions=pos, widths=0.32, showfliers=False, patch_artist=True,
                            medianprops=dict(color="black", lw=0.8), boxprops=dict(lw=0.5),
                            whiskerprops=dict(lw=0.5), capprops=dict(lw=0.5))
            for b in bp["boxes"]:
                b.set_facecolor(col)
                b.set_alpha(0.75)
        ax.set_xticks(range(len(groups)))
        ax.set_xticklabels([f"{short[o]}\n{int((st == o).sum())}" for o in groups], fontsize=6)
        ax.axhline(0, color="0.6", lw=0.4)
        ax.set_title(gene, fontstyle="italic")
        if i == 0:
            ax.set_ylabel("Level vs WT tumours (z)")
            panel_label(ax, "C", x=-0.55)
    handles = [plt.Rectangle((0, 0), 1, 1, color=LAYER[k], alpha=0.75) for k in ["RNA", "Protein"]]
    fig.legend(handles, ["mRNA", "protein"], loc="lower center", ncol=2, bbox_to_anchor=(0.5, -0.01), fontsize=7)
    save_fig(fig, P.figures / "Fig7_genomic_cis")


if __name__ == "__main__":
    main()
