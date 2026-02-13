"""Step 06 - Cohort overview: Table 1 and Figure 1 (study design, data
availability, RNA-protein match levels, overall survival).

    python scripts/python/06_cohort_overview.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines.plotting import add_at_risk_counts
from lifelines.utils import median_survival_times
from matplotlib.patches import FancyBboxPatch

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from cancerpg import io  # noqa: E402
from cancerpg.config import load_config, paths  # noqa: E402
from cancerpg.plotting import GREY, LAYER, TISSUE, panel_label, save_fig, set_style  # noqa: E402

DAYS_PER_MONTH = 30.4375


def n_pct(s: pd.Series, total: int) -> str:
    return f"{int(s)} ({100 * s / total:.0f}%)"


def km_median(t, e) -> tuple[float, float, float]:
    k = KaplanMeierFitter().fit(t, e)
    ci = median_survival_times(k.confidence_interval_)
    return k.median_survival_time_, float(ci.iloc[0, 0]), float(ci.iloc[0, 1])


def table1(df: pd.DataFrame, mm: pd.DataFrame) -> pd.DataFrame:
    n = len(df)
    rows = [("Patients", "", str(n))]
    q = df.age_years.quantile([0.5, 0.25, 0.75])
    rows.append(("Age at diagnosis, years", "median (IQR)", f"{q.iloc[0]:.0f} ({q.iloc[1]:.0f}-{q.iloc[2]:.0f})"))
    for v, c in df.sex.fillna("not reported").value_counts().items():
        rows.append(("Sex", v, n_pct(c, n)))
    for v in ["I", "II", "III", "IV"]:
        rows.append(("AJCC pathological stage", v, n_pct((df.stage_group == v).sum(), n)))
    rows.append(("AJCC pathological stage", "not reported", n_pct(df.stage_group.isna().sum(), n)))
    for v, c in df.ajcc_edition.fillna("not reported").value_counts().items():
        rows.append(("AJCC edition", v, n_pct(c, n)))
    for v, c in df.tumor_grade.fillna("not reported").value_counts().sort_index().items():
        rows.append(("Tumour grade", v, n_pct(c, n)))
    for v, c in df.residual_disease.fillna("not reported").value_counts().sort_index().items():
        rows.append(("Resection margin", v, n_pct(c, n)))
    for v, c in df.tumor_site.fillna("not reported").value_counts().items():
        rows.append(("Tumour site", v, n_pct(c, n)))
    q = df.tumor_size_cm.quantile([0.5, 0.25, 0.75])
    rows.append(("Tumour size, cm", "median (IQR)", f"{q.iloc[0]:.1f} ({q.iloc[1]:.1f}-{q.iloc[2]:.1f})"))
    rows.append(("Perineural invasion", "present", n_pct((df.perineural_invasion == "Yes").sum(), n)))
    s = df.dropna(subset=["os_days"])
    rows.append(("Overall survival", "patients with follow-up", str(len(s))))
    rows.append(("Overall survival", "deaths", n_pct(int(s.os_event.sum()), len(s))))
    med, lo, hi = km_median(s.os_days / DAYS_PER_MONTH, s.os_event)
    rows.append(("Overall survival", "median, months (95% CI)", f"{med:.1f} ({lo:.1f}-{hi:.1f})"))
    fmed, flo, fhi = km_median(s.os_days / DAYS_PER_MONTH, 1 - s.os_event)
    rows.append(("Follow-up (reverse Kaplan-Meier)", "median, months (95% CI)", f"{fmed:.1f} ({flo:.1f}-{fhi:.1f})"))
    rows.append(("Recurrence or progression", "reported", n_pct(int(df.recurrence.sum()), n)))
    m = mm.reindex(df.case_submitter_id)
    nm = int(m.KRAS.notna().sum())
    rows.append(("Somatic mutation (WXS)", "patients with WXS", str(nm)))
    for g in ["KRAS", "TP53", "CDKN2A", "SMAD4"]:
        rows.append(("Somatic mutation (WXS)", g, n_pct(int(m[g].sum()), nm)))
    return pd.DataFrame(rows, columns=["characteristic", "level", "value"])


def draw_design(ax, n: dict) -> None:
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis("off")

    def box(x, y, w, h, text, fc):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.08,rounding_size=0.25",
                                    fc=fc, ec="0.3", lw=0.6))
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=6.6, linespacing=1.25)

    def arrow(x0, y0, x1, y1):
        ax.annotate("", xy=(x1, y1), xytext=(x0, y0), arrowprops=dict(arrowstyle="-|>", lw=0.7, color="0.3"))

    box(0.1, 6.9, 4.6, 2.9, f"GDC  -  CPTAC-3, pancreas\nRNA-seq STAR counts (n={n['rna']})\n"
        f"WXS somatic mutations (n={n['maf']})\nWGS gene copy number (n={n['cn']})\n"
        f"clinical follow-up (n={n['cases']} cases)", "#E5F2E5")
    box(5.3, 6.9, 4.6, 2.9, f"PDC  -  PDC000270\nTMT11 proteome, 25 plexes\nCDAP gene-level log2 ratios\n"
        f"{n['proteins']:,} proteins x {n['prot']} samples", "#EFE5F2")
    box(1.2, 3.75, 7.6, 2.3, "SQLite sample database: case > sample > aliquot\n"
        "pooled proteome aliquots resolved to tissue pieces\nQC: duodenal mucosa, TMT plex, purity",
        "#F2F2F2")
    box(0.1, 0.2, 9.8, 2.75, f"Matched cohort: {n['mt']} tumours + {n['mn']} normal-adjacent tissues\n"
        f"(RNA and protein from the same tissue)\n\nRNA-protein concordance | tumour vs normal | pathways\n"
        f"subtypes | copy number and mutations | overall survival", "#FBE9E7")
    arrow(2.4, 6.9, 3.6, 6.1)
    arrow(7.6, 6.9, 6.4, 6.1)
    arrow(5.0, 3.75, 5.0, 3.0)


def main() -> None:
    cfg = load_config()
    P = paths(cfg)
    set_style()
    an = io.analysis_samples()
    cases = io.cases()
    mm = pd.read_csv(P.processed / "mutation_matrix.tsv", sep="\t", index_col=0)
    cn = pd.read_csv(P.processed / "cn_log2.tsv.gz", sep="\t", index_col=0, nrows=1)
    prot = pd.read_csv(P.processed / "protein_norm.tsv.gz", sep="\t", index_col=0, usecols=[0])

    pdac = cases[cases.primary_diagnosis == "Infiltrating duct carcinoma, NOS"]
    mt = an[(an.tissue == "Tumor") & an.matched]
    mn = an[(an.tissue == "Normal") & an.matched]
    coh = pdac[pdac.case_submitter_id.isin(mt.case_submitter_id)].copy()

    t1 = table1(coh, mm)
    io.write_table(t1, "table1_cohort.tsv")
    print(t1.to_string(index=False))

    # ---------------------------------------------------------------- availability
    qc = an[an.qc_pass]
    avail = pd.DataFrame({
        "RNA tumour": [qc[(qc.tissue == "Tumor") & (qc.has_rna == 1)].case_submitter_id.nunique()],
        "Protein tumour": [qc[(qc.tissue == "Tumor") & (qc.has_protein == 1)].case_submitter_id.nunique()],
        "RNA normal": [qc[(qc.tissue == "Normal") & (qc.has_rna == 1)].case_submitter_id.nunique()],
        "Protein normal": [qc[(qc.tissue == "Normal") & (qc.has_protein == 1)].case_submitter_id.nunique()],
        "WXS mutations": [int(mm.index.isin(pdac.case_submitter_id).sum())],
        "WGS copy number": [int(pd.Index(cn.columns).isin(pdac.case_submitter_id).sum())],
        "Overall survival": [int(pdac.os_days.notna().sum())],
    }).T[0]
    io.write_table(avail.rename("cases").reset_index().rename(columns={"index": "data"}), "cohort_data_availability.tsv")

    counts = {"rna": int((an.has_rna == 1).sum()), "maf": int(mm.index.isin(pdac.case_submitter_id).sum()),
              "cn": int(pd.Index(cn.columns).isin(pdac.case_submitter_id).sum()), "cases": len(pdac),
              "prot": int((an.has_protein == 1).sum()), "proteins": len(prot), "mt": len(mt), "mn": len(mn)}

    # ---------------------------------------------------------------- figure 1
    fig = plt.figure(figsize=(7.2, 6.2))
    gs = fig.add_gridspec(2, 3, height_ratios=[0.75, 1], width_ratios=[1, 0.62, 1.45], hspace=0.35, wspace=0.62)
    ax = fig.add_subplot(gs[0, :])
    draw_design(ax, counts)
    ax.set_xlim(-0.6, 10.6)
    panel_label(ax, "A", x=-0.06, y=0.96)

    ax = fig.add_subplot(gs[1, 0])
    cols = [LAYER["RNA"], LAYER["Protein"], LAYER["RNA"], LAYER["Protein"], GREY, GREY, "#444444"]
    ax.barh(avail.index[::-1], avail.values[::-1], color=cols[::-1], height=0.65)
    for i, v in enumerate(avail.values[::-1]):
        ax.text(v + 2, i, str(v), va="center", fontsize=6.5)
    ax.set_xlabel("PDAC cases")
    ax.set_title("Data available (after QC)")
    ax.set_xlim(0, avail.max() * 1.18)
    panel_label(ax, "B", x=-0.95)

    ax = fig.add_subplot(gs[1, 1])
    lv = (an[an.matched].groupby(["tissue", "match_level"]).size().unstack(fill_value=0)
          .reindex(["Tumor", "Normal"]))
    lv = lv.reindex(columns=["same_sample", "pooled_sample"], fill_value=0)
    ax.bar(lv.index, lv["same_sample"], color="#4D4D4D", label="same tissue piece", width=0.6)
    ax.bar(lv.index, lv["pooled_sample"], bottom=lv["same_sample"], color="#BABABA",
           label="proteome pool contains\nthe RNA piece", width=0.6)
    for i, t in enumerate(lv.index):
        ax.text(i, lv.loc[t].sum() + 2, str(int(lv.loc[t].sum())), ha="center", fontsize=7)
    ax.set_ylabel("Matched RNA-protein samples")
    ax.set_title("Match level")
    ax.set_ylim(0, 150)
    ax.legend(loc="upper center", bbox_to_anchor=(0.45, -0.13), fontsize=5.8, handlelength=1)
    panel_label(ax, "C", x=-0.5)

    ax = fig.add_subplot(gs[1, 2])
    s = coh.dropna(subset=["os_days"])
    k = KaplanMeierFitter(label="n").fit(s.os_days / DAYS_PER_MONTH, s.os_event)
    k.plot_survival_function(ax=ax, color="#B2182B", ci_alpha=0.15, show_censors=True,
                             censor_styles={"ms": 4, "marker": "|"})
    med, lo, hi = km_median(s.os_days / DAYS_PER_MONTH, s.os_event)
    ax.axhline(0.5, color="0.6", lw=0.5, ls=":")
    ax.text(0.98, 0.95, f"n = {len(s)}, deaths = {int(s.os_event.sum())}\nmedian {med:.1f} months\n"
            f"(95% CI {lo:.1f}-{hi:.1f})", transform=ax.transAxes, ha="right", va="top", fontsize=7)
    ax.set_xlabel("Months from diagnosis")
    ax.set_ylabel("Overall survival")
    ax.set_title("Overall survival")
    ax.set_xlim(0, 72)
    ax.set_xticks(range(0, 73, 12))
    ax.get_legend().remove()
    add_at_risk_counts(k, ax=ax, rows_to_show=["At risk"], xticks=list(range(0, 73, 12)), fontsize=6.5)
    panel_label(ax, "D", x=-0.1)

    save_fig(fig, P.figures / "Fig1_cohort")
    print(f"\nmatched tumours {len(mt)}, normals {len(mn)}; survival cohort {len(s)} ({int(s.os_event.sum())} deaths)")


if __name__ == "__main__":
    main()
