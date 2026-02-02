from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

TISSUE = {"Tumor": "#B2182B", "Normal": "#2166AC"}
LAYER = {"RNA": "#1B7837", "Protein": "#762A83"}
GREY = "#7F7F7F"


def set_style() -> None:
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
        "font.size": 8,
        "axes.titlesize": 9,
        "axes.titleweight": "bold",
        "axes.labelsize": 8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "legend.frameon": False,
        "legend.fontsize": 7,
        "svg.fonttype": "none",
        "svg.hashsalt": "cancer-proteogenomics",
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.05,
    })


def panel_label(ax, label: str, x: float = -0.12, y: float = 1.06) -> None:
    ax.text(x, y, label, transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom", ha="left")


def save_fig(fig, base: Path) -> Path:
    base = Path(base)
    base.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(base.with_suffix(".svg"), metadata={"Date": None, "Creator": None})
    fig.savefig(base.with_suffix(".png"), dpi=300, facecolor="white", metadata={"Software": None})
    fig.savefig(base.with_suffix(".jpg"), dpi=300, facecolor="white", pil_kwargs={"quality": 95})
    plt.close(fig)
    return base
