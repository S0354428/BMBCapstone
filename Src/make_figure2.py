"""
make_figure2.py  --  Figure 2 (first draft).

    A  BSA standard curve: blank-corrected A595 vs [BSA] (mg/mL), linear trendline,
       equation and R2.
    B  Protein concentration (mg/mL) for Crude, E1, E2 and FT, induced vs uninduced.
       BDL = below detection.

Run:  python make_figure2.py     -> ../Figures/figure2.png / .pdf / .svg
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

from fit_standards import equation_text, fit_standards
from load_bradford import load_standards
from quantify import quantify

OUT_DIR = Path(__file__).resolve().parent.parent / "Figures"

# Okabe-Ito colorblind-safe palette
C_INDUCED = "#0072B2"
C_UNINDUCED = "#E69F00"
C_FIT = "#222222"
C_POINTS = "#9A9A9A"

GROUPS = ["Crude", "E1", "E2", "FT"]
GROUP_LABELS = {"Crude": "Crude\nsupernatant", "E1": "E1", "E2": "E2", "FT": "Flow-\nthrough"}


def panel_a(ax):
    f = fit_standards()
    std = load_standards()
    x, y = f["x"], f["y"]

    # individual replicates (light) and means (solid)
    reps = std[["rep1", "rep2", "rep3"]].to_numpy() - f["blank"]
    ax.scatter(np.repeat(x, 3), reps.flatten(), s=14, color=C_POINTS, zorder=2,
               label="Individual replicates")
    ax.scatter(x, y, s=38, color=C_INDUCED, edgecolor="white", linewidth=0.8,
               zorder=3, label="Mean of replicates")

    xx = np.linspace(0, 1.0, 50)
    ax.plot(xx, f["slope"] * xx + f["intercept"], color=C_FIT, linewidth=1.6,
            zorder=1, label="Linear trendline")

    ax.text(0.04, 0.97, f"{equation_text(f)}\n$R^2$ = {f['r2']:.3f}",
            transform=ax.transAxes, va="top", ha="left", fontsize=9.5,
            bbox=dict(boxstyle="round,pad=0.35", facecolor="white", edgecolor="#BBBBBB"))

    ax.set_title("BSA standard curve (Bradford assay)", fontsize=11, fontweight="bold")
    ax.set_xlabel("BSA concentration (mg/mL)", fontsize=10)
    ax.set_ylabel("Blank-corrected absorbance (595 nm)", fontsize=10)
    ax.set_xlim(-0.04, 1.06)
    ax.set_ylim(-0.1, 2.0)
    ax.legend(loc="lower right", fontsize=8, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)


def panel_b(ax, rep_df):
    width = 0.38
    xpos = np.arange(len(GROUPS))
    ymax = 0

    for offset, cond, color in ((-width / 2, "Induced", C_INDUCED),
                                (width / 2, "Uninduced", C_UNINDUCED)):
        for i, grp in enumerate(GROUPS):
            r = rep_df[(rep_df["condition"] == cond) & (rep_df["group"] == grp)].iloc[0]
            x = xpos[i] + offset
            if r["status"] == "below detection":
                ax.text(x, 0.12, "BDL", ha="center", va="bottom", fontsize=8.5, color="#444444")
                continue
            ax.bar(x, r["conc_mg_mL"], width=width * 0.94, color=color, edgecolor="white",
                   linewidth=0.6, zorder=2)
            ax.errorbar(x, r["conc_mg_mL"], yerr=r["sd_mg_mL"], color="#333333",
                        capsize=3, linewidth=1.1, zorder=3)
            top = r["conc_mg_mL"] + r["sd_mg_mL"]
            ax.text(x, top + 0.12, f"{r['conc_mg_mL']:.2f}", ha="center", va="bottom", fontsize=8.5)
            ymax = max(ymax, top)

    ax.set_title("Protein concentration of crude lysate and purification fractions",
                 fontsize=11, fontweight="bold")
    ax.set_ylabel("Protein concentration (mg/mL)", fontsize=10)
    ax.set_xticks(xpos)
    ax.set_xticklabels([GROUP_LABELS[g] for g in GROUPS], fontsize=9.5)
    ax.set_ylim(0, ymax * 1.18)
    ax.spines[["top", "right"]].set_visible(False)

    handles = [
        Patch(facecolor=C_INDUCED, label="Induced"),
        Patch(facecolor=C_UNINDUCED, label="Uninduced"),
    ]
    ax.legend(handles=handles, loc="upper right", fontsize=8.5, frameon=False)
    ax.text(0.99, 0.84, "BDL = below detection", transform=ax.transAxes,
            ha="right", va="top", fontsize=8, color="#444444")


def build_figure():
    _, rep_df = quantify()
    fig, (axA, axB) = plt.subplots(
        1, 2, figsize=(11.5, 4.8), gridspec_kw={"width_ratios": [1, 1.35], "wspace": 0.28})
    panel_a(axA)
    panel_b(axB, rep_df)
    axA.text(-0.14, 1.08, "A", transform=axA.transAxes, fontsize=15, fontweight="bold")
    axB.text(-0.12, 1.08, "B", transform=axB.transAxes, fontsize=15, fontweight="bold")
    fig.subplots_adjust(left=0.07, right=0.985, top=0.88, bottom=0.14)
    return fig


if __name__ == "__main__":
    OUT_DIR.mkdir(exist_ok=True)
    fig = build_figure()
    fig.savefig(OUT_DIR / "figure2.png", dpi=300)
    fig.savefig(OUT_DIR / "figure2.pdf")
    fig.savefig(OUT_DIR / "figure2.svg")
    print("Saved to", OUT_DIR)
