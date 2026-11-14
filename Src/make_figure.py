"""
make_figure.py  --  Step 3: Figure 1 (first draft).

Two stacked panels sharing one x-axis:
    A  DNA     (Parent / Designed / Sequenced), nt 592-624
    B  Protein (same three rows), codons 198-208
The mutation site (nt 607-609 / codon 203) is shaded in both panels.

Run:  python make_figure.py        -> writes ../Figures/figure1.png / .pdf / .svg
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

from load_data import load_all
from windowing import get_window

# ---------------------------------------------------------------- settings
OUT_DIR = Path(__file__).resolve().parent.parent / "Figures"

# Row labels on the left (placeholders: edit freely)
ROW_LABELS = {
    "parent": "Parent\n(template)",
    "designed": "Designed\n(T203Y)",
    "sequenced": "Sequenced\n(clone 1-1)",
}
ROWS = ["parent", "designed", "sequenced"]  # top to bottom

# Okabe-Ito colorblind-safe palette
C_NEUTRAL = "#4D4D4D"   # bases that match the parent / not at the site
C_SITE_WT = "#0072B2"   # parental codon (ACC / Thr) at the site
C_SITE_MUT = "#D55E00"  # designed mutant codon (TAT / Tyr) at the site
C_BAND = "#FFF3B0"      # shading behind the mutation site

MONO = "DejaVu Sans Mono"


def site_color(key):
    return C_SITE_MUT if key == "designed" else C_SITE_WT


def build_figure(w):
    n = len(w["dna"]["parent"])           # 33 nt
    n_codons = n // 3                      # 11
    site_lo = w["site_start"] - w["start"]  # 0-based index of first site base
    site_hi = w["site_end"] - w["start"]
    site_codon_idx = site_lo // 3

    fig, (axA, axB) = plt.subplots(
        2, 1, figsize=(8.2, 4.6), sharex=False,
        gridspec_kw={"height_ratios": [1, 1], "hspace": 0.55},
    )

    xmax = n + 11  # extra room on the right for call-outs
    for ax in (axA, axB):
        ax.set_xlim(-0.6, xmax)
        ax.set_ylim(-0.6, 2.6)
        ax.set_yticks([2, 1, 0])
        ax.set_yticklabels([ROW_LABELS[k] for k in ROWS], fontsize=9)
        ax.tick_params(axis="y", length=0)
        for s in ("top", "right", "left", "bottom"):
            ax.spines[s].set_visible(False)
        # shaded band through the mutation site
        ax.axvspan(site_lo - 0.5, site_hi + 0.5, color=C_BAND, zorder=0)

    # ------------------------------------------------ Panel A: DNA
    for r, key in enumerate(ROWS):
        y = 2 - r
        for i, base in enumerate(w["dna"][key]):
            at_site = site_lo <= i <= site_hi
            axA.text(i, y, base, ha="center", va="center", family=MONO,
                     fontsize=11,
                     color=site_color(key) if at_site else C_NEUTRAL,
                     fontweight="bold" if at_site else "normal")

    # nucleotide ruler on top of panel A
    nt_ticks = [w["start"], 597, 602, w["site_start"], 612, 617, 622]
    axA.xaxis.set_ticks_position("top")
    axA.set_xticks([t - w["start"] for t in nt_ticks])
    axA.set_xticklabels([str(t) for t in nt_ticks], fontsize=8)
    axA.xaxis.set_label_position("top")
    axA.set_xlabel("Position in GFP coding sequence (nt)", fontsize=9, labelpad=6)

    # call-outs on the right of panel A
    axA.text(n + 0.6, 1, "TAT  (expected)", va="center", fontsize=8.5,
             color=C_SITE_MUT, fontweight="bold")
    axA.text(n + 0.6, 0, "ACC  (not mutated)", va="center", fontsize=8.5,
             color=C_SITE_WT, fontweight="bold")

    # ------------------------------------------------ Panel B: protein
    for r, key in enumerate(ROWS):
        y = 2 - r
        for c, aa in enumerate(w["protein"][key]):
            x = c * 3 + 1  # middle base of the codon
            at_site = c == site_codon_idx
            axB.text(x, y, aa, ha="center", va="center", family=MONO,
                     fontsize=12,
                     color=site_color(key) if at_site else C_NEUTRAL,
                     fontweight="bold" if at_site else "normal")

    # codon-number ruler on top of panel B
    axB.xaxis.set_ticks_position("top")
    axB.set_xticks([c * 3 + 1 for c in range(n_codons)])
    axB.set_xticklabels([str(w["first_codon"] + c) for c in range(n_codons)], fontsize=8)
    axB.tick_params(axis="x", labelbottom=False, labeltop=True)
    axB.xaxis.set_label_position("top")
    axB.set_xlabel("Codon / residue number", fontsize=9, labelpad=6)

    axB.text(n + 0.6, 1, "Tyr203  (expected)", va="center", fontsize=8.5,
             color=C_SITE_MUT, fontweight="bold")
    axB.text(n + 0.6, 0, "Thr203  (unchanged)", va="center", fontsize=8.5,
             color=C_SITE_WT, fontweight="bold")

    # panel letters
    axA.text(-0.17, 1.28, "A", transform=axA.transAxes, fontsize=13, fontweight="bold")
    axB.text(-0.17, 1.28, "B", transform=axB.transAxes, fontsize=13, fontweight="bold")

    # legend
    handles = [
        Patch(facecolor=C_BAND, label="Mutation site (nt 607-609 / codon 203)"),
        Patch(facecolor=C_SITE_WT, label="Parental (ACC, Thr)"),
        Patch(facecolor=C_SITE_MUT, label="Designed mutant (TAT, Tyr)"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, fontsize=8.5,
               bbox_to_anchor=(0.5, 0.0))
    fig.subplots_adjust(left=0.17, right=0.98, top=0.88, bottom=0.1)
    return fig


if __name__ == "__main__":
    OUT_DIR.mkdir(exist_ok=True)
    fig = build_figure(get_window(load_all()))
    fig.savefig(OUT_DIR / "figure1.png", dpi=300)
    fig.savefig(OUT_DIR / "figure1.pdf")
    fig.savefig(OUT_DIR / "figure1.svg")
    print("Saved to", OUT_DIR)
