#!/usr/bin/env python3
"""
rrn_operon_diagram.py - schematic of a typical E. coli rRNA (rrn) operon.

A teaching figure (not derived from data): the 16S-23S-5S gene layout transcribed
as one unit, with the tRNA gene(s) in the internal spacer. Sets up the copy-number
section (this ~5 kb unit is present in several copies per genome) and the Flye
repeat-collapse. Writes analysis/rrn_operon_diagram.png.

Usage:  rrn_operon_diagram.py [--out analysis/rrn_operon_diagram.png]
"""
import argparse
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, FancyArrowPatch


def gene_arrow(ax, x0, x1, y, h, color, label, sublabel=""):
    """Draw a right-pointing gene as a labelled pentagon arrow from x0 to x1."""
    head = min(180, (x1 - x0) * 0.35)
    body = x1 - x0 - head
    pts = [(x0, y - h/2), (x0 + body, y - h/2), (x1, y),
           (x0 + body, y + h/2), (x0, y + h/2)]
    ax.add_patch(Polygon(pts, closed=True, facecolor=color, edgecolor="black", lw=1.0))
    ax.text((x0 + x1)/2, y, label, ha="center", va="center",
            fontsize=9, fontweight="bold", color="white")
    if sublabel:
        ax.text((x0 + x1)/2, y - h/2 - 0.18, sublabel, ha="center", va="top", fontsize=8)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="analysis/rrn_operon_diagram.png")
    args = ap.parse_args()

    fig, ax = plt.subplots(figsize=(11, 2.6))
    y, h = 0, 0.5

    # backbone DNA line
    ax.plot([-350, 5450], [y, y], color="0.6", lw=1.2, zorder=0)

    # promoters (two tandem promoters P1/P2), drawn as bent arrows at the start
    for px in (-300, -190):
        ax.add_patch(FancyArrowPatch((px, y - 0.45), (px + 90, y - 0.45),
                                     arrowstyle="-|>", mutation_scale=12, color="black", lw=1.3))
    ax.text(-245, y - 0.78, "P1 P2\npromoters", ha="center", va="top", fontsize=8)

    # genes (approx E. coli coordinates, bp)
    gene_arrow(ax, 0,    1542, y, h, "#1f77b4", "16S", "rrs (~1.5 kb)")
    # tRNA gene(s) in the internal transcribed spacer
    ax.add_patch(Polygon([(1660, y-0.18),(1740, y-0.18),(1780, y),(1740, y+0.18),(1660, y+0.18)],
                         closed=True, facecolor="#2ca02c", edgecolor="black", lw=0.8))
    ax.text(1720, y + 0.30, "tRNA", ha="center", va="bottom", fontsize=8, color="#2ca02c")
    ax.text(1720, y - 0.40, "ITS", ha="center", va="top", fontsize=8, style="italic")
    gene_arrow(ax, 2050, 4954, y, h, "#ff7f0e", "23S", "rrl (~2.9 kb)")
    gene_arrow(ax, 5020, 5180, y, 0.34, "#d62728", "5S", "rrf")

    # terminator (stem-loop) after 5S
    ax.text(5300, y, "∩", ha="center", va="center", fontsize=20)
    ax.text(5300, y - 0.40, "term.", ha="center", va="top", fontsize=8)

    # direction-of-transcription arrow
    ax.annotate("", xy=(5250, y + 0.9), xytext=(0, y + 0.9),
                arrowprops=dict(arrowstyle="-|>", color="0.4", lw=1.2))
    ax.text(2600, y + 1.0, "transcribed as one ~5 kb pre-rRNA", ha="center", va="bottom",
            fontsize=8, color="0.4")

    # scale bar (1 kb)
    ax.plot([3800, 4800], [y - 1.15, y - 1.15], color="black", lw=2)
    ax.text(4300, y - 1.3, "1 kb", ha="center", va="top", fontsize=8)

    ax.set_xlim(-450, 5550)
    ax.set_ylim(-1.7, 1.5)
    ax.axis("off")
    ax.set_title("A typical E. coli rrn (rRNA) operon  —  ~5 kb, present in several copies per genome",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(args.out, dpi=140)
    print(f"-> {args.out}")


if __name__ == "__main__":
    main()
