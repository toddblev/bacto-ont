#!/usr/bin/env python3
"""
gc_vs_size.py - scatter of genome size vs GC% for the tutorial's reference genomes.

One picture that separates the strains by two simple, intuitive numbers. Genome
SIZE gives a clean ladder (L. casei < B. subtilis < E. coli DH5alpha < E. coli O157);
GC% separates the high-GC E. coli from the lower-GC Firmicutes. The two E. coli sit
almost on top of each other - you can't tell them apart here, which motivates the
read-mapping comparison that follows.

Reads the strain list from samples.tsv; computes size + GC straight from each
reference FASTA. Annotates each point with its known rrn operon count.

Usage:
    gc_vs_size.py [--samples samples.tsv] [--out analysis/gc_vs_size.png]
"""
import argparse
import csv

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import SeqIO


def genome_stats(fasta):
    """Return (total_bp, gc_percent) over all contigs in a FASTA."""
    g = c = a = t = 0
    for rec in SeqIO.parse(fasta, "fasta"):
        s = str(rec.seq).upper()
        g += s.count("G"); c += s.count("C"); a += s.count("A"); t += s.count("T")
    total = g + c + a + t
    return total, 100.0 * (g + c) / total if total else 0.0


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--samples", default="samples.tsv")
    ap.add_argument("--out", default="analysis/gc_vs_size.png")
    args = ap.parse_args()

    with open(args.samples) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    data = []
    for r in rows:
        size, gc = genome_stats(r["reference"])
        data.append((r["strain"], size / 1e6, gc, int(r["rrn_copies"])))
        print(f"{r['strain']:35s} {size/1e6:5.2f} Mb   GC {gc:4.1f}%   rrn {r['rrn_copies']}")

    plt.figure(figsize=(7.5, 5.5))
    for label, mb, gc, rrn in data:
        plt.scatter(mb, gc, s=90, zorder=3)
        plt.annotate(f"{label}\n({rrn} rrn operons)", (mb, gc),
                     textcoords="offset points", xytext=(8, 6), fontsize=8)
    gcs = [d[2] for d in data]
    mbs = [d[1] for d in data]
    plt.ylim(min(gcs) - 1.5, max(gcs) + 2.0)   # headroom so top labels are not clipped
    plt.xlim(min(mbs) - 0.6, max(mbs) + 1.2)
    plt.xlabel("genome size (Mb)")
    plt.ylabel("GC content (%)")
    plt.title("Four bacterial genomes: size vs GC%")
    plt.grid(alpha=0.3, zorder=0)
    plt.tight_layout()
    plt.savefig(args.out, dpi=130)
    print(f"-> {args.out}")


if __name__ == "__main__":
    main()
