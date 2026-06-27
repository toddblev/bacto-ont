#!/usr/bin/env python3
"""
cross_mapping_matrix.py - map every barcode against every reference genome.

Builds a barcode x reference table of "% of reads that map", which:
  * identifies each barcode (its reads map overwhelmingly to ONE genome -> the diagonal),
  * shows genus specificity (L. casei and B. subtilis map ~only to themselves), and
  * shows that the two E. coli cross-map to BOTH E. coli references (same species) -
    so "% mapped" alone cannot separate strains of one species.

Reads the barcode list + references from samples.tsv. Expects each barcode's reads
at <reads-dir>/<barcode>.fastq (the notebook stages them there after demux).

Usage:
    cross_mapping_matrix.py [--samples samples.tsv] [--reads-dir reads]
                            [--out analysis/cross_mapping_matrix.png] [--threads 8]
"""
import argparse
import csv
import os
import subprocess

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def count_reads(fastq):
    out = subprocess.run(["wc", "-l", fastq], capture_output=True, text=True).stdout
    return int(out.split()[0]) // 4


def pct_mapped(fastq, ref, threads):
    """Primary-mapped reads / total reads, as a percentage."""
    mm = subprocess.Popen(
        ["minimap2", "-ax", "map-ont", "--secondary=no", "-t", str(threads), ref, fastq],
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    # -F 0x904 = drop unmapped(4) + secondary(0x100) + supplementary(0x800) = primary mapped
    n = subprocess.run(["samtools", "view", "-c", "-F", "0x904", "-"],
                       stdin=mm.stdout, capture_output=True, text=True).stdout
    mm.wait()
    total = count_reads(fastq)
    return 100.0 * int(n) / total if total else 0.0


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--samples", default="samples.tsv")
    ap.add_argument("--reads-dir", default="reads")
    ap.add_argument("--out", default="analysis/cross_mapping_matrix.png")
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()

    with open(args.samples) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    bc_labels = [r["strain"] for r in rows]
    ref_labels = [r["strain"] for r in rows]
    refs = [r["reference"] for r in rows]
    fastqs = [os.path.join(args.reads_dir, r["barcode"] + ".fastq") for r in rows]

    M = np.zeros((len(rows), len(rows)))
    print("mapping each barcode against each reference ...")
    for i, fq in enumerate(fastqs):
        for j, ref in enumerate(refs):
            M[i, j] = pct_mapped(fq, ref, args.threads)
        best = ref_labels[int(M[i].argmax())]
        print(f"  {bc_labels[i]:35s} -> {best}   ({'  '.join(f'{x:4.1f}' for x in M[i])})")

    fig, ax = plt.subplots(figsize=(7.5, 6))
    im = ax.imshow(M, cmap="viridis", vmin=0, vmax=100)
    ax.set_xticks(range(len(ref_labels)))
    ax.set_xticklabels(ref_labels, rotation=35, ha="right", fontsize=8)
    ax.set_yticks(range(len(bc_labels)))
    ax.set_yticklabels(bc_labels, fontsize=8)
    ax.set_xlabel("mapped against reference")
    ax.set_ylabel("reads from barcode")
    ax.set_title("Cross-mapping: % of each barcode's reads that map")
    for i in range(len(rows)):
        for j in range(len(rows)):
            ax.text(j, i, f"{M[i, j]:.0f}", ha="center", va="center",
                    color="white" if M[i, j] < 55 else "black", fontsize=9)
    fig.colorbar(im, label="% reads mapped", shrink=0.8)
    fig.tight_layout()
    fig.savefig(args.out, dpi=130)
    print(f"-> {args.out}")


if __name__ == "__main__":
    main()
