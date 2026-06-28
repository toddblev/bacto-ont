#!/usr/bin/env python3
"""
coverage_track.py - windowed read-depth track along a genome, from a sorted BAM.

Shows how ONT read depth varies across a bacterial chromosome. Two teaching points:

  * The genome is sampled fairly evenly (depth ~ flat around the median), so the
    *area* under the track is the strain's sequencing yield.
  * The multi-copy rRNA (rrn) operons do NOT appear as tall spikes here: ONT reads
    plus `minimap2 --secondary=no` spread the operon reads across the identical
    copies, so each operon locus sits near the median. Copy number is recovered
    separately with rrn_copy_number.py (the "collapsed reference" trick).

An optional MAPQ-filtered overlay (--mapq 20) removes ambiguous (MAPQ 0) reads;
the rrn loci then turn into *dips*, which visually proves the multi-mapping.

Usage:
    coverage_track.py BAM OUT.png [--window 1000] [--title "..."] [--mapq 20]
                      [--marks bp,bp,...] [--operon OP.fna --ref GENOME.fna]

BAM must be coordinate-sorted and indexed (samtools sort + samtools index).

Pass --operon (a single extracted rrn operon) together with --ref (the genome it
was mapped to) to draw a vertical line at every genomic copy: the operon is aligned
back to the genome with minimap2, which recovers all (near-identical) copies. Or give
explicit positions with --marks.
"""
import argparse
import gzip
import os
import subprocess
import tempfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def mosdepth_windows(bam, window, mapq=0):
    """Run mosdepth and return {chrom: (centers, depths)} as numpy arrays."""
    tmp = tempfile.mkdtemp(prefix="covtrack_")
    pre = os.path.join(tmp, "md")
    cmd = ["mosdepth", "--by", str(window), "--no-per-base", "--fast-mode", "-t", "4"]
    if mapq:
        cmd += ["-Q", str(mapq)]
    cmd += [pre, bam]
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    chroms = {}
    with gzip.open(pre + ".regions.bed.gz", "rt") as fh:
        for line in fh:
            c, s, e, d = line.rstrip("\n").split("\t")
            chroms.setdefault(c, []).append(((int(s) + int(e)) / 2.0, float(d)))
    out = {}
    for c, rows in chroms.items():
        arr = np.array(rows)
        out[c] = (arr[:, 0], arr[:, 1])
    return out


def operon_positions(operon_fa, ref_fa, chrom, min_frac=0.5):
    """Locate every copy of a single operon on `chrom`.

    The operon is ~identical across its genomic copies, so aligning it back to the
    genome (allowing secondary hits) lands one alignment per copy. Returns the copy
    midpoints (bp), sorted; overlapping alignments of one copy are merged."""
    paf = subprocess.run(
        ["minimap2", "-x", "asm5", "-N", "60", "-p", "0.05", "--secondary=yes",
         ref_fa, operon_fa],
        capture_output=True, text=True).stdout
    mids = []
    for line in paf.splitlines():
        f = line.split("\t")
        if len(f) < 11 or f[5] != chrom:
            continue
        qlen, tstart, tend, alnlen = int(f[1]), int(f[7]), int(f[8]), int(f[10])
        if alnlen < min_frac * qlen:           # skip short partial hits
            continue
        mids.append((tstart + tend) / 2.0)
    mids.sort()
    merged = []
    for m in mids:
        if not merged or m - merged[-1] >= 3000:   # > operon length apart = distinct copy
            merged.append(m)
    return merged


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("bam")
    ap.add_argument("out_png")
    ap.add_argument("--window", type=int, default=1000, help="window size in bp (default 1000)")
    ap.add_argument("--title", default=None)
    ap.add_argument("--mapq", type=int, default=0,
                    help="if >0, also draw a MAPQ-filtered overlay at this threshold (e.g. 20)")
    ap.add_argument("--marks", default="",
                    help="comma-separated bp positions to mark (e.g. rrn operon locations)")
    ap.add_argument("--operon", default=None,
                    help="single operon FASTA; its genomic copies are located and marked")
    ap.add_argument("--ref", default=None,
                    help="genome FASTA the operon is aligned back to (needed with --operon)")
    args = ap.parse_args()

    tracks = mosdepth_windows(args.bam, args.window, mapq=0)
    # pick the longest contig (the chromosome), ignore small plasmids for the main track
    chrom = max(tracks, key=lambda c: len(tracks[c][0]))
    x, y = tracks[chrom]
    median = float(np.median(y))

    marks = [float(p) for p in args.marks.split(",") if p.strip()]
    if args.operon and args.ref:
        marks = operon_positions(args.operon, args.ref, chrom)
        print(f"  located {len(marks)} {os.path.basename(args.operon)} copies on {chrom}")

    plt.figure(figsize=(12, 3.2))
    plt.plot(x / 1e6, y, lw=0.6, color="#1f77b4", label="all reads")
    plt.axhline(median, color="k", lw=0.8, ls="--", label=f"median {median:.1f}x")

    if args.mapq:
        ftracks = mosdepth_windows(args.bam, args.window, mapq=args.mapq)
        if chrom in ftracks:
            fx, fy = ftracks[chrom]
            plt.plot(fx / 1e6, fy, lw=0.6, color="#d62728",
                     label=f"MAPQ>={args.mapq} only (rrn -> dips)")

    ymax = max(median * 4, float(np.percentile(y, 99.5)))

    # mark operon loci with light vertical lines
    if marks:
        for k, m in enumerate(marks):
            plt.axvline(m / 1e6, color="0.55", lw=0.8, ls=":",
                        label=f"rrn operons (n={len(marks)})" if k == 0 else None)

    plt.xlim(x[0] / 1e6, x[-1] / 1e6)
    plt.ylim(0, ymax)
    plt.xlabel(f"position on {chrom} (Mb)")
    plt.ylabel(f"depth ({args.window} bp windows)")
    plt.title(args.title or f"Read depth across {chrom}")
    plt.legend(loc="upper right", fontsize=8)
    plt.tight_layout()
    plt.savefig(args.out_png, dpi=130)
    print(f"{chrom}: {len(x)} windows, median depth {median:.1f}x  ->  {args.out_png}")


if __name__ == "__main__":
    main()
