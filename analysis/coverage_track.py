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

BAM must be coordinate-sorted and indexed (samtools sort + samtools index).
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


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("bam")
    ap.add_argument("out_png")
    ap.add_argument("--window", type=int, default=1000, help="window size in bp (default 1000)")
    ap.add_argument("--title", default=None)
    ap.add_argument("--mapq", type=int, default=0,
                    help="if >0, also draw a MAPQ-filtered overlay at this threshold (e.g. 20)")
    args = ap.parse_args()

    tracks = mosdepth_windows(args.bam, args.window, mapq=0)
    # pick the longest contig (the chromosome), ignore small plasmids for the main track
    chrom = max(tracks, key=lambda c: len(tracks[c][0]))
    x, y = tracks[chrom]
    median = float(np.median(y))

    plt.figure(figsize=(12, 3.2))
    plt.plot(x / 1e6, y, lw=0.6, color="#1f77b4", label="all reads")
    plt.axhline(median, color="k", lw=0.8, ls="--", label=f"median {median:.1f}x")

    if args.mapq:
        ftracks = mosdepth_windows(args.bam, args.window, mapq=args.mapq)
        if chrom in ftracks:
            fx, fy = ftracks[chrom]
            plt.plot(fx / 1e6, fy, lw=0.6, color="#d62728",
                     label=f"MAPQ>={args.mapq} only (rrn -> dips)")

    plt.xlim(x[0] / 1e6, x[-1] / 1e6)
    plt.ylim(0, max(median * 4, float(np.percentile(y, 99.5))))
    plt.xlabel(f"position on {chrom} (Mb)")
    plt.ylabel(f"depth ({args.window} bp windows)")
    plt.title(args.title or f"Read depth across {chrom}")
    plt.legend(loc="upper right", fontsize=8)
    plt.tight_layout()
    plt.savefig(args.out_png, dpi=130)
    print(f"{chrom}: {len(x)} windows, median depth {median:.1f}x  ->  {args.out_png}")


if __name__ == "__main__":
    main()
