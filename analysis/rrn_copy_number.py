#!/usr/bin/env python3
"""
rrn_copy_number.py - estimate rRNA (rrn) operon copy number from ONT coverage depth.

The "collapsed reference" trick:
  1. genome-median depth  = the single-copy depth, measured from the full-genome BAM.
  2. map ALL reads to ONE extracted rrn operon  ->  reads from every genomic copy
     pile onto the single locus.
  3. copy number  ~=  operon pile-up depth / genome-median depth.

Why not just read it off the genome? Long ONT reads + `minimap2 --secondary=no`
spread the operon reads across the identical copies, so on the full genome each
operon sits near the median (no tall spike). Collapsing to one operon recovers the
count. The estimate slightly UNDER-counts when reads are shorter than the ~5 kb
operon (partial coverage + edge effects) - a good discussion point.

Usage:
  rrn_copy_number.py --fastq R.fastq --operon operon.fna --full-bam full.sorted.bam
                     [--label "B. subtilis"] [--expected 10] [--threads 8]

`operon.fna` is a single ~5 kb rrn operon pre-extracted from the reference
(see extract_rrn_operon.py). `full.sorted.bam` is the reads mapped to the whole
genome (coordinate-sorted + indexed).
"""
import argparse
import gzip
import os
import subprocess
import tempfile

import numpy as np


def genome_median_depth(full_bam, window=1000):
    """Median windowed depth on the longest contig (single-copy proxy)."""
    tmp = tempfile.mkdtemp(prefix="rrn_gm_")
    pre = os.path.join(tmp, "md")
    subprocess.run(["mosdepth", "--by", str(window), "--no-per-base", "--fast-mode",
                    "-t", "4", pre, full_bam], check=True, capture_output=True, text=True)
    by_chrom = {}
    with gzip.open(pre + ".regions.bed.gz", "rt") as fh:
        for line in fh:
            c, s, e, d = line.rstrip("\n").split("\t")
            by_chrom.setdefault(c, []).append(float(d))
    chrom = max(by_chrom, key=lambda c: len(by_chrom[c]))
    return float(np.median(by_chrom[chrom])), chrom


def operon_depth(fastq, operon_fna, threads=8, edge=400):
    """Map reads to the single operon; return mean depth over its central region."""
    tmp = tempfile.mkdtemp(prefix="rrn_op_")
    bam = os.path.join(tmp, "op.bam")
    # map reads to the single operon (one place to go -> all copies pile up)
    p1 = subprocess.Popen(["minimap2", "-ax", "map-ont", "--secondary=no",
                           "-t", str(threads), operon_fna, fastq],
                          stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    subprocess.run(["samtools", "sort", "-o", bam, "-"], stdin=p1.stdout,
                   check=True, capture_output=True)
    p1.wait()
    subprocess.run(["samtools", "index", bam], check=True, capture_output=True)
    # operon length
    faidx = subprocess.run(["samtools", "faidx", operon_fna], capture_output=True, text=True)
    with open(operon_fna + ".fai") as fh:
        name, length = fh.readline().split("\t")[:2]
    length = int(length)
    region = f"{name}:{edge}-{length - edge}"
    dep = subprocess.run(["samtools", "depth", "-a", "-r", region, bam],
                         capture_output=True, text=True).stdout
    vals = [int(l.split("\t")[2]) for l in dep.splitlines()]
    return (float(np.mean(vals)) if vals else 0.0), length


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fastq", required=True)
    ap.add_argument("--operon", required=True, help="single rrn operon FASTA (~5 kb)")
    ap.add_argument("--full-bam", required=True, help="reads mapped to whole genome, sorted+indexed")
    ap.add_argument("--label", default="")
    ap.add_argument("--expected", type=int, default=None, help="known copy number, for comparison")
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()

    gmed, chrom = genome_median_depth(args.full_bam)
    odep, olen = operon_depth(args.fastq, args.operon, threads=args.threads)
    cn = odep / gmed if gmed else float("nan")

    lab = args.label or os.path.basename(args.fastq)
    print(f"== rrn operon copy number: {lab} ==")
    print(f"  genome-median depth (single copy) : {gmed:5.1f}x   ({chrom})")
    print(f"  single-operon pile-up depth       : {odep:5.1f}x   (operon {olen} bp)")
    print(f"  estimated copy number  = {odep:.1f} / {gmed:.1f} = {cn:.1f}"
          + (f"   (expected {args.expected})" if args.expected else ""))
    return cn


if __name__ == "__main__":
    main()
