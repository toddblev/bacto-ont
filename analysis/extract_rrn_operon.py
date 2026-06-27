#!/usr/bin/env python3
"""
extract_rrn_operon.py - pull ONE complete rRNA (rrn) operon out of a genome.

Reads rRNA features from a GFF3, clusters them into operons (16S-23S-5S that sit
within a few kb of each other), and writes the first complete operon as a small
~5 kb FASTA. This collapsed single-operon reference is what rrn_copy_number.py
maps all reads onto to recover operon copy number.

Run once per genome (instructor prep); the small output FASTA is shipped with the
tutorial so students don't need the genome annotation at runtime.

Usage:
    extract_rrn_operon.py REF.fna GFF3 OUT.fna [--name LABEL] [--maxgap 2000]
"""
import argparse
import subprocess
import sys


def load_rrna(gff):
    feats = []
    with open(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "rRNA":
                continue
            attr = f[8]
            prod = ""
            for kv in attr.split(";"):
                if kv.startswith("product="):
                    prod = kv[8:]
            feats.append((f[0], int(f[3]), int(f[4]), f[6], prod))
    feats.sort(key=lambda x: (x[0], x[1]))
    return feats


def first_operon(feats, maxgap):
    """First complete operon = first (+)-strand 16S, walking forward to its 5S.

    Stopping at the first 5S guarantees the span is ONE operon, never two tandem
    ones (which sit only ~1 kb apart in e.g. B. subtilis).
    """
    for i, f in enumerate(feats):
        if "16S" not in f[4] or f[3] != "+":
            continue
        chrom = f[0]
        op = [f]
        for g in feats[i + 1:]:
            if g[0] != chrom or g[1] - op[-1][2] > maxgap:
                break
            op.append(g)
            if "5S" in g[4]:
                return op
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("ref")
    ap.add_argument("gff")
    ap.add_argument("out")
    ap.add_argument("--name", default="rrn_operon")
    ap.add_argument("--maxgap", type=int, default=2000)
    args = ap.parse_args()

    feats = load_rrna(args.gff)
    if not feats:
        sys.exit(f"no rRNA features found in {args.gff}")
    # operon count = number of 16S rRNA genes (one per operon) - robust to tandem operons
    n_ops = sum(1 for f in feats if "16S" in f[4])
    operon = first_operon(feats, args.maxgap)
    if operon is None:
        sys.exit("no (+)-strand 16S->5S operon found")

    chrom = operon[0][0]
    start = min(f[1] for f in operon)
    end = max(f[2] for f in operon)
    region = f"{chrom}:{start}-{end}"
    print(f"{args.name}: {n_ops} rrn operons total in genome; "
          f"extracting first = {region} ({end - start + 1} bp)")

    seq = subprocess.run(["samtools", "faidx", args.ref, region],
                         check=True, capture_output=True, text=True).stdout
    body = "".join(seq.splitlines()[1:])
    with open(args.out, "w") as out:
        out.write(f">{args.name} {region} ({len(body)} bp, 1 of {n_ops} operons)\n")
        for i in range(0, len(body), 70):
            out.write(body[i:i + 70] + "\n")
    print(f"  -> wrote {args.out}")


if __name__ == "__main__":
    main()
