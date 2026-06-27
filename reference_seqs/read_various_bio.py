#!/usr/bin/env python3

from Bio import SeqIO
import argparse
import sys

# option flag picks the format; choices=[...] rejects anything other than fasta/genbank
parser = argparse.ArgumentParser(description="Summarize a FASTA or GenBank file.")
parser.add_argument("filename", help="sequence file to summarize")
parser.add_argument("-f", "--format", choices=["fasta", "genbank"], required=True,
                    help="format of the file (fasta or genbank)")

# get the particular arguments that user invoked on the commandline
args = parser.parse_args()
filename = args.filename
filetype = args.format

# nucleotide counter for all sequence records parsed below
total = 0

# parse up front so that file or option format choices will fail cleanly
try:
    records = list(SeqIO.parse(filename, filetype))
except ValueError:
    sys.exit(f"Error: could not read '{filename}' as {filetype}. Is it really {filetype}?")
if not records:
    sys.exit(f"Error: no {filetype} records in '{filename}'. Is it really {filetype}?")

# only report the format once we know that the file was parsed correctly
print(f"This file is in {filetype} format")

print("---------------------------")

# read all records in file to obtain ID and SEQ per record
for seq_record in records:
    print(seq_record.id + ", length " + str(len(seq_record.seq)) + " nt")
    total = total + len(seq_record.seq)

print("---------------------------")
print("Total length parsed: " + str(total) + " nt")
