#!/usr/bin/env python3

from Bio import SeqIO
import sys

# assign command-line argument to variable "filename"
filename = sys.argv[1]

# nucleotide counter for all genbank records parsed below
total = 0

# parse up front so the wrong format fails here with a clear message
try:
    records = list(SeqIO.parse(filename, "genbank"))
except ValueError:
    sys.exit(f"Error: could not read '{filename}' as genbank. Is it really genbank?")
if not records:
    sys.exit(f"Error: no genbank records in '{filename}'. Is it really genbank?")

# read all genbank-formatted records in file to obtain ID and SEQ per record
for seq_record in records:
    print(seq_record.id + ", length " + str(len(seq_record.seq)) + " nt")
    total = total + len(seq_record.seq)

print("---------------------------")
print("Total length parsed: " + str(total) + " nt")
