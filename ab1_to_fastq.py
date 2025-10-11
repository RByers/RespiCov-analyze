#!/usr/bin/env python3
"""
ab1_to_fastq.py — convert ABI/AB1 chromatograms to FASTQ using Biopython.

Usage:
  ab1_to_fastq.py FILE1.ab1 [FILE2.ab1 ...]
Notes:
  - Writes one FASTQ per input, same basename with .fastq extension.
  - Preserves ABI Phred qualities.
  - Accepts .ab1 or .abi.
"""

import os
import sys
from Bio import SeqIO

def usage() -> None:
    print(__doc__.strip(), file=sys.stderr)

def out_path(in_path: str) -> str:
    base, _ = os.path.splitext(in_path)
    return base + ".fastq"

def convert_one(in_path: str) -> int:
    try:
        out = out_path(in_path)
        SeqIO.convert(in_path, "abi", out, "fastq")
        print(f"Wrote {out}")
        return 0
    except Exception as e:
        print(f"Error converting {in_path}: {e}", file=sys.stderr)
        return 1

def main(argv) -> int:
    if not argv:
        usage()
        return 2
    rc = 0
    for p in argv:
        rc |= convert_one(p)
    return rc

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
