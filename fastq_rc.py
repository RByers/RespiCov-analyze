#!/usr/bin/env python3
# Reverse complement a FASTQ file in place, preserving quality scores
from Bio import SeqIO
import sys
import os

fastq_file = sys.argv[1]
records = []
for rec in SeqIO.parse(fastq_file, "fastq"):
    rec.seq = rec.seq.reverse_complement()
    rec.letter_annotations["phred_quality"] = rec.letter_annotations["phred_quality"][::-1]
    records.append(rec)
SeqIO.write(records, fastq_file, "fastq")