from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .clustering import cluster_fastq


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Cluster ONT amplicon reads by pileup consensus.")
    parser.add_argument("reference", help="Reference sequence (FASTA file or literal sequence).")
    parser.add_argument("fastq", help="Input FASTQ file containing ONT reads.")
    parser.add_argument(
        "-t",
        "--threshold",
        type=float,
        default=0.05,
        help="Noise threshold for denoising pileup (default: 0.05).",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        default=None,
        help="Directory for cluster FASTQ files (defaults to FASTQ directory).",
    )
    parser.add_argument(
        "-m",
        "--max-clusters",
        type=int,
        default=100,
        help="Maximum number of clusters to create (default: 100).",
    )
    parser.add_argument(
        "-c",
        "--consensus-fasta",
        type=Path,
        default=None,
        help="Optional FASTA file to write cluster consensus sequences.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    result = cluster_fastq(
        reference=args.reference,
        fastq_path=args.fastq,
        threshold=args.threshold,
        max_clusters=args.max_clusters,
        output_dir=args.output_dir,
        consensus_fasta=args.consensus_fasta,
    )
    for line in result.summary_lines():
        print(line)
    return 0


if __name__ == "__main__":
    sys.exit(main())
