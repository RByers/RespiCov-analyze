# OntCluster

OntCluster provides streaming utilities for piling up Oxford Nanopore amplicon
reads against a reference, denoising the resulting pileup, and clustering reads
into consensus-backed groups. The library is designed for use inside Jupyter
notebooks as well as from the command line via `python -m ontcluster`.

## Design Overview

- **Pileup construction** – `MakePileup` aligns each read to the reference with
  Biopython's `PairwiseAligner`, walking the alignment to accumulate counts for
  base substitutions, deletions, and insertion sequences at each reference
  position. Reads that fail to cover the full reference are counted as
  unaligned and ignored.
- **Noise suppression** – `DenoisePileup` applies a frequency threshold per
  reference position so downstream processing retains only confident variation
  signals.
- **Read filtering** – `ReadPileupIntersect` projects individual reads onto the
  denoised pileup, discarding variations not seen in the denoised consensus to
  form per-read consensus sequences.
- **Clustering workflow** – `cluster_fastq` streams the input FASTQ three
  times: first to build and denoise the pileup, second to discover unique
  consensus sequences, and third to assign each read to the highest-scoring
  cluster while writing cluster-specific FASTQ files.

All FASTQ processing is performed in a streaming fashion so large read sets do
not need to fit in memory.

## Requirements

- Python 3.9 or newer.
- Biopython 1.85 or newer (the code relies on the modern PairwiseAligner
  `coordinates` interface).
- Pytest (optional, for running the included unit and integration tests).

## Running the Tests

```bash
PYTHONPATH=. pytest ontcluster/tests
```

## Command-Line Usage

```bash
python -m ontcluster reference.fasta reads.fastq \
  --threshold 0.05 --max-clusters 200 --consensus-fasta cluster-consensus.fa
```

Cluster FASTQ files will be written alongside the input FASTQ (or inside
`--output-dir` if supplied). When `--consensus-fasta` is provided, each
cluster's consensus sequence is also written to the specified FASTA file.
The tool prints a short summary to standard output.
