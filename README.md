# RespiCoV-analyze

I am a [molecular biology hobbyist](https://lab.rbyers.ca). This repo started as some Python [scripts](RCUtils.py) and [Jupyter notebooks](RCMatchPrimers-RC1.ipynb) for analyzing my Oxford Nanopore sequencing results following the protocol in [RespiCoV: Simultaneous identification of Severe Acute Respiratory Syndrome Coronavirus 2 (SARS-CoV-2) and 46 respiratory tract viruses and bacteria by amplicon-based Oxford-Nanopore MinION sequencing](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0264855). It has since evolved into a somewhat random collection of all my scripts for looking at my virus sequence data. I don't expect any of this to be usefully reusable by anyone else. Probably I should be embarassed by the poor quality and dumb errors in my analysis, but I make it public anyway just in case it's ever useful to share.

Files of note:
 - [RCUtils.py](RCUtils.py): My main collection of utility functions shared between notebooks
 - [RCMatchPrimers-RC1.ipynb](RCMatchPrimers-RC1.ipynb): My analysis of my first RespiCov sequencing run, matching against the primers I used and learning to understand and ignore off-target amplification
 - [RCMatchPrimers-RC2.ipynb](RCMatchPrimers-RC2.ipynb): A modified copy of the above for my second run (on a Flongle flow cell) which did not work well.
 - [RCMatchSeq.ipynb](RCMatchSeq.ipynb): A different approach analyzing the reads according to the consensus sequences I'd generated from them (all Rhinovirus).
 - [AssaySeqBinding-HRV.ipynb](AssaySeqBinding-HRV.ipynb): Analyze how well my Rhinovirus PCR primers and probes match the sequences I've found in my samples.
 - [Q8AnalyzeHRV.ipynb](Q8AnalyzeHRV.ipynb): Dig into Rhinovirus results from my Q8 run (not RespiCov). Barcoding failures made this run especially problematic.
 - [TrimPrimers.ipynb](TrimPrimers.ipynb): A little utility for trimming sequences to remove primers and any other junk outside.
 - [myseqs](myseqs): The sequences I've generated from samples in my lab (both myself and via sequencing services)
 - [PrimalScheme](PrimalScheme): Tiled amplicon primer sequences I've generated from the [[PrimalScheme tool](PrimalScheme)](https://primalscheme.com/) tool.
 - [refseq](refseq): Downloaded reference genomes related to viruses I'm studying.

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
