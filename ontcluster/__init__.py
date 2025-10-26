"""
OntCluster library entry point.

This package provides streaming utilities for piling up Oxford Nanopore amplicon
reads against a reference, denoising the resulting pileup, and clustering reads
based on the consensus variation set.
"""

from .pileup import (
    BaseVariation,
    DeletionVariation,
    InsertionVariation,
    Pileup,
    MakePileup,
    DenoisePileup,
    ReadPileupIntersect,
)
from .clustering import ClusterResult, cluster_fastq

__all__ = [
    "BaseVariation",
    "DeletionVariation",
    "InsertionVariation",
    "Pileup",
    "MakePileup",
    "DenoisePileup",
    "ReadPileupIntersect",
    "ClusterResult",
    "cluster_fastq",
]
