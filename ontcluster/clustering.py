from __future__ import annotations

import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, Optional

from Bio.Align import PairwiseAligner
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .pileup import DenoisePileup, MakePileup, Pileup, ReadPileupIntersect, _default_aligner


@dataclass
class ClusterResult:
    total_reads: int
    unaligned_reads: int
    variation_count: int
    clusters: Dict[int, Seq]
    cluster_counts: Dict[int, int] = field(default_factory=dict)
    output_files: Dict[int, Path] = field(default_factory=dict)
    consensus_fasta: Optional[Path] = None

    def summary_lines(self) -> Iterable[str]:
        yield f"Processed reads: {self.total_reads}"
        yield f"Unaligned reads: {self.unaligned_reads}"
        yield f"Denoised variations: {self.variation_count}"
        yield f"Total clusters: {len(self.clusters)}"
        for cluster_id in sorted(self.cluster_counts):
            count = self.cluster_counts[cluster_id]
            filename = self.output_files.get(cluster_id)
            yield f"Cluster {cluster_id}: {count} reads -> {filename}"
        if self.consensus_fasta:
            yield f"Consensus FASTA: {self.consensus_fasta}"


def _resolve_reference(reference: str | Path | SeqRecord | Seq) -> Seq:
    if isinstance(reference, SeqRecord):
        return reference.seq
    if isinstance(reference, Seq):
        return reference
    path = Path(reference)
    if path.exists():
        return SeqIO.read(path, "fasta").seq
    return Seq(str(reference))


def _open_fastq(path: Path) -> Iterable[SeqRecord]:
    return SeqIO.parse(str(path), "fastq")


def _fastq_records(path: Path) -> Iterable[SeqRecord]:
    for record in SeqIO.parse(str(path), "fastq"):
        seq = str(record.seq).upper()
        if len(seq) < 50:
            continue
        if seq.count("N") > len(seq) / 2:
            continue
        yield record


def cluster_fastq(
    reference: str | Path | SeqRecord | Seq,
    fastq_path: str | Path,
    threshold: float,
    max_clusters: int = 100,
    output_dir: Optional[str | Path] = None,
    consensus_fasta: Optional[str | Path] = None,
    aligner: Optional[PairwiseAligner] = None,
) -> ClusterResult:
    aligner = aligner or _default_aligner()
    ref_seq = _resolve_reference(reference)
    fastq_path = Path(fastq_path)
    out_dir = Path(output_dir) if output_dir else fastq_path.parent
    out_dir.mkdir(parents=True, exist_ok=True)
    consensus_path: Optional[Path] = None
    consensus_handle = None
    if consensus_fasta:
        consensus_path = Path(consensus_fasta)
        if not consensus_path.is_absolute():
            consensus_path = out_dir / consensus_path
        consensus_path.parent.mkdir(parents=True, exist_ok=True)
        consensus_handle = open(consensus_path, "w")

    try:
        # Pass 1: build and denoise pileup
        pileup = MakePileup(ref_seq, _fastq_records(fastq_path), aligner=aligner)
        denoised = DenoisePileup(pileup, threshold)
        variation_count = denoised.variation_count()
        print(f"Processed reads: {pileup.total_reads}", file=sys.stderr)
        print(f"Unaligned reads: {pileup.unaligned_reads}", file=sys.stderr)
        print(f"Denoised variations: {variation_count}", file=sys.stderr)

        # Pass 2: discover clusters
        cluster_map: Dict[str, int] = {}
        clusters: Dict[int, Seq] = {}
        for record in _fastq_records(fastq_path):
            if len(record.seq) < len(ref_seq):
                continue
            consensus = ReadPileupIntersect(denoised, record, aligner=aligner)
            if consensus is None:
                continue
            key = str(consensus)
            if key not in cluster_map:
                if len(cluster_map) >= max_clusters:
                    raise ValueError(
                        f"Maximum cluster count {max_clusters} exceeded while processing {record.id}"
                    )
                cluster_id = len(cluster_map) + 1
                cluster_map[key] = cluster_id
                clusters[cluster_id] = consensus
                if consensus_handle:
                    consensus_handle.write(f">cluster{cluster_id}\n{key}\n")
        print(f"Total clusters: {len(clusters)}", file=sys.stderr)

        if not clusters:
            return ClusterResult(
                total_reads=pileup.total_reads,
                unaligned_reads=pileup.unaligned_reads,
                variation_count=variation_count,
                clusters={},
                consensus_fasta=consensus_path,
            )

        base_name = fastq_path.stem
        output_files: Dict[int, Path] = {}
        handles: Dict[int, any] = {}
        try:
            for cluster_id in clusters:
                filename = out_dir / f"{base_name}-C{cluster_id}.fastq"
                output_files[cluster_id] = filename
                handles[cluster_id] = open(filename, "w")

            # Pass 3: assign reads to clusters by score-only alignment
            cluster_counts: Dict[int, int] = {c: 0 for c in clusters}
            for record in _fastq_records(fastq_path):
                best_cluster = None
                best_score = float("-inf")
                read_seq = str(record.seq)
                for seq_str, cluster_id in cluster_map.items():
                    score = aligner.score(seq_str, read_seq)
                    if score > best_score:
                        best_score = score
                        best_cluster = cluster_id
                if best_cluster is None:
                    continue
                SeqIO.write(record, handles[best_cluster], "fastq")
                cluster_counts[best_cluster] += 1
        finally:
            for handle in handles.values():
                handle.close()

        for cluster_id, path in output_files.items():
            count = cluster_counts.get(cluster_id, 0)
            print(f"Cluster {cluster_id}: {count} reads -> {path}", file=sys.stderr)

        return ClusterResult(
            total_reads=pileup.total_reads,
            unaligned_reads=pileup.unaligned_reads,
            variation_count=variation_count,
            clusters=clusters,
            cluster_counts=cluster_counts,
            output_files=output_files,
            consensus_fasta=consensus_path,
        )
    finally:
        if consensus_handle:
            consensus_handle.close()
