"""
Utilities for denoising Oxford Nanopore alignments.

This module exposes a streaming denoiser that reads a BAM file, evaluates
per-position depth and allele frequencies against a reference genome, and
emits filtered reads in FASTQ format. Variants that do not meet the supplied
frequency threshold are coerced back to the reference (with minimum quality),
and regions with insufficient depth are excluded from the output (splitting
reads as needed). All processing is performed in two passes over the BAM to
avoid holding the full read set in memory.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import Dict, List

from Bio import SeqIO
import bamnostic


# Map IUPAC bases to their complements (upper case only)
_COMPLEMENT_MAP = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "R": "Y",
    "Y": "R",
    "S": "S",
    "W": "W",
    "K": "M",
    "M": "K",
    "B": "V",
    "V": "B",
    "D": "H",
    "H": "D",
    "N": "N",
}

_MIN_PHRED = 0  # minimum quality assigned when copying reference sequence


@dataclass
class _CoverageTables:
    """Convenience container for aggregated depth and allele counts."""

    depth: Dict[str, Dict[int, int]]
    bases: Dict[str, Dict[int, Counter]]
    inserts: Dict[str, Dict[int, Counter]]
    deletes: Dict[str, Dict[int, Counter]]


def _complement_base(base: str) -> str:
    """Return the complement for a single nucleotide (upper-case)."""
    return _COMPLEMENT_MAP.get(base.upper(), "N")


def _reverse_complement(seq: str) -> str:
    """Return the reverse complement for a nucleotide sequence."""
    return "".join(_complement_base(base) for base in reversed(seq))


def _initialise_tables() -> _CoverageTables:
    """Initialise nested lookup tables for coverage statistics."""
    depth = defaultdict(lambda: defaultdict(int))
    bases = defaultdict(lambda: defaultdict(Counter))
    inserts = defaultdict(lambda: defaultdict(Counter))
    deletes = defaultdict(lambda: defaultdict(Counter))
    return _CoverageTables(depth=depth, bases=bases, inserts=inserts, deletes=deletes)


def _collect_depth_and_frequencies(
    bam_path: str, reference: Dict[str, str]
) -> _CoverageTables:
    """First pass: aggregate depth and allele counts across the BAM file."""
    tables = _initialise_tables()

    with bamnostic.AlignmentFile(bam_path, "rb") as bam:
        for alignment in bam:
            if alignment.is_unmapped:
                continue

            contig = alignment.reference_name
            if contig not in reference:
                raise ValueError(
                    f"Read mapped to {contig}, which is absent from reference fasta."
                )

            sequence = getattr(alignment, "query_sequence", None)
            qualities = getattr(alignment, "query_qualities", None)
            if sequence is None or qualities is None:
                continue

            sequence = sequence.upper()
            is_reverse = alignment.is_reverse
            read_index = 0
            ref_index = alignment.reference_start

            for op, length in alignment.cigartuples:
                if op in (0, 7, 8):  # alignment match/mismatch
                    for offset in range(length):
                        read_pos = read_index + offset
                        if read_pos >= len(sequence):
                            continue
                        ref_pos = ref_index + offset
                        read_base = sequence[read_pos]
                        allele_base = (
                            read_base if not is_reverse else _complement_base(read_base)
                        )
                        tables.depth[contig][ref_pos] += 1
                        tables.bases[contig][ref_pos][allele_base] += 1
                    read_index += length
                    ref_index += length

                elif op == 1:  # insertion relative to reference
                    inserted = sequence[read_index : read_index + length]
                    ins_key = inserted if not is_reverse else _reverse_complement(inserted)
                    tables.inserts[contig][ref_index][ins_key] += 1
                    read_index += length

                elif op == 2:  # deletion relative to reference
                    deleted = reference[contig][ref_index : ref_index + length]
                    tables.deletes[contig][ref_index][deleted] += 1
                    for offset in range(length):
                        tables.depth[contig][ref_index + offset] += 1
                    ref_index += length

                elif op == 3:  # N (reference skip)
                    ref_index += length

                elif op == 4:  # soft clip
                    read_index += length

                elif op in (5, 6):  # hard clip or padding
                    continue

                else:
                    raise ValueError(f"Unhandled CIGAR operation: {op}")

    return tables


def _write_fastq_record(
    handle, read_id: str, seq_bases: List[str], quals: List[int]
) -> None:
    """Write a FASTQ record given prepared bases and phred qualities."""
    if not seq_bases:
        return
    sequence = "".join(seq_bases)
    quality = "".join(chr(q + 33) for q in quals)
    handle.write(f"@{read_id}\n{sequence}\n+\n{quality}\n")


def denoise_nanopore_reads(
    bam_path: str,
    reference_fasta: str,
    output_fastq: str,
    min_variant_frequency: float,
    min_read_depth: int,
) -> None:
    """
    Denoise Oxford Nanopore reads by applying depth and allele-frequency filters.

    Parameters
    ----------
    bam_path : str
        Path to the input BAM file containing aligned reads.
    reference_fasta : str
        Reference genome in FASTA format. May include multiple contigs.
    output_fastq : str
        Destination FASTQ file for denoised reads.
    min_variant_frequency : float
        Minimum allele frequency required to retain a variant.
    min_read_depth : int
        Minimum depth required to retain aligned positions. Regions below this
        depth are excised from the output (splitting reads as needed).
    """
    if not 0.0 <= min_variant_frequency <= 1.0:
        raise ValueError("min_variant_frequency must lie between 0.0 and 1.0")
    if min_read_depth < 0:
        raise ValueError("min_read_depth must be non-negative")

    with open(reference_fasta, "r", encoding="utf-8") as handle:
        reference = {
            record.id: str(record.seq).upper()
            for record in SeqIO.parse(handle, "fasta")
        }
    if not reference:
        raise ValueError(
            f"No sequences found in reference fasta: {reference_fasta}"
        )
    tables = _collect_depth_and_frequencies(bam_path, reference)

    with bamnostic.AlignmentFile(bam_path, "rb") as bam, open(
        output_fastq, "w"
    ) as writer:
        for alignment in bam:
            if alignment.is_unmapped:
                continue

            contig = alignment.reference_name
            if contig not in reference:
                continue  # already raised once in first pass; skip defensively

            read_sequence = getattr(alignment, "query_sequence", None)
            read_qualities = getattr(alignment, "query_qualities", None)
            if read_sequence is None:
                continue

            read_sequence = read_sequence.upper()
            quality_values = (
                list(read_qualities) if read_qualities is not None else [0] * len(read_sequence)
            )
            is_reverse = alignment.is_reverse
            read_index = 0
            ref_index = alignment.reference_start

            base_id = f"{alignment.query_name}_OD"
            segment_counter = 0
            segment_seq: List[str] = []
            segment_qual: List[int] = []

            def flush_segment() -> None:
                nonlocal segment_seq, segment_qual, segment_counter
                if not segment_seq:
                    return
                suffix = "" if segment_counter == 0 else f":{segment_counter + 1}"
                record_id = f"{base_id}{suffix}"
                _write_fastq_record(writer, record_id, segment_seq, segment_qual)
                segment_seq = []
                segment_qual = []
                segment_counter += 1

            for op, length in alignment.cigartuples:
                if op in (0, 7, 8):  # alignment match/mismatch
                    for offset in range(length):
                        read_pos = read_index + offset
                        if read_pos >= len(read_sequence):
                            continue
                        ref_pos = ref_index + offset
                        base_depth = tables.depth[contig].get(ref_pos, 0)
                        if base_depth < min_read_depth:
                            flush_segment()
                            continue

                        read_base = read_sequence[read_pos]
                        phred = quality_values[read_pos]
                        ref_base = reference[contig][ref_pos]
                        allele_base = (
                            read_base if not is_reverse else _complement_base(read_base)
                        )
                        allele_counts = tables.bases[contig].get(ref_pos)
                        allele_support = (
                            allele_counts.get(allele_base, 0) if allele_counts else 0
                        )
                        is_reference = allele_base == ref_base
                        if is_reference:
                            keep_base = True
                        else:
                            freq = allele_support / base_depth if base_depth else 0.0
                            keep_base = freq >= min_variant_frequency

                        if keep_base:
                            out_base = read_base
                            out_phred = phred
                        else:
                            ref_oriented = (
                                ref_base if not is_reverse else _complement_base(ref_base)
                            )
                            out_base = ref_oriented
                            out_phred = _MIN_PHRED

                        segment_seq.append(out_base)
                        segment_qual.append(out_phred)

                    read_index += length
                    ref_index += length

                elif op == 1:  # insertion
                    inserted_seq = read_sequence[read_index : read_index + length]
                    inserted_qual = quality_values[read_index : read_index + length]
                    anchor_depth = tables.depth[contig].get(
                        ref_index, tables.depth[contig].get(ref_index - 1, 0)
                    )
                    if anchor_depth < min_read_depth:
                        flush_segment()
                    else:
                        ins_key = (
                            inserted_seq
                            if not is_reverse
                            else _reverse_complement(inserted_seq)
                        )
                        ins_counter = tables.inserts[contig].get(ref_index)
                        support = ins_counter.get(ins_key, 0) if ins_counter else 0
                        freq = support / anchor_depth if anchor_depth else 0.0
                        if freq >= min_variant_frequency:
                            segment_seq.extend(inserted_seq)
                            segment_qual.extend(inserted_qual)
                    read_index += length

                elif op == 2:  # deletion
                    deleted_seq = reference[contig][ref_index : ref_index + length]
                    depth_span = [
                        tables.depth[contig].get(ref_index + offset, 0)
                        for offset in range(length)
                    ]
                    if any(depth < min_read_depth for depth in depth_span):
                        flush_segment()
                    else:
                        deletion_counter = tables.deletes[contig].get(ref_index)
                        support = (
                            deletion_counter.get(deleted_seq, 0)
                            if deletion_counter
                            else 0
                        )
                        depth_at_start = tables.depth[contig].get(ref_index, 0)
                        freq = support / depth_at_start if depth_at_start else 0.0
                        if freq < min_variant_frequency:
                            restored = (
                                deleted_seq
                                if not is_reverse
                                else _reverse_complement(deleted_seq)
                            )
                            segment_seq.extend(restored)
                            segment_qual.extend([_MIN_PHRED] * len(restored))
                    ref_index += length

                elif op == 3:  # N (skip)
                    ref_index += length

                elif op == 4:  # soft clip
                    read_index += length

                elif op in (5, 6):  # hard clip / padding
                    continue

                else:
                    raise ValueError(f"Unhandled CIGAR operation: {op}")

            flush_segment()


__all__ = ["denoise_nanopore_reads"]
