from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Tuple, Union

from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


@dataclass(frozen=True)
class BaseVariation:
    base: str

    _VALID_BASES = {"A", "C", "G", "T", "N"}

    def __post_init__(self) -> None:
        if len(self.base) != 1:
            raise ValueError("base variation requires a single character value")
        base = self.base.upper()
        if base not in self._VALID_BASES:
            raise ValueError(f"unrecognized base: {self.base}")
        object.__setattr__(self, "base", base)


@dataclass(frozen=True)
class DeletionVariation:
    pass


@dataclass(frozen=True)
class InsertionVariation:
    sequence: str

    def __post_init__(self) -> None:
        if not self.sequence:
            raise ValueError("insertion variation requires a sequence value")
        object.__setattr__(self, "sequence", self.sequence.upper())


VariationType = Union[BaseVariation, DeletionVariation, InsertionVariation]


class Pileup:
    """Pileup representation of variations against a reference sequence."""

    def __init__(self, reference: Seq):
        if not isinstance(reference, Seq):
            raise TypeError("reference must be a Bio.Seq.Seq instance")
        self.reference = reference
        self.variation_counts: List[Dict[VariationType, int]] = [
            {} for _ in range(len(self.reference))
        ]
        self.total_reads: int = 0
        self.unaligned_reads: int = 0

    def __len__(self) -> int:
        return len(self.reference)

    def _ensure_index(self, index: int) -> None:
        if index < 0 or index >= len(self.reference):
            raise IndexError("variation index out of range")

    def _add_variation(self, index: int, variation: VariationType, count: int = 1) -> None:
        self._ensure_index(index)
        slot = self.variation_counts[index]
        slot[variation] = slot.get(variation, 0) + count

    def increment_base(self, index: int, base: str) -> None:
        variation = BaseVariation(base)
        self._add_variation(index, variation)

    def increment_deletion(self, index: int) -> None:
        variation = DeletionVariation()
        self._add_variation(index, variation)

    def increment_insertion(self, index: int, sequence: str) -> None:
        variation = InsertionVariation(sequence)
        self._add_variation(index, variation)

    def variation_count(self) -> int:
        return sum(len(slot) for slot in self.variation_counts)

    def filtered(self, threshold: float) -> "Pileup":
        filtered = Pileup(self.reference)
        filtered.total_reads = self.total_reads
        filtered.unaligned_reads = self.unaligned_reads
        for index, slot in enumerate(self.variation_counts):
            for variation, count in slot.items():
                # Keep variations that meet the frequency threshold.
                # Note that we don't take the total variations at this position as that would
                # double-count insertions.
                if count / self.total_reads >= threshold:
                    filtered.variation_counts[index][variation] = count
        return filtered

    def get_variation_count(self, index: int, variation: VariationType) -> int:
        self._ensure_index(index)
        return self.variation_counts[index].get(variation, 0)

    def _has_variation(self, index: int, variation: VariationType) -> bool:
        return self.get_variation_count(index, variation) > 0

    def base_allowed(self, index: int, base: str) -> bool:
        base = base.upper()
        if index >= len(self.reference):
            return False
        if base == self.reference[index].upper():
            return True
        return self._has_variation(index, BaseVariation(base))

    def deletion_allowed(self, index: int) -> bool:
        if index >= len(self.reference):
            return False
        return self._has_variation(index, DeletionVariation())

    def insertion_allowed(self, index: int, sequence: str) -> bool:
        if index >= len(self.reference):
            return False
        return self._has_variation(index, InsertionVariation(sequence))

    def to_debug_string(self) -> str:
        lines: List[str] = [
            f"Pileup(total_reads={self.total_reads}, unaligned_reads={self.unaligned_reads})"
        ]
        for index, slot in enumerate(self.variation_counts):
            if not slot:
                continue
            entry = ", ".join(
                f"{_variation_label(variation)}={count}"
                for variation, count in sorted(
                    slot.items(), key=lambda item: _variation_sort_key(item[0])
                )
            )
            lines.append(f"{index}: {entry}")
        return "\n".join(lines)


def _variation_label(variation: VariationType) -> str:
    if isinstance(variation, BaseVariation):
        return f"base:{variation.base}"
    if isinstance(variation, DeletionVariation):
        return "deletion"
    if isinstance(variation, InsertionVariation):
        return f"ins:{variation.sequence}"
    return "unknown"


def _variation_sort_key(variation: VariationType) -> Tuple[int, str]:
    if isinstance(variation, BaseVariation):
        return (0, variation.base)
    if isinstance(variation, DeletionVariation):
        return (1, "")
    if isinstance(variation, InsertionVariation):
        return (2, variation.sequence)
    return (3, "")


def _default_aligner() -> PairwiseAligner:
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -2.0
    aligner.extend_gap_score = -0.5
    return aligner


def _validate_sequence(sequence: str) -> None:
    for char in sequence.upper():
        if char not in BaseVariation._VALID_BASES:
            raise ValueError(f"unrecognized base: {char}")

def _alignment_coordinates(alignment):
    return alignment.coordinates


def _insertion_anchor(ref_position: int, ref_len: int) -> int:
    if ref_len <= 0:
        return 0
    if ref_position >= ref_len:
        return ref_len - 1
    return ref_position


def _process_alignment_with_coordinates(
    pileup: Pileup, read_seq: str, coords
) -> None:
    if coords is None:
        raise ValueError("Alignment coordinates missing")
    read_seq = read_seq.upper()
    ref_len = len(pileup)
    num_steps = coords.shape[1] - 1
    for idx in range(num_steps):
        ref_start = int(coords[0, idx])
        ref_end = int(coords[0, idx + 1])
        read_start = int(coords[1, idx])
        read_end = int(coords[1, idx + 1])
        ref_delta = ref_end - ref_start
        read_delta = read_end - read_start
        if ref_delta > 0 and read_delta > 0:
            span = ref_delta
            for offset in range(span):
                ref_index = ref_start + offset
                read_index = read_start + offset
                pileup.increment_base(ref_index, read_seq[read_index])
        elif ref_delta > 0:
            for ref_index in range(ref_start, ref_end):
                pileup.increment_deletion(ref_index)
        elif read_delta > 0:
            insertion_seq = read_seq[read_start:read_end]
            if insertion_seq:
                anchor = _insertion_anchor(ref_start, ref_len)
                pileup.increment_insertion(anchor, insertion_seq)
        # If both deltas zero, nothing to do


def MakePileup(
    reference: Seq | SeqRecord | str,
    reads: Iterable[SeqRecord],
    aligner: Optional[PairwiseAligner] = None,
) -> Pileup:
    if isinstance(reference, SeqRecord):
        ref_seq = reference.seq
    else:
        ref_seq = Seq(str(reference))
    pileup = Pileup(ref_seq)
    aligner = aligner or _default_aligner()

    for read in reads:
        pileup.total_reads += 1
        read_seq = str(read.seq)
        _validate_sequence(read_seq)
        alignment = aligner.align(str(pileup.reference), read_seq)
        if not alignment or alignment.score <= 0:
            pileup.unaligned_reads += 1
            continue
        aligned = alignment[0]
        coords = _alignment_coordinates(aligned)
        if coords is None or coords.shape[1] == 0:
            pileup.unaligned_reads += 1
            continue
        if int(coords[0, 0]) != 0 or int(coords[0, -1]) != len(pileup):
            pileup.unaligned_reads += 1
            continue
        _process_alignment_with_coordinates(pileup, read_seq, coords)
    return pileup


def DenoisePileup(pileup: Pileup, threshold: float) -> Pileup:
    if not 0 <= threshold <= 1:
        raise ValueError("threshold must be between 0 and 1")
    return pileup.filtered(threshold)


def ReadPileupIntersect(
    pileup: Pileup,
    read: SeqRecord,
    aligner: Optional[PairwiseAligner] = None,
) -> Optional[Seq]:
    aligner = aligner or _default_aligner()
    alignment = aligner.align(str(pileup.reference), str(read.seq))
    if not alignment or alignment.score <= 0:
        return None
    if len(read.seq) < len(pileup):
        return None
    aligned = alignment[0]
    coords = _alignment_coordinates(aligned)
    if coords is None or coords.shape[1] == 0:
        return None
    if int(coords[0, 0]) != 0 or int(coords[0, -1]) != len(pileup):
        return None
    return _intersect_with_coordinates(pileup, str(read.seq), coords)

def _intersect_with_coordinates(pileup: Pileup, read_seq: str, coords) -> Seq:
    read_seq = read_seq.upper()
    result: List[str] = []
    ref_len = len(pileup)
    num_steps = coords.shape[1] - 1
    for idx in range(num_steps):
        ref_start = int(coords[0, idx])
        ref_end = int(coords[0, idx + 1])
        read_start = int(coords[1, idx])
        read_end = int(coords[1, idx + 1])
        ref_delta = ref_end - ref_start
        read_delta = read_end - read_start
        if ref_delta > 0 and read_delta > 0:
            span = ref_delta
            for offset in range(span):
                ref_index = ref_start + offset
                read_base = read_seq[read_start + offset]
                if pileup.base_allowed(ref_index, read_base):
                    result.append(read_base)
                else:
                    result.append(pileup.reference[ref_index].upper())
        elif ref_delta > 0:
            for ref_index in range(ref_start, ref_end):
                if not pileup.deletion_allowed(ref_index):
                    result.append(pileup.reference[ref_index].upper())
        elif read_delta > 0:
            insertion_seq = read_seq[read_start:read_end]
            if insertion_seq:
                anchor = _insertion_anchor(ref_start, ref_len)
                if pileup.insertion_allowed(anchor, insertion_seq):
                    result.append(insertion_seq)
    return Seq("".join(result))
