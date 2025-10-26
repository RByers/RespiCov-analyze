import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from ontcluster.pileup import (
    BaseVariation,
    DeletionVariation,
    InsertionVariation,
    MakePileup,
    DenoisePileup,
    Pileup,
    ReadPileupIntersect,
)


def test_makepileup_counts_matches_and_snps():
    reference = Seq("ACGT")
    reads = [SeqRecord(Seq("ACGT")), SeqRecord(Seq("AGGT"))]
    pileup = MakePileup(reference, reads)

    assert pileup.total_reads == 2
    assert pileup.unaligned_reads == 0
    assert pileup.option_count() == 1
    assert pileup.get_variation_count(0, BaseVariation("A")) == 2
    assert pileup.get_variation_count(1, BaseVariation("C")) == 1
    assert pileup.get_variation_count(1, BaseVariation("G")) == 1
    assert pileup.get_variation_count(2, BaseVariation("G")) == 2
    assert pileup.get_variation_count(3, BaseVariation("T")) == 2


def test_makepileup_detects_deletions_and_insertions():
    reference = Seq("ACGT")
    reads = [
        SeqRecord(Seq("ACGT")),
        SeqRecord(Seq("AGT")),  # deletion of C
        SeqRecord(Seq("ACGTT")),  # insertion of T after T
        SeqRecord(Seq("ATTGT")),  # insertion and substitution at the same spot
    ]
    pileup = MakePileup(reference, reads)

    assert pileup.option_count() == 4
    assert pileup.get_variation_count(1, DeletionVariation()) == 1
    assert pileup.get_variation_count(3, InsertionVariation("T")) == 1
    assert pileup.get_variation_count(1, InsertionVariation("T")) == 1
    assert pileup.get_variation_count(1, BaseVariation("T")) == 1
    assert pileup.get_variation_count(0, BaseVariation("A")) == 4
    assert pileup.get_variation_count(1, BaseVariation("C")) == 2
    assert pileup.get_variation_count(2, BaseVariation("G")) == 4
    assert pileup.get_variation_count(3, BaseVariation("T")) == 4


def test_pileup_clustering():
    reference = Seq("GAAAAG")
    reads = [
        SeqRecord(Seq("GCATAG")),   # cluster 1 GCAAAG: A->C at 1
        SeqRecord(Seq("GCAACGAG")),
        SeqRecord(Seq("GCAAGG")),
        SeqRecord(Seq("GCAATG")),    
        SeqRecord(Seq("GAAGATG")),   # cluster 2 GAAGAAG: G insertion at 3
        SeqRecord(Seq("TAAGAAG")),
        SeqRecord(Seq("GAAGAATTG")),
        SeqRecord(Seq("AGAAG")),

    ]
    pileup = MakePileup(reference, reads)
    print("DEBUG PILEUP:\n" + pileup.to_debug_string())
    denoised = DenoisePileup(pileup, threshold=0.3)
    print("DEBUG DENOISED PILEUP:\n" + denoised.to_debug_string())
    assert denoised.get_variation_count(1, BaseVariation("C")) == 4
    assert denoised.get_variation_count(1, BaseVariation("A")) == 4
    assert denoised.get_variation_count(3, InsertionVariation("G")) == 3
    assert denoised.get_variation_count(3, BaseVariation("A")) == 7


def test_read_pileup_intersect_returns_allowed_variants():
    reference = Seq("ACGT")
    read = SeqRecord(Seq("AGGTT"))
    pileup = MakePileup(reference, [SeqRecord(reference), read])
    assert pileup.option_count() == 2
    consensus = ReadPileupIntersect(pileup, read)
    assert consensus == Seq("AGGTT")


def test_read_pileup_intersect_filters_disallowed_variants():
    reference = Seq("ACGT")
    pileup = MakePileup(reference, [SeqRecord(reference)])
    assert pileup.option_count() == 0
    assert pileup.get_variation_count(1, BaseVariation("C")) == 1
    read = SeqRecord(Seq("AGGT"))
    consensus = ReadPileupIntersect(pileup, read)
    assert consensus == Seq("ACGT")


def test_read_pileup_intersect_requires_full_coverage():
    pileup = Pileup(Seq("ACGT"))
    assert pileup.option_count() == 0
    read = SeqRecord(Seq("ACG"))
    assert ReadPileupIntersect(pileup, read) is None


def test_makepileup_rejects_unknown_bases():
    reference = Seq("A")
    bad_read = SeqRecord(Seq("AZ"))
    with pytest.raises(ValueError):
        MakePileup(reference, [bad_read])
