import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from ontcluster.clustering import cluster_fastq


def test_cluster_fastq_creates_expected_clusters(tmp_path):
    reference = Seq("ACGT" * 15)
    def make_record(seq: str, identifier: str) -> SeqRecord:
        record = SeqRecord(Seq(seq), id=identifier, description="")
        record.letter_annotations["phred_quality"] = [40] * len(seq)
        return record

    records = [
        make_record("ACGT" * 15, "read1"),
        make_record("ACGT" * 15, "read2"),
        make_record("AGGT" * 15, "read3"),
        make_record("AGGT" * 15, "read4"),
    ]
    fastq_path = tmp_path / "reads.fastq"
    SeqIO.write(records, fastq_path, "fastq")

    consensus_path = tmp_path / "consensus.fa"
    result = cluster_fastq(
        reference,
        fastq_path,
        threshold=0.1,
        max_clusters=10,
        output_dir=tmp_path,
        consensus_fasta=consensus_path,
    )

    assert result.total_reads == 4
    assert result.unaligned_reads == 0
    assert len(result.clusters) == 2
    assert result.consensus_fasta == consensus_path
    assert consensus_path.exists()
    consensus_records = list(SeqIO.parse(consensus_path, "fasta"))
    assert len(consensus_records) == 2

    for cluster_id, seq in result.clusters.items():
        output_file = result.output_files[cluster_id]
        assert output_file.exists()
        written_records = list(SeqIO.parse(output_file, "fastq"))
        assert len(written_records) == result.cluster_counts[cluster_id]
        written_sequences = {str(rec.seq) for rec in written_records}
        if str(seq) == "ACGT" * 15:
            assert written_sequences == {"ACGT" * 15}
        elif str(seq) == "AGGT" * 15:
            assert written_sequences == {"AGGT" * 15}
        else:
            raise AssertionError("Unexpected cluster sequence")


def test_cluster_fastq_max_clusters(tmp_path):
    reference = Seq("ACGT" * 20)

    def make_record(seq: str, identifier: str) -> SeqRecord:
        record = SeqRecord(Seq(seq), id=identifier, description="")
        record.letter_annotations["phred_quality"] = [40] * len(seq)
        return record

    records = [
        make_record(str(reference), "read1"),
        make_record("TGCA" * 20, "read2"),
        make_record("CAGT" * 20, "read3"),
    ]
    fastq_path = tmp_path / "max.fastq"
    SeqIO.write(records, fastq_path, "fastq")

    with pytest.raises(ValueError) as excinfo:
        cluster_fastq(
            reference,
            fastq_path,
            threshold=0.05,
            max_clusters=1,
            output_dir=tmp_path,
        )
    assert "Maximum cluster count" in str(excinfo.value)
