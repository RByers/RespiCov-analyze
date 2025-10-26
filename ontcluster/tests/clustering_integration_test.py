from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from ontcluster.clustering import cluster_fastq


def test_cluster_fastq_creates_expected_clusters(tmp_path):
    reference = Seq("ACGT")
    def make_record(seq: str, identifier: str) -> SeqRecord:
        record = SeqRecord(Seq(seq), id=identifier, description="")
        record.letter_annotations["phred_quality"] = [40] * len(seq)
        return record

    records = [
        make_record("ACGT", "read1"),
        make_record("ACGT", "read2"),
        make_record("AGGT", "read3"),
        make_record("AGGT", "read4"),
    ]
    fastq_path = tmp_path / "reads.fastq"
    SeqIO.write(records, fastq_path, "fastq")

    result = cluster_fastq(reference, fastq_path, threshold=0.1, output_dir=tmp_path)

    assert result.total_reads == 4
    assert result.unaligned_reads == 0
    assert len(result.clusters) == 2

    for cluster_id, seq in result.clusters.items():
        output_file = result.output_files[cluster_id]
        assert output_file.exists()
        written_records = list(SeqIO.parse(output_file, "fastq"))
        assert len(written_records) == result.cluster_counts[cluster_id]
        written_sequences = {str(rec.seq) for rec in written_records}
        if str(seq) == "ACGT":
            assert written_sequences == {"ACGT"}
        elif str(seq) == "AGGT":
            assert written_sequences == {"AGGT"}
        else:
            raise AssertionError("Unexpected cluster sequence")
