import os
import struct
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import bamnostic
from bamnostic import bgzf

from ont_denoise import denoise_nanopore_reads


BASE_CODES = {
    "A": 1,
    "C": 2,
    "M": 3,
    "G": 4,
    "R": 5,
    "S": 6,
    "V": 7,
    "T": 8,
    "W": 9,
    "Y": 10,
    "H": 11,
    "K": 12,
    "D": 13,
    "B": 14,
    "N": 15,
}


def _reg2bin(beg: int, end: int) -> int:
    if end <= beg:
        end = beg + 1
    end -= 1
    if beg >> 14 == end >> 14:
        return 4681 + (beg >> 14)
    if beg >> 17 == end >> 17:
        return 585 + (beg >> 17)
    if beg >> 20 == end >> 20:
        return 73 + (beg >> 20)
    if beg >> 23 == end >> 23:
        return 9 + (beg >> 23)
    if beg >> 26 == end >> 26:
        return 1 + (beg >> 26)
    return 0


def _pack_cigar(cigar):
    return b"".join(struct.pack("<I", (length << 4) | op) for op, length in cigar)


def _encode_seq(seq: str) -> bytes:
    seq = seq.upper()
    buf = bytearray((len(seq) + 1) // 2)
    for i, base in enumerate(seq):
        code = BASE_CODES.get(base, 15)
        if i % 2 == 0:
            buf[i // 2] = code << 4
        else:
            buf[i // 2] |= code
    return bytes(buf)


def _encode_qual(quals):
    return bytes(quals)


def _write_single_read_bam(path: Path, ref_name: str, ref_seq: str, read_seq: str):
    reference_lengths = [len(ref_seq)]
    header_text = f"@HD\tVN:1.6\n@SQ\tSN:{ref_name}\tLN:{len(ref_seq)}\n"

    header_bytes = bytearray()
    header_bytes += b"BAM\x01"
    header_bytes += struct.pack("<i", len(header_text))
    header_bytes += header_text.encode()
    header_bytes += struct.pack("<i", 1)
    name_bytes = ref_name.encode() + b"\x00"
    header_bytes += struct.pack("<i", len(name_bytes))
    header_bytes += name_bytes
    header_bytes += struct.pack("<i", len(ref_seq))

    cigar = [(0, len(read_seq))]
    ref_start = 2
    bin_mq_nl = (_reg2bin(ref_start, ref_start + len(read_seq)) << 16) | (60 << 8) | (
        len("read1") + 1
    )
    flag_nc = (0 << 16) | len(cigar)
    quals = [30] * len(read_seq)

    record = bytearray()
    record += struct.pack("<i", 0)  # ref id
    record += struct.pack("<i", ref_start)
    record += struct.pack("<I", bin_mq_nl)
    record += struct.pack("<I", flag_nc)
    record += struct.pack("<I", len(read_seq))
    record += struct.pack("<i", -1)
    record += struct.pack("<i", -1)
    record += struct.pack("<i", 0)
    record += b"read1\x00"
    record += _pack_cigar(cigar)
    record += _encode_seq(read_seq)
    record += _encode_qual(quals)
    packed = struct.pack("<i", len(record)) + record

    with bgzf.BgzfWriter(path, "wb") as writer:
        writer.write(bytes(header_bytes))
        writer.write(packed)


class SingleReadDenoiseTest(unittest.TestCase):
    def test_single_read_matches_reference(self):
        read_seq = "ACGTA"
        ref_seq = "TT" + read_seq + "TT"

        with TemporaryDirectory() as tmpdir:
            base = Path(tmpdir)
            ref_path = base / "ref.fa"
            bam_path = base / "simple.bam"
            out_fastq = base / "out.fastq"

            ref_path.write_text(f">contig1\n{ref_seq}\n", encoding="utf-8")
            _write_single_read_bam(bam_path, "contig1", ref_seq, read_seq)

            denoise_nanopore_reads(
                bam_path=str(bam_path),
                reference_fasta=str(ref_path),
                output_fastq=str(out_fastq),
                min_variant_frequency=0.5,
                min_read_depth=1,
            )

            self.assertTrue(out_fastq.exists())
            with out_fastq.open("r", encoding="utf-8") as handle:
                lines = [line.strip() for line in handle]

            self.assertEqual(len(lines), 4)
            self.assertEqual(lines[0], "@read1_OD")
            self.assertEqual(lines[1], read_seq)
            self.assertEqual(lines[2], "+")
            self.assertEqual(lines[3], "".join(chr(30 + 33) for _ in read_seq))


if __name__ == "__main__":
    unittest.main()
