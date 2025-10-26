"""
Quick smoke test for the ont-denoise library.

This script runs the denoiser against the Influenza A dataset provided in
../Q8/infa/ using conservative filter thresholds and reports the number of
reads emitted.
"""

from pathlib import Path

def main() -> None:
    from ont_denoise import denoise_nanopore_reads

    base_dir = Path("../Q8/infa")
    bam_path = base_dir / "unc-INFA-H1N1.bam"
    reference = base_dir / "ref.fa"
    output_fastq = Path("unc-INFA-H1N1-denoised.fastq")

    denoise_nanopore_reads(
        bam_path=str(bam_path),
        reference_fasta=str(reference),
        output_fastq=str(output_fastq),
        min_variant_frequency=0.1,
        min_read_depth=5,
    )

    if output_fastq.exists():
        with output_fastq.open("r", encoding="utf-8") as handle:
            line_count = sum(1 for _ in handle)
        read_count = line_count // 4
        print(f"Wrote {read_count} denoised reads to {output_fastq}")
    else:
        raise FileNotFoundError(f"Expected output FASTQ not found: {output_fastq}")


if __name__ == "__main__":
    main()
