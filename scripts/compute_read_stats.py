# Dependencies: biopython, pandas, numpy (all provided by conda_env.yml)
# =============================================================================

import argparse
import gzip
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO

# ---------------------------------------------------------------------------
# Logging configuration
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def compute_gc_content(sequence: str) -> float:
    """
    Calculate GC content percentage for a DNA/RNA sequence.

    Parameters
    ----------
    sequence : str
        Nucleotide sequence (upper or lower case).

    Returns
    -------
    float
        GC content as a percentage (0–100). Returns 0.0 for empty sequences.
    """
    seq_upper = sequence.upper()
    length = len(seq_upper)
    if length == 0:
        return 0.0
    gc_count = seq_upper.count("G") + seq_upper.count("C")
    return round((gc_count / length) * 100, 4)


def compute_mean_quality(phred_scores: list) -> float:
    """
    Compute the arithmetic mean of a list of Phred quality scores.

    Bio.SeqIO already converts the ASCII-encoded quality string into integer
    Phred scores (Q), so no additional conversion is required here.

    Parameters
    ----------
    phred_scores : list of int
        Phred quality scores as integers.

    Returns
    -------
    float
        Mean Phred quality score, rounded to 4 decimal places.
        Returns 0.0 for reads with no quality scores.
    """
    if not phred_scores:
        return 0.0
    return round(float(np.mean(phred_scores)), 4)


def open_fastq(path: Path):
    """
    Open a FASTQ file transparently, supporting both plain-text and gzip.

    Parameters
    ----------
    path : Path
        Path to the FASTQ file (.fastq or .fastq.gz).

    Returns
    -------
    file-like object
        An open file handle in text mode.
    """
    suffix = path.suffix.lower()
    if suffix == ".gz":
        logger.info("Detected gzip-compressed FASTQ.")
        return gzip.open(path, "rt")
    return open(path, "r")


# ---------------------------------------------------------------------------
# Core pipeline function
# ---------------------------------------------------------------------------

def parse_fastq(fastq_path: Path) -> pd.DataFrame:
    """
    Iterate over all reads in a FASTQ file and compute per-read statistics.

    Parameters
    ----------
    fastq_path : Path
        Path to the input FASTQ file.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns: read_id, gc_content_pct, read_length,
        mean_quality_score.
    """
    records = []
    read_count = 0
    skipped_count = 0

    logger.info(f"Parsing FASTQ file: {fastq_path}")

    with open_fastq(fastq_path) as handle:
        for record in SeqIO.parse(handle, "fastq"):
            read_count += 1

            sequence = str(record.seq)
            read_length = len(sequence)

            # Skip degenerate reads (empty or no quality information)
            if read_length == 0:
                logger.warning(f"Skipping empty read: {record.id}")
                skipped_count += 1
                continue

            phred_scores = record.letter_annotations.get("phred_quality", [])
            if not phred_scores:
                logger.warning(f"Skipping read with no quality scores: {record.id}")
                skipped_count += 1
                continue

            gc_pct = compute_gc_content(sequence)
            mean_qual = compute_mean_quality(phred_scores)

            records.append(
                {
                    "read_id": record.id,
                    "gc_content_pct": gc_pct,
                    "read_length": read_length,
                    "mean_quality_score": mean_qual,
                }
            )

            # Progress log every 10,000 reads to avoid silent hangs on large files
            if read_count % 10_000 == 0:
                logger.info(f"  Processed {read_count:,} reads so far …")

    logger.info(
        f"Finished parsing. Total reads: {read_count:,} | "
        f"Processed: {len(records):,} | Skipped: {skipped_count:,}"
    )

    df = pd.DataFrame(
        records,
        columns=["read_id", "gc_content_pct", "read_length", "mean_quality_score"],
    )
    return df


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Compute per-read statistics (GC content, length, mean quality) "
            "from a FASTQ file and write results to a CSV."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--input", "-i",
        required=True,
        type=Path,
        help="Path to input FASTQ file (.fastq or .fastq.gz).",
    )
    parser.add_argument(
        "--output", "-o",
        required=True,
        type=Path,
        help="Path for the output CSV file.",
    )
    return parser.parse_args()


def main() -> None:
    """Main entry point for the read statistics script."""
    args = parse_args()

    # Validate input file exists
    if not args.input.exists():
        logger.error(f"Input file not found: {args.input}")
        sys.exit(1)

    # Ensure output directory exists
    args.output.parent.mkdir(parents=True, exist_ok=True)

    # Parse FASTQ and compute statistics
    df = parse_fastq(args.input)

    if df.empty:
        logger.error("No valid reads found in the FASTQ file. Exiting.")
        sys.exit(1)

    # Write CSV output
    df.to_csv(args.output, index=False)
    logger.info(f"Read statistics written to: {args.output}")
    logger.info(f"  Rows in CSV: {len(df):,}")

    # Print a brief summary to stdout for the Snakemake log
    logger.info("── Quick summary ──────────────────────────────────────────")
    logger.info(f"  Read count     : {len(df):,}")
    logger.info(f"  Mean length    : {df['read_length'].mean():,.1f} bp")
    logger.info(f"  Median length  : {df['read_length'].median():,.1f} bp")
    logger.info(f"  Mean GC %%      : {df['gc_content_pct'].mean():.2f} %%")
    logger.info(f"  Mean quality   : {df['mean_quality_score'].mean():.2f}")
    logger.info("───────────────────────────────────────────────────────────")


if __name__ == "__main__":
    main()
