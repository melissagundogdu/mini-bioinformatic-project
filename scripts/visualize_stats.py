# Dependencies: pandas, numpy, matplotlib, seaborn (all in conda_env.yml)
# =============================================================================

import argparse
import logging
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")          # Non-interactive backend for server/pipeline use
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
import seaborn as sns

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
# Plot aesthetics (consistent across all figures)
# ---------------------------------------------------------------------------
PALETTE_PRIMARY   = "#2E86AB"   # steel blue  – bars / histograms
PALETTE_KDE       = "#E84855"   # vivid red   – KDE line
PALETTE_MEDIAN    = "#F4A261"   # warm orange – median reference line
PALETTE_N50       = "#6A4C93"   # violet      – N50 reference line (length plot)
FIGURE_DPI        = 150
FIGURE_SIZE       = (10, 6)
FONT_FAMILY       = "DejaVu Sans"

sns.set_theme(style="whitegrid", font=FONT_FAMILY)
plt.rcParams.update(
    {
        "axes.titlesize": 15,
        "axes.titleweight": "bold",
        "axes.labelsize": 12,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "figure.dpi": FIGURE_DPI,
    }
)


# ---------------------------------------------------------------------------
# Statistical helper functions
# ---------------------------------------------------------------------------

def calculate_n50(lengths: pd.Series) -> float:
    """
    Calculate the N50 of a read-length distribution.

    N50 is the length L such that reads of length >= L account for at least
    50 % of the total sequenced bases. It is a standard metric in long-read
    sequencing to summarise assembly or read contiguity.

    Parameters
    ----------
    lengths : pd.Series
        Series of integer read lengths.

    Returns
    -------
    float
        N50 value in base pairs.
    """
    sorted_lengths = np.sort(lengths.values)[::-1]   # descending
    cumulative     = np.cumsum(sorted_lengths)
    half_total     = cumulative[-1] / 2.0
    n50_idx        = np.searchsorted(cumulative, half_total)
    return float(sorted_lengths[n50_idx])


def summary_statistics(df: pd.DataFrame) -> dict:
    """
    Compute descriptive statistics for all three metrics and the read N50.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with columns: gc_content_pct, read_length, mean_quality_score.

    Returns
    -------
    dict
        Nested dict keyed by metric name, with sub-keys: mean, median, std,
        min, max.  Read length also contains 'n50'.
    """
    stats = {}
    metrics = {
        "GC Content (%)":        "gc_content_pct",
        "Read Length (bp)":      "read_length",
        "Mean Quality Score":    "mean_quality_score",
    }
    for label, col in metrics.items():
        s = df[col]
        entry = {
            "mean":   round(s.mean(),   2),
            "median": round(s.median(), 2),
            "std":    round(s.std(),    2),
            "min":    round(s.min(),    2),
            "max":    round(s.max(),    2),
        }
        if col == "read_length":
            entry["n50"] = round(calculate_n50(s), 0)
        stats[label] = entry
    return stats


def format_summary_report(stats: dict, total_reads: int) -> str:
    """
    Format the statistics dictionary as a human-readable text report.

    Parameters
    ----------
    stats : dict
        Output of summary_statistics().
    total_reads : int
        Total number of reads in the dataset.

    Returns
    -------
    str
        Formatted multi-line string.
    """
    lines = [
        "=" * 60,
        "  LONG-READ QC SUMMARY STATISTICS",
        "=" * 60,
        f"  Total reads analysed : {total_reads:,}",
        "",
    ]
    for metric, values in stats.items():
        lines.append(f"  {metric}")
        lines.append(f"    Mean   : {values['mean']:>12,.2f}")
        lines.append(f"    Median : {values['median']:>12,.2f}")
        lines.append(f"    Std    : {values['std']:>12,.2f}")
        lines.append(f"    Min    : {values['min']:>12,.2f}")
        lines.append(f"    Max    : {values['max']:>12,.2f}")
        if "n50" in values:
            lines.append(f"    N50    : {values['n50']:>12,.0f}")
        lines.append("")
    lines.append("=" * 60)
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Plotting functions
# ---------------------------------------------------------------------------

def _add_vline_with_label(ax, value: float, color: str, label: str) -> None:
    """Draw a vertical reference line and annotate it."""
    ax.axvline(value, color=color, linestyle="--", linewidth=1.6, label=label)


def plot_gc_content(df: pd.DataFrame, outpath: Path) -> None:
    """
    Histogram + KDE for GC content distribution.

    A vertical dashed line marks the median and a shaded band indicates the
    typical expected GC range for bacterial genomes (40–65 %) as a visual
    reference.

    Parameters
    ----------
    df       : DataFrame containing 'gc_content_pct' column.
    outpath  : Path where the PNG will be saved.
    """
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)

    sns.histplot(
        data=df,
        x="gc_content_pct",
        bins=50,
        kde=True,
        color=PALETTE_PRIMARY,
        edgecolor="white",
        linewidth=0.4,
        kde_kws={},
        ax=ax,
    )

    median_gc = df["gc_content_pct"].median()
    _add_vline_with_label(ax, median_gc, PALETTE_MEDIAN, f"Median: {median_gc:.1f} %")

    # Typical bacterial GC reference band
    ax.axvspan(40, 65, alpha=0.07, color="green", label="Typical GC range (40–65 %)")

    ax.set_title("GC Content Distribution")
    ax.set_xlabel("GC Content (%)")
    ax.set_ylabel("Number of Reads")
    ax.legend(framealpha=0.85)
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter("%.0f%%"))

    fig.tight_layout()
    fig.savefig(outpath, dpi=FIGURE_DPI)
    plt.close(fig)
    logger.info(f"Saved GC content plot → {outpath}")


def plot_read_length(df: pd.DataFrame, outpath: Path, n50: float) -> None:
    """
    Log-scaled histogram + KDE for read length distribution.

    Log scaling is used because ONT read lengths often span several orders of
    magnitude. Vertical lines mark the median and N50.

    Parameters
    ----------
    df      : DataFrame containing 'read_length' column.
    outpath : Path where the PNG will be saved.
    n50     : Pre-calculated N50 value to annotate.
    """
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)

    # Compute log10 of lengths for binning; back-transform axes labels
    log_lengths = np.log10(df["read_length"].clip(lower=1))

    # Build bins in log space for even visual distribution
    bins = np.linspace(log_lengths.min(), log_lengths.max(), 60)

    ax.hist(
        log_lengths,
        bins=bins,
        color=PALETTE_PRIMARY,
        edgecolor="white",
        linewidth=0.4,
        label="Read count",
    )

    # KDE over log-transformed data
    from scipy.stats import gaussian_kde  # local import; scipy in conda_env.yml
    kde = gaussian_kde(log_lengths, bw_method="scott")
    x_range = np.linspace(log_lengths.min(), log_lengths.max(), 400)
    # Scale KDE to histogram area
    bin_width = bins[1] - bins[0]
    scale = len(log_lengths) * bin_width
    ax.plot(x_range, kde(x_range) * scale, color=PALETTE_KDE, linewidth=2, label="KDE")

    # Reference lines
    median_len = df["read_length"].median()
    _add_vline_with_label(ax, np.log10(median_len),  PALETTE_MEDIAN, f"Median: {median_len:,.0f} bp")
    _add_vline_with_label(ax, np.log10(n50),         PALETTE_N50,   f"N50: {n50:,.0f} bp")

    # Replace log-scale tick labels with human-readable bp values
    tick_vals = [2, 3, 4, 5, 6]   # 100, 1k, 10k, 100k, 1M
    ax.set_xticks(tick_vals)
    ax.set_xticklabels([f"{10**v:,.0f}" for v in tick_vals])

    ax.set_title("Read Length Distribution (log scale)")
    ax.set_xlabel("Read Length (bp, log scale)")
    ax.set_ylabel("Number of Reads")
    ax.legend(framealpha=0.85)

    fig.tight_layout()
    fig.savefig(outpath, dpi=FIGURE_DPI)
    plt.close(fig)
    logger.info(f"Saved read length plot → {outpath}")


def plot_mean_quality(df: pd.DataFrame, outpath: Path) -> None:
    """
    Histogram + KDE for mean Phred quality score distribution.

    A vertical reference line marks Q10 (90 % accuracy per base, a common
    ONT minimum) and Q20 (99 % accuracy, high-quality threshold).

    Parameters
    ----------
    df      : DataFrame containing 'mean_quality_score' column.
    outpath : Path where the PNG will be saved.
    """
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)

    sns.histplot(
        data=df,
        x="mean_quality_score",
        bins=50,
        kde=True,
        color=PALETTE_PRIMARY,
        edgecolor="white",
        linewidth=0.4,
        kde_kws={},
        ax=ax,
    )

    median_q = df["mean_quality_score"].median()
    _add_vline_with_label(ax, median_q, PALETTE_MEDIAN, f"Median: {median_q:.1f}")
    # Standard ONT quality thresholds
    ax.axvline(10, color="#555555", linestyle=":", linewidth=1.4, label="Q10 threshold (90 % accuracy)")
    ax.axvline(20, color="#2ca02c", linestyle=":", linewidth=1.4, label="Q20 threshold (99 % accuracy)")

    ax.set_title("Mean Read Quality Score Distribution")
    ax.set_xlabel("Mean Phred Quality Score (Q)")
    ax.set_ylabel("Number of Reads")
    ax.legend(framealpha=0.85)

    fig.tight_layout()
    fig.savefig(outpath, dpi=FIGURE_DPI)
    plt.close(fig)
    logger.info(f"Saved quality score plot → {outpath}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Generate distribution plots and summary statistics from a "
            "per-read stats CSV produced by compute_read_stats.py."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--input", "-i",
        required=True,
        type=Path,
        help="Path to input CSV (output of compute_read_stats.py).",
    )
    parser.add_argument(
        "--outdir", "-d",
        required=True,
        type=Path,
        help="Directory where output PNG files will be saved.",
    )
    parser.add_argument(
        "--prefix", "-p",
        default="sample",
        help="Prefix for output filenames (e.g., sample name).",
    )
    return parser.parse_args()


def main() -> None:
    """Main entry point for the visualisation script."""
    args = parse_args()

    if not args.input.exists():
        logger.error(f"Input CSV not found: {args.input}")
        sys.exit(1)

    args.outdir.mkdir(parents=True, exist_ok=True)

    # ── Load data ──────────────────────────────────────────────────────────
    logger.info(f"Loading statistics from: {args.input}")
    df = pd.read_csv(args.input)
    logger.info(f"Loaded {len(df):,} reads.")

    # ── Compute statistics ─────────────────────────────────────────────────
    stats        = summary_statistics(df)
    summary_text = format_summary_report(stats, total_reads=len(df))

    # Print to stdout (captured in Snakemake log)
    print(summary_text)

    # Save to text file
    txt_path = args.outdir / f"{args.prefix}_summary_stats.txt"
    txt_path.write_text(summary_text)
    logger.info(f"Summary statistics saved → {txt_path}")

    # ── Generate plots ─────────────────────────────────────────────────────
    n50 = stats["Read Length (bp)"]["n50"]

    plot_gc_content(
        df,
        outpath=args.outdir / f"{args.prefix}_gc_content.png",
    )
    plot_read_length(
        df,
        outpath=args.outdir / f"{args.prefix}_read_length.png",
        n50=n50,
    )
    plot_mean_quality(
        df,
        outpath=args.outdir / f"{args.prefix}_mean_quality.png",
    )

    logger.info("All visualisations complete.")


if __name__ == "__main__":
    main()

