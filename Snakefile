# =============================================================================
# Snakefile: Mini-Bioinformatics Pipeline for Long-Read QC
# Author: Bioinformatics Pipeline
# Description: Quality Control pipeline for Oxford Nanopore (long-read) data.
#              Runs NanoPlot for standard QC and a custom Python script for
#              per-read statistics (GC content, length, mean quality).
# Usage: snakemake --cores 4 --use-conda
# =============================================================================

import os

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
configfile: "config/config.yaml"

SAMPLE   = config["sample_name"]
FASTQ    = config["fastq_input"]
OUTDIR   = config["output_dir"]

# ---------------------------------------------------------------------------
# Target rule – defines the final outputs the pipeline must produce
# ---------------------------------------------------------------------------
rule all:
    input:
        # NanoPlot QC report
        expand("{outdir}/nanoplot/{sample}/NanoPlot-report.html", outdir=OUTDIR, sample=SAMPLE),
        # Custom per-read statistics CSV
        expand("{outdir}/read_stats/{sample}_read_stats.csv",     outdir=OUTDIR, sample=SAMPLE),
        # Visualisation plots
        expand("{outdir}/plots/{sample}_gc_content.png",          outdir=OUTDIR, sample=SAMPLE),
        expand("{outdir}/plots/{sample}_read_length.png",         outdir=OUTDIR, sample=SAMPLE),
        expand("{outdir}/plots/{sample}_mean_quality.png",        outdir=OUTDIR, sample=SAMPLE),
        # Summary statistics text file
        expand("{outdir}/plots/{sample}_summary_stats.txt",       outdir=OUTDIR, sample=SAMPLE),


# ---------------------------------------------------------------------------
# Rule 1: NanoPlot – long-read specific QC
# ---------------------------------------------------------------------------
rule nanoplot_qc:
    """
    Run NanoPlot on the raw FASTQ file to produce an HTML QC report.
    NanoPlot is purpose-built for Oxford Nanopore long-read data.
    """
    input:
        fastq = FASTQ
    output:
        report = f"{OUTDIR}/nanoplot/{SAMPLE}/NanoPlot-report.html"
    params:
        outdir = f"{OUTDIR}/nanoplot/{SAMPLE}",
        prefix = SAMPLE
    conda:
        "config/conda_env.yml"
    threads: 4
    log:
        f"{OUTDIR}/logs/nanoplot_{SAMPLE}.log"
    shell:
        """
        mkdir -p {params.outdir}
        NanoPlot \
            --fastq {input.fastq} \
            --outdir {params.outdir} \
            --prefix {params.prefix} \
            --threads {threads} \
            --plots dot \
            --N50 \
            --title "{SAMPLE} NanoPlot QC" \
            > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# Rule 2: Custom per-read statistics
# ---------------------------------------------------------------------------
rule compute_read_stats:
    """
    Run the custom Python script to calculate per-read GC content (%),
    read length, and mean Phred quality score. Output is a CSV file.
    """
    input:
        fastq = FASTQ
    output:
        csv = f"{OUTDIR}/read_stats/{SAMPLE}_read_stats.csv"
    params:
        outdir = f"{OUTDIR}/read_stats"
    conda:
        "config/conda_env.yml"
    log:
        f"{OUTDIR}/logs/read_stats_{SAMPLE}.log"
    shell:
        """
        mkdir -p {params.outdir}
        python scripts/compute_read_stats.py \
            --input  {input.fastq} \
            --output {output.csv} \
            > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# Rule 3: Visualisation
# ---------------------------------------------------------------------------
rule visualize_stats:
    """
    Read the per-read CSV and generate distribution plots for GC content,
    read length (log-scaled), and mean quality score. Also prints summary
    statistics (mean, median, N50) to a text file.
    """
    input:
        csv = f"{OUTDIR}/read_stats/{SAMPLE}_read_stats.csv"
    output:
        gc_plot      = f"{OUTDIR}/plots/{SAMPLE}_gc_content.png",
        len_plot     = f"{OUTDIR}/plots/{SAMPLE}_read_length.png",
        qual_plot    = f"{OUTDIR}/plots/{SAMPLE}_mean_quality.png",
        summary_txt  = f"{OUTDIR}/plots/{SAMPLE}_summary_stats.txt"
    params:
        outdir = f"{OUTDIR}/plots",
        prefix = SAMPLE
    conda:
        "config/conda_env.yml"
    log:
        f"{OUTDIR}/logs/visualize_{SAMPLE}.log"
    shell:
        """
        mkdir -p {params.outdir}
        python scripts/visualize_stats.py \
            --input  {input.csv} \
            --outdir {params.outdir} \
            --prefix {params.prefix} \
            > {log} 2>&1
        """
