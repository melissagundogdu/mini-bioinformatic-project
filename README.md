# 🧬 Mini-Bioinformatics Long-Read QC Pipeline

A reproducible, production-grade quality control pipeline for Oxford Nanopore
long-read sequencing data. Built with **Snakemake**, **NanoPlot**, and custom
Python analysis scripts.

---

## 📁 Repository Structure

```
mini-bioinformatics-pipeline/
├── Snakefile                    # Snakemake workflow definition
├── Dockerfile                   # Docker image for full reproducibility
├── email_to_professor.md        # Non-technical summary for Professor Kılıç
│
├── config/
│   ├── config.yaml              # Pipeline configuration (input/output paths)
│   └── conda_env.yml            # Conda environment with all dependencies
│
├── scripts/
│   ├── compute_read_stats.py    # Custom per-read statistics (GC, length, quality)
│   └── visualize_stats.py       # Distribution plots + summary statistics
│
├── data/                        # ← Place your FASTQ file here
│   └── barcode77.fastq
│
└── results/                     # Auto-created by the pipeline
    ├── nanoplot/barcode77/      # NanoPlot HTML report + plots
    ├── read_stats/              # Per-read CSV file
    ├── plots/                   # Distribution PNGs + summary TXT
    └── logs/                    # Snakemake rule logs
```

---

## ⚡ Quick Start (Single Command)

### Option A — Conda (Recommended)

**1. Install Miniconda** (skip if already installed):
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

**2. Clone the repository and place your data:**
```bash
git clone https://github.com/melissagundogdu/mini-bioinformatic-project
cd mini-bioinformatics-pipeline
mkdir -p data
cp /path/to/barcode77.fastq data/
```

**3. Run the entire pipeline with one command:**
```bash
snakemake --cores 4 --use-conda --conda-frontend conda
```

That's it. Snakemake will automatically create the `longread-qc` Conda
environment on the first run (takes ~5 minutes) and then execute all three
pipeline steps.

---

### Option B — Docker

**1. Build the image:**
```bash
docker build -t longread-qc:1.0 .
```

**2. Run the pipeline:**
```bash
docker run --rm \
  -v "$(pwd)/data:/pipeline/data" \
  -v "$(pwd)/results:/pipeline/results" \
  longread-qc:1.0
```

---

## 🔧 Configuration

Edit `config/config.yaml` to change the sample name or input path:

```yaml
sample_name: "barcode77"            # Used as prefix for all output files
fastq_input: "data/barcode77.fastq" # Supports .fastq and .fastq.gz
output_dir:  "results"
```

---

## 📊 Pipeline Steps

### Step 1 — NanoPlot QC (`nanoplot_qc`)

Runs [NanoPlot](https://github.com/wdecoster/NanoPlot) — the gold-standard QC
tool for Oxford Nanopore data — to produce an interactive HTML report.

**Output:** `results/nanoplot/barcode77/barcode77-NanoPlot-report.html`

### Step 2 — Per-Read Statistics (`compute_read_stats`)

Custom Python script using `Bio.SeqIO` to calculate for **every individual read**:

| Column               | Description                                |
|----------------------|--------------------------------------------|
| `read_id`            | Sequence identifier from the FASTQ header  |
| `gc_content_pct`     | GC bases / total bases × 100              |
| `read_length`        | Number of bases                            |
| `mean_quality_score` | Arithmetic mean of Phred Q scores          |

**Output:** `results/read_stats/barcode77_read_stats.csv`

### Step 3 — Visualisation (`visualize_stats`)

Generates three publication-quality distribution plots and a summary statistics
text file:

| Plot | Description |
|------|-------------|
| `barcode77_gc_content.png`  | Histogram + KDE; median and typical GC band annotated |
| `barcode77_read_length.png` | Log-scaled histogram + KDE; median and N50 annotated  |
| `barcode77_mean_quality.png`| Histogram + KDE; Q10 and Q20 thresholds annotated     |

**Printed statistics:** read count, mean, median, std, min, max, and **N50**
for all three metrics.

---

## 📦 Dependencies

All dependencies are managed automatically via the Conda environment
(`config/conda_env.yml`):

| Package      | Version | Purpose                        |
|--------------|---------|--------------------------------|
| snakemake    | 7.32.4  | Workflow orchestration         |
| nanoplot     | 1.42.0  | Long-read QC                   |
| python       | 3.11    | Runtime                        |
| biopython    | 1.83    | FASTQ parsing via SeqIO        |
| pandas       | 2.2.1   | CSV I/O and data manipulation  |
| numpy        | 1.26.4  | Numerical operations           |
| matplotlib   | 3.8.3   | Plotting backend               |
| seaborn      | 0.13.2  | Statistical visualisation      |
| scipy        | 1.13.0  | KDE for read-length plot       |

---

---

## 📈 Results — barcode77 Sample

Pipeline was executed on `barcode77.fastq` (Oxford Nanopore, R10.4.1 chemistry).

| Metric                | Value     | Assessment          |
|-----------------------|-----------|---------------------|
| Total reads           | 81,011    | ✅ Good coverage    |
| Median read length    | 547 bp    | ⚠️ Short but normal |
| N50 read length       | 1,761 bp  | ✅ Sufficient       |
| Median GC content     | 53.5%     | ✅ Expected range   |
| Median quality score  | Q17.31    | ✅ Above Q10 threshold |

**Conclusion:** Data quality is sufficient to proceed to alignment with Minimap2 (`map-ont` preset).

---

### Recommendation Logic

- **Proceed to alignment** if median Q ≥ Q10 AND N50 ≥ 1,000 bp.
- **High-confidence variant calling** requires median Q ≥ Q20.
- **Investigate contamination** if the GC distribution is bimodal or very wide.
- **Consider re-sequencing** if > 30 % of reads are below Q7 or shorter than 200 bp.

---

## 📧 Communication

See [`email_to_professor.md`](email_to_professor.md) for a plain-language
summary of results and recommendations prepared for Professor Kılıç.

---

## 🛠 Running Individual Steps

```bash
# Only run NanoPlot QC
snakemake --cores 4 --use-conda results/nanoplot/barcode77/NanoPlot-report.html

# Only compute per-read statistics
snakemake --cores 4 --use-conda results/read_stats/barcode77_read_stats.csv

# Only generate plots (requires CSV from step above)
snakemake --cores 4 --use-conda results/plots/barcode77_gc_content.png

# Dry-run (see what would be executed without running)
snakemake --cores 4 --use-conda --dry-run

# Generate a visual DAG of the pipeline
snakemake --dag | dot -Tpng > pipeline_dag.png
```

---


