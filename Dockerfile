FROM mambaorg/micromamba:1.5.6

LABEL description="Long-read QC pipeline for Oxford Nanopore data"

USER root
RUN apt-get update && apt-get install -y --no-install-recommends procps \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /pipeline

COPY config/conda_env.yml /pipeline/conda_env.yml
COPY Snakefile config/ scripts/ /pipeline/

RUN micromamba install --yes --name base --file /pipeline/conda_env.yml \
    && micromamba clean --all --yes

ENV PATH="/opt/conda/bin:$PATH"
RUN mkdir -p /pipeline/data /pipeline/results

ENTRYPOINT ["snakemake", "--cores", "4"]
CMD ["--config", "fastq_input=/pipeline/data/barcode77.fastq"]
