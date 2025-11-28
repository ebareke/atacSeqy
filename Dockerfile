# ============================================================
# Dockerfile for atacSeqy — Autonomous ATAC-seq Processing Engine
# Base: micromamba + bioconda stack
# ============================================================

FROM mambaorg/micromamba:latest

LABEL maintainer="Dr. Eric Bareke <ebareke@github>"
LABEL description="Container image for atacSeqy ATAC-seq pipeline"

# Ensure micromamba auto-activates base env in Dockerfile RUNs
ENV MAMBA_DOCKERFILE_ACTIVATE=1

# Working directory inside container
WORKDIR /workspace

# ------------------------------------------------------------
# 1. Install core tools via micromamba (bioconda + conda-forge)
# ------------------------------------------------------------
RUN micromamba install -y -n base -c conda-forge -c bioconda \
    bash \
    bwa \
    samtools \
    bedtools \
    macs2 \
    deeptools \
    yq \
    fastp \
    fastqc \
    multiqc \
    python=3.11 \
    r-base \
    r-essentials \
    && micromamba clean -afy

# ------------------------------------------------------------
# 2. Install R/Bioconductor packages (ArchR, chromVAR, DESeq2)
#    You can comment this block out if you want a lighter image.
# ------------------------------------------------------------
RUN R -e 'install.packages("BiocManager", repos = "https://cloud.r-project.org"); \
          BiocManager::install(c("ArchR", "chromVAR", "SummarizedExperiment", "DESeq2"), ask = FALSE, update = FALSE)'

# ------------------------------------------------------------
# 3. Copy atacSeqy repository into container
# ------------------------------------------------------------
# Assumes Docker build context is the repo root.
COPY . /workspace

# Make sure scripts are executable
RUN chmod +x /workspace/run.sh || true && \
    chmod +x /workspace/validate_config.sh || true

# ------------------------------------------------------------
# 4. Environment variables (optional)
# ------------------------------------------------------------
# Set a default output directory (can be overridden at runtime)
ENV ATACSEQY_OUTDIR=/workspace/results

# Add base environment to PATH explicitly (micromamba already handles this,
# but this makes it clearer and robust for some runners).
ENV PATH=/opt/conda/bin:$PATH

# ------------------------------------------------------------
# 5. Default command
# ------------------------------------------------------------
# By default drop into a shell; user can run:
#   bash run.sh --config config.yaml --samples samples.csv
#
# Example with data mounting:
#   docker run --rm -it \
#     -v /path/to/data:/data \
#     ghcr.io/ebareke/atacseqy:latest \
#     bash run.sh --config config.yaml --samples samples.csv --cluster local
#
CMD ["bash"]
