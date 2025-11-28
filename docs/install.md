# 🛠️ Installation Guide — atacSeqy

This guide explains how to install **atacSeqy** and its dependencies on:
- Local machines (macOS / Linux)
- Conda or mamba environments
- HPC clusters (SLURM / PBS)
- Docker or Singularity containers

> Maintainer: **Dr. Eric Bareke**, Majewski Lab, Human Genetics, McGill University

---

# 1. 📦 Requirements

## Operating System
- Linux (recommended)
- macOS (Intel or Apple Silicon)

## Required Tools
| Tool | Version | Purpose |
|------|---------|---------|
| `bash` | 4+ | Pipeline execution |
| `bwa` | 0.7+ | Alignment |
| `samtools` | 1.10+ | BAM operations |
| `bedtools` | 2.29+ | BED operations |
| `macs2` | 2.2+ | Peak calling |
| `deeptools` | 3.5+ | Fingerprinting, coverage |
| `yq` | 4+ | YAML parsing |
| `R` | 4.2+ | ArchR, chromVAR, DESeq2 (optional) |

Optional but recommended:
- MultiQC
- FastQC / Atropos / Fastp
- ArchR, chromVAR, DESeq2

---

# 2. 🚀 Quick Installation Using Conda (Recommended)

## Step 1 — Create the environment
```bash
mamba create -n atacseqy -c conda-forge -c bioconda \
  bwa samtools bedtools macs2 deeptools yq r-base r-essentials \
  fastqc fastp multiqc
```

## Step 2 — Activate the environment
```bash
mamba activate atacseqy
```

## Step 3 — Install R packages (Optional)
To enable ArchR, DESeq2, chromVAR:

```r
install.packages("BiocManager")
BiocManager::install(c("ArchR", "chromVAR", "SummarizedExperiment", "DESeq2"))
```

You’re ready to run the pipeline.

---

# 3. 🧬 Installing on an HPC Cluster

Most HPC clusters do not allow `conda activate` inside job scripts.

## Method A — Load modules

Example:
```bash
module load bwa/0.7.17
module load samtools/1.12
module load bedtools/2.30
module load macs2/2.2
module load deeptools/3.5
module load R/4.2
```

Configure the module list in your SLURM or PBS scripts.

## Method B — Use micromamba

```bash
module load micromamba
micromamba create -n atacseqy ...
```

---

# 4. 🧩 Installing With Docker

If you prefer containerization:

## Run container
```bash
docker run -it --rm \
  -v "$PWD":/data \
  ghcr.io/ebareke/atacseqy:latest
```

## Build container manually
```bash
docker build -t atacseqy .
```

---

# 5. 🧱 Installing With Singularity (HPC Safe)

```bash
singularity pull atacseqy.sif docker://ghcr.io/ebareke/atacseqy:latest
```

To run:
```bash
singularity exec atacseqy.sif bash run.sh --config config.yaml --samples samples.csv
```

---

# 6. 📁 Directory Setup

You should create the following minimal structure:
```
atacSeqy/
│── run.sh
│── config.yaml
│── samples.csv
└── results/   (auto-generated)
```

---

# 7. 🧪 Test the Installation

Verify all dependencies:

```bash
which bwa samtools bedtools macs2 yq
```

Run empty dry-run test:
```bash
bash run.sh --config config.yaml --samples samples.csv --dryrun
```

If you see commands printed with no errors, you're ready.

---

# 8. ❗ Troubleshooting Installation

### ❌ `macs2: command not found`
Install from bioconda:
```bash
mamba install -c bioconda macs2
```

### ❌ `yq: not found`
```bash
mamba install -c conda-forge yq
```

### ❌ Missing R dependencies
Open R and install via BiocManager:
```r
BiocManager::install(c("ArchR", "chromVAR", "DESeq2"))
```

### ❌ Permission errors on HPC
Use scratch space:
```
/home/$USER/projects/atacSeqy/
```

---

# 🎉 Installation Complete
You are now ready to run the pipeline.

See:
- `docs/usage.md` — Running the pipeline
- `docs/configuration.md` — Configuring species & analysis
- `docs/faq.md` — Additional help

