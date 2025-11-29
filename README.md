# 🧬 atacSeqy

### Autonomous ATAC-seq Processing Engine

A high-performance, species-aware, cluster-ready ATAC-seq pipeline engineered for **HPC**, **cloud**, and **local systems**.\
Designed to deliver **end-to-end chromatin accessibility analysis** with minimal human intervention.

---

<p align="center">
  <img src="docs/img/logo.svg" width="260" alt="atacSeqy logo" />
</p>


---

## ✨ Features

- **FASTQ → Peaks → Consensus → QC → ArchR → chromVAR → DESeq2**
- **Fully automated** execution & validation
- **Multi-species** support (Human, Mouse, Drosophila; easily extendable)
- **SLURM / PBS / Local** execution modes
- High-clarity **QC**: FRiP, TSS, Mito %, Insert Size, Fingerprinting
- **MACS2** peak calling with blacklist removal
- Optional **ArchR** project creation
- Automated **MultiQC** summarization
- Lightweight built-in **test dataset** under `tests/test_dataset/`

---

## 📁 Repository Structure

```text
atacSeqy/
├── run.sh
├── validate_config.sh
├── config.yaml
├── config.multispecies.yaml
├── tests/
│   └── test_dataset/
├── docs/
│   ├── usage.md
│   ├── configuration.md
│   ├── install.md
│   ├── faq.md
│   ├── architecture.md
│   ├── workflow.svg
│   ├── tutorial.md
│   └── methods_for_publication.md
└── .github/
    └── workflows/
        ├── ci.yml
        ├── lint.yml
        ├── docker.yml
        ├── release.yml
        └── pages.yml
```

---

## ⚙️ Installation

### Using Conda / Mamba

```bash
mamba create -n atacseqy -c conda-forge -c bioconda \
  bwa samtools bedtools macs2 deeptools yq fastp fastqc multiqc
mamba activate atacseqy
```

### Optional R ecosystem (ArchR, chromVAR, DESeq2)

```r
install.packages("BiocManager", repos = "https://cloud.r-project.org")
BiocManager::install(c("ArchR", "chromVAR", "DESeq2"))
```

### Docker

```bash
docker pull ghcr.io/ebareke/atacseqy:latest
```

Or build from the provided `Dockerfile`:

```bash
docker build -t atacseqy:latest .
```

---

## 🧬 Quick Start

```bash
bash run.sh \
  --config config.yaml \
  --samples samples.csv \
  --species human \
  --threads 16 \
  --cluster local
```

Dry-run (no execution, print commands only):

```bash
bash run.sh --config config.yaml --samples samples.csv --dryrun
```

Run on SLURM:

```bash
bash run.sh --config config.yaml --samples samples.csv --cluster slurm
```

Run on PBS:

```bash
bash run.sh --config config.yaml --samples samples.csv --cluster pbs
```

---

## 🧪 Sample Sheet Example

```csv
sample_id,fastq1,fastq2,bam,group,replicate,species
CTRL_1,fastq/C1_R1.fq.gz,fastq/C1_R2.fq.gz,,Control,1,human
TREAT_1,fastq/T1_R1.fq.gz,fastq/T1_R2.fq.gz,,Treatment,1,human
```

For BAM-based input:

```csv
sample_id,fastq1,fastq2,bam,group,replicate,species
BAM_1,,,bams/S1.bam,Control,1,human
```

---

## 📊 Pipeline Outputs

```text
results/
├── alignment/
├── fragments/
├── peaks/
├── consensus/
├── qc/
├── archr/          # optional
└── multiqc/
```

Includes:

- ATAC-shifted fragments
- MACS2 peak calls
- QC metrics (Mito %, FRiP, TSS)
- Consensus peak sets
- Peak count matrices
- Optional **ArchR** project with embeddings and motif deviations

---

## 🔮 Workflow Overview

You can reference the SVG workflow diagram:

md
![atacSeqy Workflow](docs/workflow.svg)


This illustrates the flow from FASTQ/BAM → QC → Alignment → Fragments → Peaks → Consensus → QC Aggregation → ArchR/chromVAR.

---

## ☑️ Validation Before Running

Use the validator to catch mistakes early:

```bash
bash validate_config.sh config.yaml samples.csv
```

It checks:

- YAML syntax
- Sample sheet integrity
- Species matching between CSV and config
- Existence of input files (FASTQ/BAM)

---

## 📚 Documentation

Full documentation is provided in the `docs/` directory:

- `docs/usage.md` — Usage guide
- `docs/configuration.md` — Configuration handbook (YAML & samples)
- `docs/install.md` — Installation instructions
- `docs/faq.md` — Frequently asked questions
- `docs/architecture.md` — Internal architecture
- `docs/tutorial.md` — End-to-end example
- `docs/methods_for_publication.md` — Publication-ready Methods

---

## 🧠 Development & Contribution

Contributions are welcome!

- See `docs/usage.md` for guidelines.
- See `CHANGELOG.md` for release history.
- Open issues and PRs at: [https://github.com/ebareke/atacSeqy/issues](https://github.com/ebareke/atacSeqy/issues)

---

## 📚 Citation

Please cite **atacSeqy** using the `CITATION.cff` file included in the repository.

> **Dr. Eric Bareke**\
> Majewski Lab, Department of Human Genetics, McGill University

---

## 🤝 License

This project is released under the **MIT License**.\
See `LICENSE` for full details.

---

## 🛰️ Support

For questions, issues, or feature requests, please open a GitHub Issue:\
👉 [https://github.com/ebareke/atacSeqy/issues](https://github.com/ebareke/atacSeqy/issues)

