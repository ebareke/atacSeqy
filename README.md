# atacSeqy — ATAC-seq Orchestration Pipeline

---

## 🌌 Overview

**atacSeqy** is a fully automated, HPC‑ready **ATAC‑seq pipeline** designed for seamless chromatin accessibility analysis from FASTQs to consensus peaks, QC dashboards, and downstream ArchR/ChromVAR integration. It is driven entirely by a single Bash script: **run.sh**.

The pipeline supports:

- Mixed **paired-end** and **single-end** input
- Per‑species configuration (genome, blacklist, TSS, mito naming)
- Automated QC (mitochondrial %, FRiP, TSS enrichment, fingerprints, insert size)
- Consensus peak generation
- Differential accessibility & UMAP visualizations
- Optional ArchR project creation and motif deviation analysis

---

## 🚀 Features

- ✔️ End‑to‑end ATAC‑seq processing using one command
- ✔️ Local, PBS, or SLURM execution
- ✔️ Species-aware YAML configuration
- ✔️ Automatic FRiP / TSS Enrichment / Fingerprint QC plots
- ✔️ MultiQC integration
- ✔️ ArchR / chromVAR optional workflows
- ✔️ Reproducible output structure

---

## 📁 Repository Structure

```text
.
├── run.sh                  # Main pipeline script
├── config.yaml             # Example species configuration
├── samples.csv             # Example sample sheet
├── docs/
│   └── atacseqy-banner.svg # Futuristic pipeline banner
└── results/                # Generated after execution
```

---

## 🔧 Prerequisites

### System Requirements

- Linux (local or HPC)
- Bash ≥ 4.x
- \~50–500GB storage depending on dataset

### Dependencies

You must have the following in your environment:

| Category                  | Tools                                                |
| ------------------------- | ---------------------------------------------------- |
| **Alignment**             | `bwa`, `samtools`                                    |
| **Filtering/Shifting**    | `bedtools`, `samtools`, `awk`                        |
| **Peak Calling**          | `macs2`                                              |
| **QC**                    | `deepTools`, `multiqc`                               |
| **Counting**              | `featureCounts`                                      |
| **YAML Parsing**          | `yq`                                                 |
| **R Packages (optional)** | `ArchR`, `chromVAR`, `DESeq2`, `pheatmap`, `ggplot2` |

---

## 🧬 Input Files

### 1. YAML Configuration (`config.yaml`)

Defines species-specific parameters:

- Genome FASTA
- Blacklist file
- Mitochondrial chromosome
- Effective genome size
- Peak mode (narrow / broad)

### 2. Sample Sheet (`samples.csv`)

Example:

```csv
sample_id,fastq1,fastq2,group,replicate,species
S1,S1_R1.fq.gz,S1_R2.fq.gz,Control,1,human
S2,S2_R1.fq.gz,S2_R2.fq.gz,Treated,1,human
```

---

## ▶️ How to Run

### **Local Execution**

```bash
bash run.sh --config config.yaml --samples samples.csv --threads 16 --cluster local
```

### **SLURM Execution**

```bash
bash run.sh --config config.yaml --samples samples.csv --threads 16 --cluster slurm
```

The script automatically:

- Generates `array.slurm`
- Submits one job per sample
- Submits a dependent QC+ArchR job

### **PBS Execution**

```bash
bash run.sh --config config.yaml --samples samples.csv --threads 16 --cluster pbs
```

---

## 📤 Output Structure

```text
results/
├── samples/
│   └── SAMPLE/
│       ├── aligned.bam
│       ├── aligned.bw
│       ├── macs2/
│       └── fragments.bedpe
├── consensus/
│   ├── consensus.bed
│   └── counts.tsv
├── qc/
│   ├── frip.tsv
│   ├── mito_ratio.tsv
│   ├── fingerprint.png
│   └── tss/
└── ArchRProject/ (optional)
```

---

## 📊 Quality Control Summary

| Metric             | Description                                    |
| ------------------ | ---------------------------------------------- |
| **FRiP score**     | Fraction of reads in peaks                     |
| **Mito ratio**     | % of reads mapping to mitochondrial chromosome |
| **TSS enrichment** | Accessibility around TSS regions               |
| **Insert size**    | Fragment length QC                             |
| **Fingerprints**   | Signal-to-noise estimates                      |

---

## 🧭 Workflow Diagram (ASCII)

```
FASTQ → QC → Alignment → Filtering → ATAC Shift → Peak Calling → QC → Consensus → ArchR
```

---

## 🎨 Repository Banner (SVG)

You may place this file inside `docs/atacseqy-banner.svg`:

```svg
<svg width="1200" height="260" viewBox="0 0 1200 260" xmlns="http://www.w3.org/2000/svg">
  <rect width="100%" height="100%" fill="#0b0f19"/>
  <text x="50" y="80" fill="#66e3ff" font-family="monospace" font-size="38">atacSeqy</text>
  <text x="50" y="120" fill="#9ca3af" font-family="monospace" font-size="18">ATAC-seq Orchestration Pipeline</text>
  <circle cx="1100" cy="60" r="32" fill="#00d4ff" opacity="0.3"/>
  <circle cx="1100" cy="60" r="12" fill="#00d4ff"/>
  <rect x="50" y="160" width="1050" height="4" fill="#1e293b"/>
  <rect x="120" y="150" width="140" height="22" rx="4" fill="#2563eb"/>
  <rect x="320" y="150" width="180" height="22" rx="4" fill="#7c3aed"/>
  <rect x="550" y="150" width="200" height="22" rx="4" fill="#14b8a6"/>
  <rect x="800" y="150" width="220" height="22" rx="4" fill="#f43f5e"/>
</svg>
```

---

## 🧪 Example Command for Full Analysis

```bash
bash run.sh \
  --config config.yaml \
  --samples samples.csv \
  --species human \
  --threads 32 \
  --cluster slurm
```

---

## 📅 Roadmap

-

---

## 📬 Contact

**Author:** Eric Bareke\
For issues or contributions, please open a GitHub Issue or Pull Request.

---

