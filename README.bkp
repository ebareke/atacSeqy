# 🌌 **atacSeqy — Autonomous ATAC‑seq Processing Engine**

---

## 🚀 **Overview**

**atacSeqy** is a next‑generation, fully automated **ATAC‑seq processing pipeline** designed for modern genomics labs, HPC environments, and large multi-sample cohorts.

It processes **FASTQ or BAM** files all the way to:

- Consensus peaks
- QC metrics (mito %, FRiP, TSS enrichment, fingerprints)
- Normalized bigWigs
- Differential accessibility
- ArchR + chromVAR analysis (optional)
- UMAP / PCA embeddings

All powered by one script: **run.sh**.

---

## ✨ **Key Features**

- ✔️ Autonomous end‑to‑end ATAC‑seq processing
- ✔️ AI‑optimized QC thresholds
- ✔️ Paired-end + single-end auto-detection
- ✔️ SLURM & PBS HPC array support
- ✔️ YAML species configuration (genome, blacklist, TSS, mito, peak mode)
- ✔️ MACS2 adaptive peak calling
- ✔️ MultiQC summary report
- ✔️ ArchR + chromVAR integrations
- ✔️ Futuristic visual identity & infographics

---

## 🧬 **Repository Structure**

```
atacSeqy/
├── run.sh                  # Main pipeline engine
├── config.yaml             # Species configuration
├── samples.csv             # Input sample sheet
├── CITATION.cff            # Citation metadata
├── LICENSE                 # MIT license
├── CONTRIBUTING.md         # Contribution guidelines
├── docs/
│   └── atacseqy-banner.svg # Futuristic repository banner
└── results/                # Outputs generated after run
```

---

## 🔧 **Prerequisites**

### System

- Linux (Ubuntu/CentOS/RHEL)
- Bash ≥ 4.0
- 50–500 GB storage recommended

### Tools

| Category     | Tools                                    |
| ------------ | ---------------------------------------- |
| Alignment    | `bwa`, `samtools`                        |
| Filtering    | `bedtools`, `awk`, `grep`                |
| QC           | `deepTools`, `multiqc`                   |
| Peak calling | `macs2`                                  |
| Counting     | `featureCounts`                          |
| YAML parsing | `yq`                                     |
| Optional (R) | `ArchR`, `chromVAR`, `DESeq2`, `ggplot2` |

---

## 📥 **Input Files**

### YAML Config (`config.yaml`)

Defines species:

- genome FASTA
- blacklist
- TSS BED
- mito chromosome
- peak mode (narrow/broad)

### Sample Sheet (`samples.csv`)

```
sample_id,fastq1,fastq2,group,replicate,species
S1,S1_R1.fq.gz,S1_R2.fq.gz,Control,1,human
```

---

## ▶️ **How to Run**

### Local

```bash
bash run.sh --config config.yaml --samples samples.csv --threads 16 --cluster local
```

### SLURM

```bash
bash run.sh --config config.yaml --samples samples.csv --threads 16 --cluster slurm
```

Automatically generates & submits:

- Array job → one sample per node
- QC + ArchR job (dependency chain)

### PBS

```bash
bash run.sh --config config.yaml --samples samples.csv --threads 16 --cluster pbs
```

---

## 📤 **Output Overview**

```
results/
├── samples/SAMPLE/
│   ├── aligned.bam
│   ├── aligned.bw
│   ├── macs2/
│   └── fragments.bedpe
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

## 📊 **QC Metrics**

| Metric         | Meaning                         |
| -------------- | ------------------------------- |
| FRiP           | Fraction of reads in peaks      |
| Mito Ratio     | % reads in mitochondrial contig |
| TSS Enrichment | Accessibility at TSS            |
| Insert Size    | Fragment length distribution    |
| Fingerprint    | Library complexity              |

---

## 🧭 **Workflow Diagram**

```
FASTQ → QC → Alignment → Filtering → ATAC Shift → Peak Calling → QC → Consensus → ArchR
```

---

## 🧪 **Example Full-Run Command**

```bash
bash run.sh \
  --config config.yaml \
  --samples samples.csv \
  --species human \
  --threads 32 \
  --cluster slurm
```

---

## 🧾 **Citation**

If you use atacSeqy, please cite:

```
Dr. Eric Bareke, Majewski Lab, Human Genetics, McGill University.
atacSeqy: Autonomous ATAC-seq Processing Engine.
GitHub: https://github.com/ebareke/atacSeqy
```

Full machine-readable version is available in **CITATION.cff**.

---

## 📜 License

This project is licensed under the **MIT License**.\
See the included `LICENSE` file.

---

## 🤝 Contributing

Guidelines are provided in `CONTRIBUTING.md`.\
We welcome:

- Bug reports
- Documentation improvements
- New species templates
- Optimization for HPC clusters

---

