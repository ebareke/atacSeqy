# 📄 Methods — atacSeqy (For Publication Supplement)

This document provides a publication-ready **Methods** section describing the atacSeqy pipeline. It can be included as Supplemental Methods for manuscripts involving ATAC-seq preprocessing, peak calling, quality control, and downstream analyses.

> **Corresponding Author:** Dr. Eric Bareke, Majewski Lab, Department of Human Genetics, McGill University

---

# 🧬 Overview of the atacSeqy Pipeline

**atacSeqy** is an automated, reproducible workflow for processing ATAC-seq (Assay for Transposase-Accessible Chromatin with sequencing) data from raw FASTQ or pre-aligned BAM files to peak-level and consensus-level chromatin accessibility profiles. The pipeline integrates industry-standard tools for alignment, quality control, peak calling, and optional downstream analyses including ArchR-based single-cell-style embeddings and chromVAR motif deviation scores.

The workflow is portable across local machines, cloud environments, and HPC clusters using SLURM or PBS job schedulers.

---

# 📥 Input Data and Metadata Structure

The pipeline accepts:
- **Paired-end or single-end FASTQ files**
- **Pre-aligned BAM files** (optional)
- A standardized **sample metadata sheet** (CSV) including sample IDs, input paths, experimental group, replicate number, and species identifier

A YAML configuration file (`config.yaml`) defines species-specific parameters including genome FASTA, BWA index prefix, effective genome size, blacklist regions, and TSS annotations.

---

# 🧹 Read Preprocessing

Read QC and preprocessing are performed using **Fastp** or **Atropos**, depending on configuration. Preprocessing includes:
- Adapter trimming
- Removal of low-quality bases
- Duplicate trimming (if enabled)
- Generation of per-sample QC reports compatible with MultiQC

---

# 🧲 Alignment and Fragment Processing

Trimmed reads are aligned to the reference genome using **BWA-MEM**, followed by:
- SAM-to-BAM conversion and sorting using **samtools**
- Removal of PCR duplicates (optional)
- Filtering of improperly paired or low-quality alignments
- Mitochondrial read quantification
- Generation of ATAC-shifted fragment BED files (Tn5-mediated shift)

---

# 📊 Quality Control Metrics

The pipeline computes multiple ATAC-seq quality control metrics:
- **Mitochondrial fraction**
- **FRiP (Fraction of Reads in Peaks)**
- **TSS enrichment score**, generated from annotated species-specific TSS regions
- **Insert size distribution** for nucleosome patterning evaluation
- **DeepTools fingerprinting** visualizations

QC outputs are automatically aggregated in a unified **MultiQC** report.

---

# 🔔 Peak Calling and Consensus Peaks

Peak calling is performed using **MACS2** with support for:
- Narrow or broad peak modes
- Species-specific effective genome sizes
- Automatic blacklist filtering

After per-sample peak calling, atacSeqy generates:
- **Union consensus peak sets** across replicates or sample groups
- Counts matrices for downstream DESeq2 analysis

---

# 📈 Differential Accessibility Analysis

If enabled, atacSeqy uses **DESeq2** to perform differential accessibility analysis:
- Peak-level normalization
- Shrinkage of fold changes
- Group comparisons defined by metadata

Results include:
- Differential peak tables
- Normalized signal matrices
- Volcano plots and quality metrics

---

# 🧬 Optional Downstream Analysis (ArchR + chromVAR)

When R dependencies are available, atacSeqy integrates:

### **ArchR**
- LSI dimensionality reduction
- Clustering and UMAP embeddings
- Gene score matrices
- Pseudo-bulk accessibility profiles

### **chromVAR**
- Motif deviation scores
- Z-score matrices for TF motif activity

These outputs facilitate integrative chromatin accessibility investigations.

---

# ⚙️ Reproducibility and Execution

The pipeline is executed using a single entry point:
```bash
bash run.sh --config config.yaml --samples samples.csv
```

The workflow supports:
- **Local execution** (single machine mode)
- **SLURM job arrays** with dependency chaining
- **PBS/Torque job arrays**

The modular design ensures full reproducibility:
- All intermediate files and logs are archived
- All parameters are recorded in JSON and YAML metadata files
- Versioned configuration files ensure experiment-level tracking

---

# 🧪 Software Versions

The following core tools are used internally (minimum versions):
- BWA ≥ 0.7
- SAMtools ≥ 1.10
- BEDTools ≥ 2.29
- MACS2 ≥ 2.2
- DeepTools ≥ 3.5
- R ≥ 4.2 (for optional modules)
- ArchR ≥ 1.0
- chromVAR ≥ 1.0

---

# 📦 Output Structure

atacSeqy generates a standardized output folder structure:
```
results/
├── alignment/
├── qc/
├── peaks/
├── consensus/
├── counts/
├── archr/          # optional
└── multiqc_report.html
```

All results are annotated with sample IDs, date stamps, and species labels.

---

# 📚 Citation

Please cite atacSeqy using the `CITATION.cff` file included in the repository.

---

# 🧾 Summary

atacSeqy provides an efficient, reproducible, species-aware pipeline for ATAC-seq analysis, automating all major preprocessing, QC, and peak-centric processing steps while offering optional advanced downstream analyses. Its flexible design allows seamless integration in lab, cloud, and HPC environments.

This Methods document is publication-ready and can be directly included in supplementary materials.

