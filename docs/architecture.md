# 🧱 atacSeqy Architecture Overview

This document explains the **internal architecture**, **module design**, and **data flow** of the **atacSeqy** pipeline. It is intended for developers, contributors, and power users who want to understand or extend the pipeline.

> **Author:** Dr. Eric Bareke — Majewski Lab, Human Genetics, McGill University

---

# 1. High-Level Architecture

atacSeqy is designed around three core principles:

1. **Modularity** – each ATAC‑seq step is isolated and replaceable.
2. **Reproducibility** – all parameters are versioned via YAML and logs.
3. **Scalability** – a unified interface works on local, SLURM, or PBS clusters.

The pipeline is orchestrated by a single entry script:

```
run.sh
```

which dynamically generates tasks, arrays, dependent jobs, and post-processing modules.

---

# 2. Modular Components

The pipeline is divided into the following functional modules:

```
┌────────────────────────────────────┐
│  1. Input & Validation             │
└────────────────────────────────────┘
┌────────────────────────────────────┐
│  2. Preprocessing & QC             │
└────────────────────────────────────┘
┌────────────────────────────────────┐
│  3. Alignment & Filtering          │
└────────────────────────────────────┘
┌────────────────────────────────────┐
│  4. Fragment Derivation            │
└────────────────────────────────────┘
┌────────────────────────────────────┐
│  5. Peak Calling (MACS2)           │
└────────────────────────────────────┘
┌────────────────────────────────────┐
│  6. Consensus Building             │
└────────────────────────────────────┘
┌────────────────────────────────────┐
│  7. QC Aggregation & MultiQC        │
└────────────────────────────────────┘
┌────────────────────────────────────┐
│  8. Optional R Modules (ArchR)     │
└────────────────────────────────────┘
```

Each module runs independently in its own workspace and writes logs to:

```
results/logs/<module>/<sample>.log
```

---

# 3. Configuration Layer

All processing is controlled via two user-provided files:

### **1. config.yaml**

Defines:

- genomes and indices
- species‑specific mito behavior
- blacklist and TSS files
- MACS2 settings
- downstream module activation (ArchR, ChromVAR, DESeq2)

### **2. samples.csv**

Defines:

- sample IDs
- FASTQ/BAM paths
- group / replicate metadata
- sample species

### Configuration Validation

Before any processing starts:

```
validate_config.sh
```

checks:

- YAML formatting
- sample sheet formatting
- missing paths
- species mismatch

This prevents wasted compute time.

---

# 4. Execution Engine

The execution model supports three modes:

## **Mode A — Local Machine**

All steps run in serial or limited parallelism (`--threads`). Used for testing or small datasets.

## **Mode B — SLURM Job Arrays**

The pipeline autogenerates:

- one job array per sample
- dependency‑linked post-processing job

Example:

```
sbatch --array=1-12 align.job
sbatch --dependency=afterok:<ARRAY_JOBID> postprocess.job
```

## **Mode C — PBS / Torque**

Equivalent implementation using:

```
qsub -t 1-12
```

with automatic dependency chaining.

---

# 5. Data Flow Diagram

```
FASTQ / BAM
   │
   ▼
Preprocessing → Read QC → Trim Reports
   │
   ▼
Alignment (BWA) → Filtering (samtools)
   │
   ▼
Fragment Derivation (Tn5 shift)
   │
   ▼
Peak Calling (MACS2)
   │
   ▼
Consensus Peaks → Count Matrices
   │
   ├──► Differential Accessibility (DESeq2)
   │
   └──► ArchR Project (optional)

MultiQC aggregates QC across entire dataset
```

---

# 6. Internal File Structure

```
results/
├── aligned/            # BAM files
├── fragments/          # Tn5 shifted fragments
├── peaks/              # Per-sample peaks
├── consensus/          # Union peaks + counts
├── qc/                 # QC metrics and tables
├── logs/               # Module logs
├── multiqc/            # MultiQC HTML
└── archr/              # ArchR project directory (optional)
```

---

# 7. Logging Layer

All modules log to their own directories:

```
results/logs/<module>/<sample>.log
```

A global summary log records:

- runtime
- tool versions
- YAML configuration snapshot
- cluster job IDs (SLURM/PBS)

---

# 8. Extendability

### How to add a new module

1. Create a function block in `run.sh`
2. Add module dependencies and input/output checks
3. Register module in the pipeline order table
4. Add configuration keys (optional)
5. Document module in `docs/`

### Replaceable Components

- MACS2 may be replaced with Genrich or HMMRATAC
- BWA can be replaced with Bowtie2
- QC modules can be expanded
- ArchR can be replaced or augmented

---

# 9. Reproducibility

atacSeqy ensures reproducibility by:

- embedding configuration snapshots into results
- writing versioned logs for all tools
- strictly deterministic ordering of tasks
- avoiding non‑deterministic temporary names

This makes the pipeline publication‑ready.

---

# 10. Summary

The atacSeqy architecture is:

- **Modular** — each step isolated and replaceable
- **Scalable** — runs smoothly on laptop or HPC
- **Reproducible** — versioned parameters and logs
- **Extensible** — ideal for contributions and lab‑specific adaptations

For a full methods description, see: `docs/methods_for_publication.md`

