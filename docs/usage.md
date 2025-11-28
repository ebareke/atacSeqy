# 🧬 atacSeqy — Usage Guide

This document explains how to run **atacSeqy** on local machines and HPC clusters using **SLURM** or **PBS** backends.

> Pipeline: atacSeqy  
> Author: Dr. Eric Bareke (Majewski Lab, Human Genetics, McGill University)

---

# 1. Quick Start

```bash
bash run.sh \
  --config config.yaml \
  --samples samples.csv \
  --species human \
  --threads 16 \
  --cluster local
```

---

# 2. Command-Line Options

| Flag | Description |
|------|-------------|
| `--config FILE` | YAML configuration file |
| `--samples FILE` | CSV sample sheet |
| `--species STR` | Species key (overridden per sample if column exists) |
| `--threads N` | Number of threads per task |
| `--outdir DIR` | Output directory (default: results) |
| `--cluster STR` | `local`, `slurm`, `pbs` |
| `--deseq` | Enable DESeq2 differential accessibility |
| `--idr` | Enable IDR replicate analysis |
| `--resume` | Run post-processing after arrays (consensus, QC, ArchR) |
| `--dryrun` | Print commands only |
| `--help` | Show help message |

---

# 3. Running Locally

```bash
bash run.sh \
  --config config.yaml \
  --samples samples.csv \
  --species human \
  --threads 8 \
  --cluster local
```

Use for:
- Small datasets  
- Testing  
- Development with `--dryrun`

---

# 4. Running on SLURM

```bash
bash run.sh \
  --config config.yaml \
  --samples samples.csv \
  --threads 32 \
  --cluster slurm \
  --deseq \
  --resume
```

**atacSeqy automatically:**
1. Creates a SLURM job array  
2. Submits one task per sample  
3. Submits a dependent post-processing job

See: `cluster/slurm_example_job.sh`

---

# 5. Running on PBS

```bash
bash run.sh \
  --config config.yaml \
  --samples samples.csv \
  --cluster pbs \
  --threads 16
```

See: `cluster/pbs_example_job.sh`

---

# 6. Dry Run (No Execution)

```bash
bash run.sh \
  --config config.yaml \
  --samples samples.csv \
  --cluster local \
  --dryrun
```

---

# 7. Typical Workflow

1. Prepare and test `config.yaml`  
2. Build `samples.csv`  
3. Validate configuration:
   ```bash
   bash validate_config.sh config.yaml samples.csv
   ```
4. Run locally (optional)  
5. Run on SLURM or PBS  
6. Examine `results/` and QC reports  
7. Explore ArchR project (if enabled)

---

# 8. Outputs

Main outputs include:
- Aligned BAM/BigWig files
- MACS2 peaks
- Consensus peak set
- ATAC-shifted fragments
- FRiP/TSS/Mito QC tables
- MultiQC report
- ArchRProject (optional)

See the main README for full output structure.

---

# 9. Support

If you encounter issues, open a GitHub Issue:  
👉 https://github.com/ebareke/atacSeqy/issues

