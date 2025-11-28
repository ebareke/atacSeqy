# ⚙️ Configuration Guide — atacSeqy

This document describes how to configure **atacSeqy**, including:
- The YAML configuration file (`config.yaml`)
- The sample sheet (`samples.csv`)
- Multi-species setups
- Best practices for reproducible ATAC‑seq analysis

> **Author:** Dr. Eric Bareke  
> **Lab:** Majewski Lab, Department of Human Genetics, McGill University

---

# 1. YAML Configuration (`config.yaml`)

The YAML configuration controls all species‑specific and global settings for atacSeqy.

A minimal example looks like:

```yaml
species:
  human:
    genome_fa: /path/to/hg38.fa
    bwa_index_prefix: /path/to/hg38
    effective_genome_size: 2913022398
    blacklist: /path/to/hg38.blacklist.bed
    tss_bed: /path/to/hg38_tss.bed
    mito_name: chrM
    mito_filter: true
    mito_threshold: 0.2
    peak_mode: narrow

analysis:
  archr_enabled: true
  chromvar_enabled: true
  idr_enabled: false
  deseq2_enabled: true

defaults:
  species: human
  threads: 16
  outdir: results
```

---

# 2. Species Configuration Block

Each species appears under the `species:` key.

### Required Fields

| Field | Description |
|-------|-------------|
| `genome_fa` | Path to genome FASTA |
| `bwa_index_prefix` | BWA index prefix |
| `effective_genome_size` | Required by MACS2 |
| `blacklist` | BED file of blacklist regions |
| `tss_bed` | BED file of TSS sites (for TSS enrichment QC) |
| `mito_name` | Mito chromosome name |

### Optional / QC Fields

| Field | Meaning |
|-------|---------|
| `mito_filter` | If true, mitochondrial QC is applied |
| `mito_threshold` | Max mitochondrial fraction |
| `peak_mode` | `narrow`, `broad`, or `auto` |

---

# 3. Analysis Block

Controls downstream modules.

```yaml
analysis:
  archr_enabled: true
  chromvar_enabled: true
  idr_enabled: false
  deseq2_enabled: true
```

### Module Roles

- **ArchR** → LSI, embeddings, clustering, gene scores
- **chromVAR** → Motif deviation scores
- **IDR** → Replicate peak consistency
- **DESeq2** → Differential accessibility

---

# 4. Defaults Block

Defines fallback values:

```yaml
defaults:
  species: human
  threads: 16
  outdir: results
```

| Key | Meaning |
|------|---------|
| `species` | Used if no species column in samples.csv |
| `threads` | Threads per task |
| `outdir` | Main output directory |

---

# 5. Sample Sheet (`samples.csv`)

This file describes all input samples and their metadata.

### Required Columns

```csv
sample_id,fastq1,fastq2,bam,group,replicate,species
```

### Example (paired‑end)

```csv
sample_id,fastq1,fastq2,bam,group,replicate,species
CTRL_1,/fastq/C1_R1.fq.gz,/fastq/C1_R2.fq.gz,,Control,1,human
CTRL_2,/fastq/C2_R1.fq.gz,/fastq/C2_R2.fq.gz,,Control,2,human
TREAT_1,/fastq/T1_R1.fq.gz,/fastq/T1_R2.fq.gz,,Treatment,1,human
```

### Single‑end example

```csv
S1,/fastq/S1.fq.gz,,,GroupA,1,human
```

### BAM‑based input

```csv
BAM_1,,,/bams/S1.bam,GroupA,1,human
```

### Column descriptions

| Column | Meaning |
|--------|---------|
| `sample_id` | Unique name |
| `fastq1` | FASTQ R1 |
| `fastq2` | FASTQ R2 (optional) |
| `bam` | Optional pre‑aligned BAM |
| `group` | Control/Treatment/etc |
| `replicate` | Biological/technical replicate |
| `species` | Must match `species:` entry in config |

---

# 6. Multi‑Species Example

```yaml
species:
  human:
    ...
  mouse:
    ...
  drosophila:
    ...
```

And in `samples.csv`:

```csv
sample_id,fastq1,fastq2,bam,group,replicate,species
H1,/h1_R1.fq.gz,/h1_R2.fq.gz,,Control,1,human
M1,/m1_R1.fq.gz,/m1_R2.fq.gz,,Treatment,1,mouse
D1,/d1_R1.fq.gz,/d1_R2.fq.gz,,Case,1,drosophila
```

The pipeline automatically pulls the correct species block.

---

# 7. Best Practices

- Use **absolute paths** for FASTA, blacklist, TSS, and FASTQ files.  
- Keep versioned copies of production configs:  
  `config_2025-01-12.yaml`, `samples_2025-01-12.csv`
- Validate format consistency using:
  ```bash
  bash validate_config.sh config.yaml samples.csv
  ```
- Ensure blacklist and TSS BED files match genome assembly.
- Use correct mitochondrial names per species.

---

# 8. Troubleshooting

### "Species not found"
Value in `samples.csv` does not match keys under `species:`.

### "FASTA not found"
`genome_fa` points to a nonexistent file.

### "Unexpected mitochondrial fraction"
Check:
- mito chromosome naming
- blacklist correctness
- sequencing quality

---

# 9. Summary

This configuration system enables:
- Multi-species workflows  
- HPC‑scale reproducibility  
- Modularity and clarity for ATAC-seq experiments

For usage instructions see: `docs/usage.md`.  
For workflow diagram see: `docs/workflow.svg`.

