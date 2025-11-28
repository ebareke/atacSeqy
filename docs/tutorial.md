# 📘 atacSeqy Tutorial — End‑to‑End Example

This tutorial walks you through a complete **ATAC‑seq analysis** using **atacSeqy** — from raw FASTQ files to peak calling, consensus peaks, QC summarization, and optional ArchR analysis.

> **Author:** Dr. Eric Bareke — Majewski Lab, Human Genetics, McGill University

---

# 1. 📁 Prepare Your Working Directory

Create a project folder:
```
mkdir atacseqy_tutorial
cd atacseqy_tutorial
```

Recommended structure:
```
atacseqy_tutorial/
├── run.sh
├── config.yaml
├── samples.csv
└── results/          # created automatically
```

Download or copy the `run.sh` script into this folder.

---

# 2. 🔧 Install Dependencies
See `docs/install.md` for full instructions. For quick setup:

```bash
mamba create -n atacseqy -c conda-forge -c bioconda \
  bwa samtools bedtools macs2 deeptools yq fastp fastqc multiqc
mamba activate atacseqy
```

Optional for downstream R analysis:
```r
BiocManager::install(c("ArchR", "chromVAR", "DESeq2"))
```

---

# 3. 🧬 Add FASTQ Files

Place your FASTQ files in a directory, for example:
```
fastq/
├── sampleA_R1.fastq.gz
├── sampleA_R2.fastq.gz
├── sampleB_R1.fastq.gz
└── sampleB_R2.fastq.gz
```

You may also use BAM files instead of FASTQ.

---

# 4. 📝 Create a Sample Sheet (`samples.csv`)

Your metadata file should look like this:

```csv
sample_id,fastq1,fastq2,bam,group,replicate,species
A,/path/fastq/sampleA_R1.fastq.gz,/path/fastq/sampleA_R2.fastq.gz,,Control,1,human
B,/path/fastq/sampleB_R1.fastq.gz,/path/fastq/sampleB_R2.fastq.gz,,Treatment,1,human
```

If using BAM inputs:
```csv
S1,,,/bams/S1.bam,Control,1,human
```

---

# 5. ⚙️ Create a Minimal Config File (`config.yaml`)

Example:
```yaml
species:
  human:
    genome_fa: /path/to/hg38.fa
    bwa_index_prefix: /path/to/hg38
    blacklist: /path/to/hg38.blacklist.bed
    tss_bed: /path/to/hg38_tss.bed
    effective_genome_size: 2913022398
    mito_name: chrM
    mito_filter: true
    mito_threshold: 0.2
    peak_mode: narrow

defaults:
  species: human
  threads: 16
  outdir: results
```

Multi‑species users: see `config.multispecies.yaml`.

---

# 6. 🔍 Validate the Configuration

Run:
```bash
bash validate_config.sh config.yaml samples.csv
```

This checks:
- YAML validity
- sample sheet format
- FASTQ/BAM existence
- species matching

No critical errors → proceed.

---

# 7. 🚀 Run atacSeqy

To run locally:
```bash
bash run.sh --config config.yaml --samples samples.csv --threads 16 --cluster local
```

To test without running anything:
```bash
bash run.sh --config config.yaml --samples samples.csv --dryrun
```

To run on SLURM:
```bash
bash run.sh --config config.yaml --samples samples.csv --cluster slurm
```

To run on PBS:
```bash
bash run.sh --config config.yaml --samples samples.csv --cluster pbs
```

---

# 8. 📂 Understanding Output Files

The pipeline generates:
```
results/
├── aligned/         # BAM + index
├── fragments/       # Tn5 shifted fragments
├── qc/              # mito, FRiP, TSS, insert size
├── peaks/           # per-sample MACS2 peaks
├── consensus/       # union peaks + counts
├── multiqc/         # HTML summary
└── archr/           # optional ArchR project
```

Key outputs:
- `*_fragments.bed.gz` — shifted fragments
- `*_peaks.narrowPeak` — MACS2 peaks
- `consensus/peak_matrix.tsv` — count matrix
- `results/multiqc_report.html` — global quality overview

---

# 9. 📈 Explore QC

Open MultiQC report:
```
results/multiqc/multiqc_report.html
```

Check per-sample:
- mitochondrial fraction
- TSS enrichment
- Insert size distribution
- FRiP scores

---

# 10. 🧬 Optional: Run ArchR

If R modules are installed and ArchR enabled in config:
```bash
Rscript archr_build.R
```

ArchR outputs include:
- UMAP embeddings
- LSI clusters
- Gene score matrices
- Motif deviation matrices

---

# 11. 🧪 Optional: Differential Accessibility (DESeq2)

Enable in config:
```yaml
deseq2_enabled: true
```

Outputs:
- Differential peak table
- MA/volcano plots
- Normalized matrices

---

# 12. 🎉 Summary

After this tutorial, you can now:
- Prepare metadata and configuration files
- Run atacSeqy locally or on a cluster
- Interpret QC and peak-level outputs
- Optionally run ArchR and DESeq2

For deeper details:
- Architecture: `docs/architecture.md`
- Configuration: `docs/configuration.md`
- Installation: `docs/install.md`
- FAQ: `docs/faq.md`

You are now ready to process full ATAC‑seq datasets using **atacSeqy**!

