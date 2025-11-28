# ❓ Frequently Asked Questions (FAQ) — atacSeqy

This page contains answers to the most common questions when running **atacSeqy: Autonomous ATAC-seq Processing Engine**.

> Maintainer: **Dr. Eric Bareke** (Majewski Lab, Human Genetics, McGill University)

---

# 🔍 GENERAL USAGE

## **1. What does atacSeqy do?**

It processes ATAC-seq FASTQ/BAM files end-to-end:

- Trimming & QC
- Alignment & filtering
- Mitochondrial QC, FRiP, TSS enrichment
- Shifted ATAC fragments
- MACS2 peak calling
- Consensus peaks
- MultiQC summary
- Optional **ArchR + chromVAR** analysis

---

## **2. What input formats are supported?**

atacSeqy accepts:

- **FASTQ** (single-end or paired-end)
- **BAM** (if already aligned, depending on pipeline mode)

---

## **3. Does the sample sheet require FASTQ columns if I use BAM?**

No. If using BAM input, leave `fastq1` and `fastq2` empty.

Example:

```csv
S1,,,/path/to/S1.bam,GroupA,1,human
```

---

# ⚙️ CONFIGURATION

## **4. How do I know my **``** name is correct?**

The value in the `species` column of `samples.csv` must match a key under:

```yaml
species:
  <KEY>:
```

If `species: human` exists in YAML, **species must be exactly **`` in CSV.

---

## **5. Where do I get blacklist and TSS BED files?**

Typical sources:

- ENCODE project (hg38, hg19, mm10)
- UCSC Genome Browser
- Pre-built annotation packages

ATAC-seqy expects them to match the genome assembly.

---

## **6. What does **``** mean?**

It is the **maximum allowed mitochondrial fraction**.

For example:

```yaml
mito_threshold: 0.2   # 20%
```

If samples exceed this, they may be flagged in QC.

---

# 🚀 EXECUTION

## **7. How do I test the pipeline without running anything?**

Use **dry run mode**:

```bash
bash run.sh --config config.yaml --samples samples.csv --dryrun
```

This prints all commands without executing them.

---

## **8. Does atacSeqy work on any SLURM cluster?**

Yes. It uses standard SLURM flags. Only the `--partition` and `--account` settings might need adjustments.

---

## **9. Does atacSeqy support PBS/Torque?**

Yes. Use:

```bash
--cluster pbs
```

The pipeline will generate PBS array scripts.

---

## **10. Why do I get "FASTQ missing" warnings?**

The validator (`validate_config.sh`) checks file existence. Warnings are normal during development or if using placeholder paths.

---

# 📊 QC & ANALYSIS

## **11. What QC metrics does atacSeqy compute?**

- Mitochondrial fraction
- FRiP
- TSS enrichment
- Fingerprinting (deepTools)
- Insert size distribution
- Library complexity

---

## **12. Are ArchR and chromVAR required?**

No. They are **optional**:

```yaml
archr_enabled: true
chromvar_enabled: true
```

If disabled, pipeline completes without R dependencies.

---

## **13. Do I need replicates for DESeq2?**

Yes, **at least two replicates per group**. Otherwise DESeq2 will skip comparisons.

---

# 🧪 TROUBLESHOOTING

## **14. atacSeqy crashes with “species not found”. Why?**

Your `species` column contains values missing in `config.yaml`.

Fix: Add the species in YAML or correct the sample sheet.

---

## **15. "Invalid FASTA path" error**

Check:

- Path points to a real `.fa`/`.fasta` file
- Files are readable (permissions)
- BWA index prefix matches prefix of `.fa.*` files

---

## **16. "MACS2 not found"**

Install via conda/mamba:

```bash
conda install -c bioconda macs2
```

Or run `setup.sh`.

---

# 🧬 MISCELLANEOUS

## **17. How do I cite atacSeqy?**

See the **CITATION.cff** file.

---

## **18. How do I contribute?**

Read the **CONTRIBUTING.md** page.

---

## **19. Can I request new features?**

Absolutely! Submit a GitHub Issue or PR.

---

# 🎉 More Questions?

Submit an issue here: 👉 [https://github.com/ebareke/atacSeqy/issues](https://github.com/ebareke/atacSeqy/issues)

