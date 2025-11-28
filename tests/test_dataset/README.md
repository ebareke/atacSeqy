# Test Dataset for atacSeqy

This directory contains a **minimal synthetic dataset** for testing the atacSeqy pipeline.

It is NOT biologically meaningful and is intended only for:

- CI / GitHub Actions tests
- Local dry-run tests
- Quick smoke tests after installation

## Contents

- `config.yaml` — minimal configuration pointing to `genome.fa`
- `samples.csv` — two tiny paired-end samples (A and B)
- `genome.fa` — toy reference genome (one small contig)
- `fastq/` — very small paired-end FASTQ files (2 reads per file)

## Recommended Usage

Dry run:

```bash
bash run.sh \
  --config tests/test_dataset/config.yaml \
  --samples tests/test_dataset/samples.csv \
  --cluster local \
  --dryrun
