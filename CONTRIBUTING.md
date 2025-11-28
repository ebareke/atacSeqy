# Contributing to atacSeqy

Thank you for your interest in contributing to **atacSeqy**!  
This project powers automated, scalable ATAC-seq processing for genomics labs and HPC environments.

> Repository: https://github.com/ebareke/atacSeqy  
> Maintainer: **Dr. Eric Bareke, Majewski Lab, Human Genetics, McGill University**

---

## 🧠 How You Can Contribute

You can contribute by:
- Fixing bugs or edge cases in `run.sh`
- Improving documentation (`README`, `docs/usage.md`, `docs/configuration.md`)
- Adding new species templates (genomes, blacklists, TSS files)
- Enhancing SLURM/PBS cluster support
- Submitting new QC metrics or downstream analysis modules
- Improving workflow diagrams or repository visuals

---

## 🔁 Contribution Workflow

### 1. Fork & Clone the Repository
```bash
git clone https://github.com/ebareke/atacSeqy.git
cd atacSeqy
```

### 2. Create a Feature Branch
```bash
git checkout -b feature/my-enhancement
```

### 3. Implement Your Changes
Guidelines:
- Avoid hard-coded absolute paths
- Keep `run.sh` modular and readable
- Maintain backward compatibility for YAML schema
- Document any new CLI flags or configuration keys

### 4. Run Linting (Recommended)
If you have `shellcheck` installed:
```bash
shellcheck run.sh
```
Fix warnings where appropriate.

### 5. Update Examples
If your changes affect sample or config structure:
- Update `config.example.yaml`
- Update `samples.example.csv`

### 6. Commit Your Changes
```bash
git commit -am "Add new motif QC module and update docs"
```

### 7. Open a Pull Request
Open a PR from your fork to the main repository:
- Describe what was changed and why
- Reference issues (e.g., `Closes #12`)
- Include new dependencies if added

---

## 🐛 Reporting Bugs
When opening an issue, please include:
- OS/distribution
- Local vs HPC execution (SLURM/PBS, cluster name)
- Versions of bwa, samtools, bedtools, macs2, deepTools
- Full command used
- Relevant error messages or logs

---

## 🌱 Requesting New Features
Feature requests should include:
- Description of the feature
- Scientific or practical use-case
- Implementation ideas if applicable
- Dependencies required

---

## 🔮 Roadmap Contributions
We especially welcome contributions to:
- New species templates (human, mouse, drosophila, zebrafish...)
- Additional QC modules (entropy, duplication models, advanced FRiP plots)
- Enhanced SLURM/PBS templates
- Improved ArchR visualizations
- Containerization (Docker/Singularity)
- HTML QC dashboards

---

## 🙏 Acknowledgements
Thank you for helping expand **atacSeqy** for the scientific community. Your contributions directly support scalable, reproducible chromatin accessibility research.

