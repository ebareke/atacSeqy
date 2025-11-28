# 📜 CHANGELOG — atacSeqy

All notable changes to **atacSeqy** will be documented in this file.  
This project follows a simple chronological changelog format.

---

## [Unreleased]
### Added
- Placeholder section for upcoming features.
- Ready for future contributions.

### Changed
- No changes yet.

### Fixed
- No fixes yet.

---

## [1.0.0] — 2025-01-01
### 🎉 Initial Stable Release
The first official release of **atacSeqy: Autonomous ATAC-seq Processing Engine**.

#### Added
- Full end‑to‑end ATAC‑seq workflow (FASTQ/BAM → consensus peaks → QC → optional ArchR/chromVAR)
- Species-aware YAML configuration system
- Mitochondrial QC, FRiP, TSS‑enrichment, fingerprinting metrics
- MACS2 peak calling with narrow/broad modes
- Consensus peak generation
- Normalized bigWigs (RPKM, CPM)
- MultiQC summary generation
- SLURM & PBS cluster mode with job arrays
- Post‑processing pipeline (`--resume`) for aggregated QC and ArchR integration
- Support for DESeq2 differential accessibility
- Example config files (minimal & multi-species)
- Validator script (`validate_config.sh`)
- Futuristic workflow diagram (SVG)
- Project documentation under `docs/`
- Project logo & banner

#### Changed
- N/A (First release)

#### Fixed
- N/A (First release)

---

## Legend
- **Added**: New features
- **Changed**: Updates or improvements
- **Fixed**: Bug fixes and patches

---

For older versions or future updates, this file will continue to serve as the authoritative change history for atacSeqy.

