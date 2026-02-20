# Changelog

All notable changes to the BRB-Seq pipeline.

## [1.0.0] - 2026-02-20

### Added
- Complete pipeline modernization from legacy scripts
- Single YAML configuration system
- 3-column sample mapping format (SampleName, Barcode, Well)
- SLURM dependency chaining for automatic multi-step execution
- Automated repooling with Excel template generation
- Chemistry presets (Alithea, PrimeSeq, TripBRB)
- Sample name validation
- Conda environment for reproducibility
- Updated genome references (GRCm39/GRCh38, Ensembl 113)
- STAR 2.7.11b indices with sjdbOverhang 99

### Fixed
- FastQC version check expecting v0.x format
- Typo in RSeQC version variable (rsqecVersion → rseqcVersion)
- Stray closing parenthesis in cutadapt polyA grep
- Inconsistent spack paths across scripts
- Incomplete temp file cleanup
- Undefined variables in MultiQC script

### Changed
- Switched from custom logs to STAR Log.final.out for metrics
- Replaced spack with conda for package management
- Consolidated mapping file from 5-8 columns to 3 columns
- Moved read paths from per-row to global config
- Automated genome path resolution via paths.txt

### Removed
- Group column from mapping file (unused)
- Manual Excel copy-paste workflow
- Multiple separate parameter files

## [0.x] - Legacy

Previous versions using bash scripts with inline parameters.
