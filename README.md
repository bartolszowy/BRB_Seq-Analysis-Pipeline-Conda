# BRB-seq Analysis Pipeline

> Production-ready pipeline for BRB-seq data analysis with automated repooling and SLURM integration

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Conda](https://img.shields.io/badge/install%20with-conda-brightgreen.svg)](https://conda.io/)

## Overview

This pipeline processes BRB-seq (Bulk RNA Barcoding and sequencing) data from raw FASTQ files through to gene count matrices. It supports both spike-in QC runs and full production runs with comprehensive quality metrics.

**What is BRB-seq?** A cost-effective bulk RNA-seq method that uses combinatorial barcoding to dramatically reduce sequencing costs while maintaining high-quality transcriptome data. [Learn more](https://doi.org/10.1186/s13059-019-1671-x)

### Key Features

✅ **Single YAML configuration** — All parameters in one file, no scattered configs  
✅ **Automated repooling** — Generates Excel templates with optimal volumes for epMotion  
✅ **Chemistry presets** — Built-in support for Alithea, PrimeSeq, TripBRB protocols  
✅ **Smart sample validation** — Catches naming errors before wasting compute time  
✅ **SLURM dependency chaining** — Submit once, 3 steps run automatically  
✅ **Modern genome references** — GRCm39 (mouse), GRCh38 (human), Ensembl 113  
✅ **Reproducible environment** — Complete conda environment with locked versions

## Quick Start

```bash
# 1. Clone repository
git clone https://github.com/bartolszowy/BRB_Seq-Analysis-Pipeline-Conda.git
cd brb-seq-pipeline

# 1.5 Download Conda onto HTCF if this is your first time using Conda
srun --pty -c 2 --mem=8G -t 01:00:00 bash
eval `spack load --sh miniconda3`

# 2. Create conda environment
conda env create -f environment.yml
conda activate brb_seq

# 3. Copy and configure templates
cp config/config_spikein_example.yaml myproject_config.yaml # Edit with your paths

# 4. Run pipeline
bash scripts/submit_pipeline.sh myproject_config.yaml
```

## Pipeline Architecture

The pipeline consists of 3 SLURM-managed steps with automatic dependency chaining:

```
┌─────────────────────────────────────────────────────────────┐
│ Step 1: Per-Sample Processing (array job, 1-N samples)     │
├─────────────────────────────────────────────────────────────┤
│  • Demultiplex by RT barcode (cutadapt)                    │
│  • Quality control (FastQC on raw reads)                   │
│  • Trim adapters and polyA tails (cutadapt)                │
│  • Quality control (FastQC on trimmed reads)               │
│  • Align reads (STAR)                                       │
│  • Count features (featureCounts)                           │
│  • Post-alignment QC (RSeQC, Qualimap) — full run only     │
└─────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────┐
│ Step 2: Aggregate QC Metrics (depends on Step 1)           │
├─────────────────────────────────────────────────────────────┤
│  • Combine FeatureCounts outputs                            │
│  • Combine log files                                        │
│  • Generate MultiQC HTML report                             │
└─────────────────────────────────────────────────────────────┘
                              ↓
         ┌────────────────────┴────────────────────┐
         │ Spike-in              │ Full run          │
         ↓                       ↓                   
┌─────────────────────┐  ┌──────────────────────────┐
│ Step 3a: Repooling  │  │ Step 3b: Multiplex Align │
│ (depends on Step 2) │  │ (depends on Step 2)      │
├─────────────────────┤  ├──────────────────────────┤
│ • Calculate volumes │  │ • STARsolo alignment     │
│ • Generate Excel    │  │ • UMI deduplication      │
│ • Export epMotion   │  │ • Count matrices         │
└─────────────────────┘  └──────────────────────────┘
```

## Two Run Modes

### Spike-in Analysis
**Purpose:** Quality assessment and repooling optimization before full sequencing run

**When to use:** After library prep, before full sequencing (typically 1-10M reads total)

**Outputs:**
- `MultiQC/{project}_QCReport.html` — Quality metrics for all samples
- `repooling/repooling_report.tsv` — Per-sample metrics table
- `repooling/epmotion_export.tsv` — Direct import to epMotion liquid handler
- `repooling/repooling_template_filled.xlsx` — Interactive Excel template
- `repooling/repooling_summary.txt` — Human-readable summary with warnings

### Full Run Analysis
**Purpose:** Complete sequencing run with comprehensive QC

**When to use:** After repooling and full sequencing (typically 200-500M reads total)

**Outputs:**
- `MultiQC/{project}_QCReport.html` — Comprehensive QC including RSeQC and Qualimap
- `Counts_Files/Raw_Counts.txt` — Gene × sample raw count matrix
- `Counts_Files/Dedup_Counts.txt` — Gene × sample UMI-deduplicated count matrix ← **Use this for downstream analysis**

## Configuration

### Minimal Example

```yaml
project:
  name: "MyProject_SpikeIn"
  directory: "/scratch/user/projects/MyProject/spike_in"
  run_type: "spike_in"  # or "full_run"
  species: "mouse"       # or "human"

input:
  samples_file: "samples.tsv"
  read1: "/data/R1.fastq.gz"
  read2: "/data/R2.fastq.gz"

chemistry:
  preset: "Alithea"  # or "PrimeSeq", "TripBRB", "custom"

slurm:
  mail: "user@wustl.edu"
```

### Sample File Format

Simple 3-column tab-separated file:

```
SampleName	Barcode	Well
Sample_WT_1	AAGTAGAGAGTA	A1
Sample_WT_2	CCGTGAGCGGAG	A2
Sample_KO_1	TTGGAGTTGGAG	A3
```

**Critical:** Sample names must contain only letters, numbers, and underscores. No hyphens, dots, or spaces!

### Chemistry Presets

| Preset | Cell Barcode | UMI | Total Length |
|--------|--------------|-----|--------------|
| **Alithea** | 1-14 (14bp) | 15-28 (14bp) | 28bp |
| **PrimeSeq** | 1-12 (12bp) | 13-28 (16bp) | 28bp |
| **TripBRB** | 1-16 (16bp) | 17-26 (10bp) | 26bp |
| **custom** | User-defined | User-defined | Variable |

## Requirements

### System Requirements
- SLURM cluster with job arrays
- ~75GB RAM per sample during alignment
- ~100GB disk space for genome indices
- Conda or Mamba package manager

### Software (installed via conda)
- **STAR 2.7.11b** — RNA-seq alignment
- **cutadapt 4.9** — Adapter trimming and demultiplexing
- **featureCounts 2.0.6** — Gene-level counting
- **FastQC 0.12.1** — Quality control
- **MultiQC 1.13** — QC aggregation
- **RSeQC 5.0.3** — Post-alignment QC
- **Qualimap 2.3** — BAM QC
- **samtools 1.21** — BAM manipulation
- **R 4.4.2** — Count matrix generation
- **Python 3.12** — Configuration and repooling

See `environment.yml` for complete list.

## Installation

### 0.5 Log into HTCF

```bash
ssh YOUR_USERNAME@login.htcf.wustl.edu
```
then enter password

### 1. Clone Repository

```bash
git clone https://github.com/bartolszowy/BRB_Seq-Analysis-Pipeline-Conda.git
cd brb-seq-pipeline
```

### 1.5 Download and Install Conda onto your HTCF account

**This is for one-time setup per HCTF account**
```bash
srun --pty -c 2 --mem=8G -t 01:00:00 bash
eval `spack load --sh miniconda3`
```
### 2. Create Conda Environment

```bash
conda env create -f environment.yml
conda activate brb_seq

# Verify installation
STAR --version  # Should show 2.7.11b
cutadapt --version  # Should show 4.9
python --version  # Should show 3.12.x
```


### 3. Build Genome Indices

**This is a one-time setup per genome.** Indices take ~1-2 hours to build and require ~30GB disk space each.

#### Mouse (GRCm39)

```bash
# Create directory
mkdir -p /lts/bmlab/YOUR_USER/Index/Mus_musculus/Ensembl/GRCm39/download
cd /lts/bmlab/YOUR_USER/Index/Mus_musculus/Ensembl/GRCm39/download/

# Download genome and annotations
wget https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
wget https://ftp.ensembl.org/pub/release-113/gtf/mus_musculus/Mus_musculus.GRCm39.113.gtf.gz
gunzip *.gz

# Convert GTF to BED (required for RSeQC in full runs)
gtfToGenePred Mus_musculus.GRCm39.113.gtf temp.genePred
genePredToBed temp.genePred Mus_musculus.GRCm39.113.bed
rm temp.genePred

# Build STAR index
cd ..
mkdir -p STAR_2.7.11b

STAR \
    --runMode genomeGenerate \
    --genomeDir STAR_2.7.11b \
    --genomeFastaFiles download/Mus_musculus.GRCm39.dna.primary_assembly.fa \
    --sjdbGTFfile download/Mus_musculus.GRCm39.113.gtf \
    --sjdbOverhang 99 \
    --runThreadN 8 \
    --limitGenomeGenerateRAM 35000000000

# Create paths.txt for automatic resolution
cat > paths.txt << EOF
STAR_INDEX	$(pwd)/STAR_2.7.11b
GTF	$(pwd)/download/Mus_musculus.GRCm39.113.gtf
BED	$(pwd)/download/Mus_musculus.GRCm39.113.bed
EOF
```

#### Human (GRCh38)

```bash
# Same process for human
mkdir -p /lts/bmlab/YOUR_USER/Index/Homo_sapiens/Ensembl/GRCh38/download
cd /lts/bmlab/YOUR_USER/Index/Homo_sapiens/Ensembl/GRCh38/download/

wget https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz
gunzip *.gz

gtfToGenePred Homo_sapiens.GRCh38.113.gtf temp.genePred
genePredToBed temp.genePred Homo_sapiens.GRCh38.113.bed
rm temp.genePred

cd ..
mkdir -p STAR_2.7.11b

STAR \
    --runMode genomeGenerate \
    --genomeDir STAR_2.7.11b \
    --genomeFastaFiles download/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    --sjdbGTFfile download/Homo_sapiens.GRCh38.113.gtf \
    --sjdbOverhang 99 \
    --runThreadN 8 \
    --limitGenomeGenerateRAM 35000000000

cat > paths.txt << EOF
STAR_INDEX	$(pwd)/STAR_2.7.11b
GTF	$(pwd)/download/Homo_sapiens.GRCh38.113.gtf
BED	$(pwd)/download/Homo_sapiens.GRCh38.113.bed
EOF
```

### 4. Update Index Base Path

Edit `scripts/submit_pipeline.sh` line 173 with your index location:

```bash
INDEX_BASE=/lts/bmlab/YOUR_USER/Index
```

## Usage

### 1. Set Up Project Directory

```bash
mkdir -p /scratch/user/projects/MyProject
cd /scratch/user/projects/MyProject

# Copy pipeline scripts
cp -r ~/brb-seq-pipeline/scripts .

# Copy and edit config
cp ~/brb-seq-pipeline/config/config_spikein_example.yaml config.yaml
nano config.yaml

# Create samples file
cat > samples.tsv << EOF
SampleName	Barcode	Well
Sample_WT_1	AAGTAGAGAGTA	A1
Sample_WT_2	CCGTGAGCGGAG	A2
Sample_KO_1	TTGGAGTTGGAG	A3
EOF
```

### 2. Submit Pipeline

```bash
conda activate brb_seq
bash scripts/submit_pipeline.sh config.yaml
```

The pipeline will:
1. Validate your configuration and sample names
2. Submit all 3 steps with dependency chaining
3. Email you when complete (or if any step fails)

### 3. Monitor Progress

```bash
# Check job status
squeue -u $USER

# Watch live progress
tail -f output/logs/step1_*.out

# Check for errors
grep -i error output/logs/*.err
```

### 4. Retrieve Results

#### Spike-in outputs:
```bash
# View HTML QC report
firefox output/MultiQC/*_QCReport.html

# Check repooling recommendations
cat output/repooling/repooling_summary.txt

# Open Excel template
open output/repooling/repooling_template_filled.xlsx
```

#### Full run outputs:
```bash
# Load count matrix in R
dedup <- read.table("output/Counts_Files/Dedup_Counts.txt", 
                    header=TRUE, row.names=1, sep="\t")

# Dimensions: genes × samples
dim(dedup)

# Continue with your analysis (DESeq2, edgeR, Seurat, etc.)
```

## Troubleshooting

### Common Issues

**"Config file not found"**
```bash
# Use absolute path or run from project directory
bash /path/to/scripts/submit_pipeline.sh /path/to/config.yaml
```

**"Invalid sample name"**
```bash
# Fix sample names - only letters, numbers, underscores
# BAD:  Sample-WT-1, Sample.KO.1
# GOOD: Sample_WT_1, Sample_KO_1
```

**"paths.txt not found"**
```bash
# Verify genome indices are built
ls /lts/bmlab/YOUR_USER/Index/Mus_musculus/Ensembl/GRCm39/paths.txt

# Check file format (tab-separated)
cat paths.txt
# Should show:
# STAR_INDEX	/path/to/STAR_2.7.11b
# GTF	/path/to/genome.gtf
# BED	/path/to/genome.bed
```

**Jobs stay pending forever**
```bash
# Check why
scontrol show job <JOB_ID>

# Look for "Reason:" field
# "Priority" = waiting in queue (normal)
# "Resources" = not enough nodes available
# "DependencyNeverSatisfied" = previous step failed
```

**Step 1 task failed**
```bash
# Find which sample failed
sacct -j <STEP1_JOB_ID> --format=JobID,State,ExitCode | grep FAILED

# Check error log
cat output/logs/step1_<JOB_ID>_<TASK_ID>.err
```

**BED file error in full run**
```bash
# RSeQC needs BED file - convert GTF to BED
cd /path/to/genome/dir
gtfToGenePred genome.gtf temp.genePred
genePredToBed temp.genePred genome.bed
rm temp.genePred

# Update paths.txt
echo -e "BED\t$(pwd)/genome.bed" >> paths.txt
```

For more issues, see [docs/TROUBLESHOOTING.md](docs/TROUBLESHOOTING.md)

## Output Files

### Directory Structure

```
project_directory/
├── MultiQC/
│   └── {project}_QCReport.html          # Main QC report
├── Counts_Files/                         # Full run only
│   ├── Raw_Counts.txt                    # Gene × sample raw counts
│   └── Dedup_Counts.txt                  # Gene × sample deduplicated counts ★
├── repooling/                            # Spike-in only
│   ├── repooling_report.tsv              # Per-sample metrics
│   ├── epmotion_export.tsv               # epMotion liquid handler import
│   ├── repooling_template_filled.xlsx    # Interactive Excel template
│   └── repooling_summary.txt             # Human-readable summary
├── STAR/
│   └── {sample}_Log.final.out            # Per-sample alignment metrics
├── FeatureCounts/
│   └── {sample}_featureCounts.txt.summary  # Per-sample counting metrics
├── fastqc/
│   └── *.html                            # Per-sample FastQC reports
└── logs/
    ├── step1_*.out                       # Per-sample processing logs
    ├── step2_*.out                       # MultiQC log
    ├── step3_*.out                       # Repooling/multiplex log
    └── submission_*.txt                  # Pipeline submission record
```


**Maintained by the Muegge Lab at Washington University in St. Louis**

For questions, contact: [Lab email/Slack]
