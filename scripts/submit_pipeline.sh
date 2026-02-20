#!/bin/bash
# =============================================================================
# submit_pipeline.sh
# Master submission script for the BRB-Seq pipeline.
#
# Usage:
#   bash submit_pipeline.sh config.yaml
#
# This script:
#   1. Parses config.yaml
#   2. Validates sample names
#   3. Auto-counts samples to set SLURM array size
#   4. Resolves reference genome paths from species
#   5. Submits all pipeline steps with SLURM dependency chaining
#
# Requirements:
#   conda activate brb_seq
# =============================================================================

set -euo pipefail

################
# Check arguments
################

if [ "$#" -lt 1 ]; then
    echo ""
    echo "ERROR: No config file provided."
    echo ""
    echo "Usage: bash submit_pipeline.sh config.yaml"
    echo ""
    exit 1
fi

CONFIG=$1

if [ ! -f "${CONFIG}" ]; then
    echo "ERROR: Config file not found: ${CONFIG}"
    exit 1
fi

################
# Activate conda environment
################

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate brb_seq

################
# Parse config.yaml with Python
################

# Use a single Python call to extract all values at once
eval $(python3 << EOF
import yaml, sys

with open("${CONFIG}") as f:
    c = yaml.safe_load(f)

# Project
print(f"PROJECT_NAME={c['project']['name']}")
print(f"PROJECT_DIR={c['project']['directory']}")
print(f"RUN_TYPE={c['project']['run_type']}")
print(f"SPECIES={c['project']['species']}")

# Input
print(f"SAMPLES_FILE={c['input']['samples_file']}")
print(f"READ1={c['input']['read1']}")
print(f"READ2={c['input']['read2']}")

# Chemistry
print(f"CHEMISTRY_PRESET={c['chemistry']['preset']}")
cb_start = c['chemistry'].get('cb_start') or ''
cb_len   = c['chemistry'].get('cb_len')   or ''
umi_start= c['chemistry'].get('umi_start') or ''
umi_len  = c['chemistry'].get('umi_len')   or ''
print(f"CB_START_CUSTOM={cb_start}")
print(f"CB_LEN_CUSTOM={cb_len}")
print(f"UMI_START_CUSTOM={umi_start}")
print(f"UMI_LEN_CUSTOM={umi_len}")

# Reference
print(f"STAR_INDEX_OVERRIDE={c['reference'].get('star_index') or ''}")
print(f"GTF_OVERRIDE={c['reference'].get('gtf') or ''}")
print(f"BED_OVERRIDE={c['reference'].get('bed') or ''}")

# Flags
print(f"DEMULTIPLEX={str(c['flags']['demultiplex']).upper()}")
print(f"REMOVE_INTERMEDIATE={str(c['flags']['remove_intermediate']).upper()}")
print(f"REMOVE_BAM={str(c['flags']['remove_bam']).upper()}")

# SLURM
print(f"CPUS={c['slurm']['cpus_per_sample']}")
print(f"MEM={c['slurm']['mem_per_sample_mb']}")
print(f"MAX_JOBS={c['slurm']['max_concurrent_jobs']}")
print(f"MAIL={c['slurm']['mail']}")
print(f"CONDA_ENV={c['slurm']['conda_env']}")

# Repooling
print(f"TARGET_READS={c['repooling']['target_reads_per_sample']}")
print(f"EVENNESS_THRESHOLD={c['repooling']['evenness_threshold']}")
EOF
)

################
# Validate project name
################

if [[ ! "${PROJECT_NAME}" =~ ^[A-Za-z0-9_]+$ ]]; then
    echo "ERROR: project.name '${PROJECT_NAME}' contains invalid characters."
    echo "  Only letters, numbers, and underscores are allowed."
    echo "  Remove any hyphens (-), dots (.), or spaces."
    exit 1
fi

################
# Validate run type and species
################

if [[ "${RUN_TYPE}" != "spike_in" && "${RUN_TYPE}" != "full_run" ]]; then
    echo "ERROR: run_type must be 'spike_in' or 'full_run'. Got: ${RUN_TYPE}"
    exit 1
fi

if [[ "${SPECIES}" != "mouse" && "${SPECIES}" != "human" ]]; then
    echo "ERROR: species must be 'mouse' or 'human'. Got: ${SPECIES}"
    exit 1
fi

################
# Validate sample names
################

echo "Validating sample names in ${SAMPLES_FILE}..."

INVALID=0
while IFS=$'\t' read -r samplename barcode well; do
    # Skip header
    [[ "${samplename}" == "SampleName" ]] && continue
    # Skip empty lines
    [[ -z "${samplename}" ]] && continue

    if [[ ! "${samplename}" =~ ^[A-Za-z0-9_]+$ ]]; then
        echo "  ERROR: Invalid sample name '${samplename}'"
        echo "         Only letters, numbers, and underscores allowed."
        echo "         Remove any hyphens (-), dots (.), or spaces."
        INVALID=1
    fi
done < "${SAMPLES_FILE}"

if [ "${INVALID}" -eq 1 ]; then
    echo ""
    echo "Fix sample names in ${SAMPLES_FILE} before submitting."
    exit 1
fi

echo "  All sample names are valid."
echo ""

################
# Auto-count samples
################

N_SAMPLES=$(tail -n +2 "${SAMPLES_FILE}" | grep -c -v '^$')
echo "Detected ${N_SAMPLES} samples in ${SAMPLES_FILE}"

if [ "${N_SAMPLES}" -eq 0 ]; then
    echo "ERROR: No samples found in ${SAMPLES_FILE}"
    exit 1
fi

################
# Resolve chemistry preset to CB/UMI parameters
################

case "${CHEMISTRY_PRESET}" in
    Alithea)
        CB_START=1; CB_LEN=14; UMI_START=15; UMI_LEN=14 ;;
    PrimeSeq)
        CB_START=1; CB_LEN=12; UMI_START=13; UMI_LEN=16 ;;
    TripBRB)
        CB_START=1; CB_LEN=16; UMI_START=17; UMI_LEN=10 ;;
    custom)
        if [[ -z "${CB_START_CUSTOM}" || -z "${CB_LEN_CUSTOM}" || \
              -z "${UMI_START_CUSTOM}" || -z "${UMI_LEN_CUSTOM}" ]]; then
            echo "ERROR: chemistry.preset is 'custom' but cb_start/cb_len/umi_start/umi_len are not all set."
            exit 1
        fi
        CB_START=${CB_START_CUSTOM}
        CB_LEN=${CB_LEN_CUSTOM}
        UMI_START=${UMI_START_CUSTOM}
        UMI_LEN=${UMI_LEN_CUSTOM} ;;
    *)
        echo "ERROR: Unknown chemistry preset '${CHEMISTRY_PRESET}'."
        echo "  Must be one of: Alithea, PrimeSeq, TripBRB, custom"
        exit 1 ;;
esac

echo "Chemistry preset: ${CHEMISTRY_PRESET}"
echo "  CB_START=${CB_START}, CB_LEN=${CB_LEN}, UMI_START=${UMI_START}, UMI_LEN=${UMI_LEN}"
echo ""

################
# Resolve reference genome paths
################

INDEX_BASE=/lts/bmlab/bmlab_1/bolszowy/Index

if [ "${SPECIES}" = "mouse" ]; then
    PATHS_FILE=${INDEX_BASE}/Mus_musculus/Ensembl/GRCm39/paths.txt
elif [ "${SPECIES}" = "human" ]; then
    PATHS_FILE=${INDEX_BASE}/Homo_sapiens/Ensembl/GRCh38/paths.txt
fi

# Use override from config if provided, otherwise read from paths.txt
if [ -n "${STAR_INDEX_OVERRIDE}" ]; then
    STAR_INDEX=${STAR_INDEX_OVERRIDE}
    GTF=${GTF_OVERRIDE}
    BED=${BED_OVERRIDE}
else
    if [ ! -f "${PATHS_FILE}" ]; then
        echo "ERROR: Reference paths file not found: ${PATHS_FILE}"
        echo "  Have the STAR indices been built? Run Build_STAR_Index.sh first."
        echo "  Or set reference.star_index and reference.gtf manually in ${CONFIG}"
        exit 1
    fi
    STAR_INDEX=$(grep "^STAR_INDEX" "${PATHS_FILE}" | cut -f2)
    GTF=$(grep "^GTF" "${PATHS_FILE}" | cut -f2)
    BED=$(grep "^BED" "${PATHS_FILE}" | cut -f2 || echo "")
fi

echo "Reference genome:"
echo "  Species:    ${SPECIES}"
echo "  STAR index: ${STAR_INDEX}"
echo "  GTF:        ${GTF}"
echo ""

################
# Validate reference paths
################

if [ -z "${STAR_INDEX}" ]; then
    echo "ERROR: STAR_INDEX is empty."
    echo "  Check formatting of ${PATHS_FILE}."
    echo "  Expected format: STAR_INDEX<TAB>/absolute/path"
    exit 1
fi

if [ ! -d "${STAR_INDEX}" ]; then
    echo "ERROR: STAR index directory does not exist:"
    echo "  ${STAR_INDEX}"
    exit 1
fi

if [ ! -f "${STAR_INDEX}/genomeParameters.txt" ]; then
    echo "ERROR: STAR index appears incomplete:"
    echo "  Missing genomeParameters.txt in ${STAR_INDEX}"
    exit 1
fi

if [ -z "${GTF}" ]; then
    echo "ERROR: GTF path is empty."
    echo "  Check formatting of ${PATHS_FILE}."
    exit 1
fi

if [ ! -f "${GTF}" ]; then
    echo "ERROR: GTF file not found:"
    echo "  ${GTF}"
    exit 1
fi

################
# Print submission summary
################

echo "======================================"
echo "BRB-Seq Pipeline Submission"
echo "======================================"
echo "Project:      ${PROJECT_NAME}"
echo "Run type:     ${RUN_TYPE}"
echo "Species:      ${SPECIES}"
echo "Samples:      ${N_SAMPLES}"
echo "Output dir:   ${PROJECT_DIR}"
echo "Config:       ${CONFIG}"
echo "======================================"
echo ""

################
# Create project directory
################

mkdir -p "${PROJECT_DIR}"

################
# Export variables for child scripts
################

export PROJECT_NAME PROJECT_DIR RUN_TYPE SPECIES
export SAMPLES_FILE READ1 READ2
export CB_START CB_LEN UMI_START UMI_LEN
export STAR_INDEX GTF BED
export DEMULTIPLEX REMOVE_INTERMEDIATE REMOVE_BAM
export CONDA_ENV TARGET_READS EVENNESS_THRESHOLD

################
# Step 1 — Per-sample processing (array job)
################

SCRIPTS_DIR="$(dirname "$0")"

JOB1=$(sbatch \
    --job-name="${PROJECT_NAME}_step1" \
    --array=1-${N_SAMPLES}%${MAX_JOBS} \
    --cpus-per-task=${CPUS} \
    --mem=${MEM} \
    --mail-type=FAIL \
    --mail-user=${MAIL} \
    --output="${PROJECT_DIR}/logs/step1_%A_%a.out" \
    --error="${PROJECT_DIR}/logs/step1_%A_%a.err" \
    --parsable \
    --export=ALL \
    ${SCRIPTS_DIR}/01_per_sample.sh)

echo "Submitted step 1 (per-sample):  job ${JOB1}  [${N_SAMPLES} tasks]"

################
# Step 2 — MultiQC aggregation
# Depends on ALL tasks in step 1 completing successfully
################

JOB2=$(sbatch \
    --job-name="${PROJECT_NAME}_step2_multiqc" \
    --cpus-per-task=2 \
    --mem=32000 \
    --mail-type=END,FAIL \
    --mail-user=${MAIL} \
    --output="${PROJECT_DIR}/logs/step2_%A.out" \
    --error="${PROJECT_DIR}/logs/step2_%A.err" \
    --dependency=afterok:${JOB1} \
    --parsable \
    --export=ALL \
    ${SCRIPTS_DIR}/02_multiqc.sh)

echo "Submitted step 2 (MultiQC):     job ${JOB2}  [depends on ${JOB1}]"

################
# Step 3 — Repooling report (spike-in) or Multiplex align (full run)
################

if [ "${RUN_TYPE}" = "spike_in" ]; then

    JOB3=$(sbatch \
        --job-name="${PROJECT_NAME}_step3_repooling" \
        --cpus-per-task=1 \
        --mem=8000 \
        --mail-type=END,FAIL \
        --mail-user=${MAIL} \
        --output="${PROJECT_DIR}/logs/step3_%A.out" \
        --error="${PROJECT_DIR}/logs/step3_%A.err" \
        --dependency=afterok:${JOB2} \
        --parsable \
        --export=ALL \
        ${SCRIPTS_DIR}/03_repooling_report.sh)

    echo "Submitted step 3 (repooling):   job ${JOB3}  [depends on ${JOB2}]"
    echo ""
    echo "Spike-in pipeline submitted. Final output will be in:"
    echo "  ${PROJECT_DIR}/MultiQC/"
    echo "  ${PROJECT_DIR}/repooling_report.tsv"

elif [ "${RUN_TYPE}" = "full_run" ]; then

    JOB3=$(sbatch \
        --job-name="${PROJECT_NAME}_step3_align" \
        --cpus-per-task=8 \
        --mem=75000 \
        --mail-type=FAIL \
        --mail-user=${MAIL} \
        --output="${PROJECT_DIR}/logs/step3_%A.out" \
        --error="${PROJECT_DIR}/logs/step3_%A.err" \
        --dependency=afterok:${JOB2} \
        --parsable \
        --export=ALL \
        ${SCRIPTS_DIR}/03_multiplex_align.sh)

    echo "Submitted step 3 (multiplex align): job ${JOB3}  [depends on ${JOB2}]"
    echo ""
    echo "Full run pipeline submitted. Final output will be in:"
    echo "  ${PROJECT_DIR}/MultiQC/"
    echo "  ${PROJECT_DIR}/Counts_Files/"

fi

################
# Write submission record
################

mkdir -p "${PROJECT_DIR}/logs"
SUBMISSION_LOG="${PROJECT_DIR}/logs/submission_$(date +%Y%m%d_%H%M%S).txt"
echo -e "Field\tValue" > ${SUBMISSION_LOG}
echo -e "Submitted\t$(date)" >> ${SUBMISSION_LOG}
echo -e "Config\t${CONFIG}" >> ${SUBMISSION_LOG}
echo -e "Project\t${PROJECT_NAME}" >> ${SUBMISSION_LOG}
echo -e "Run_Type\t${RUN_TYPE}" >> ${SUBMISSION_LOG}
echo -e "Species\t${SPECIES}" >> ${SUBMISSION_LOG}
echo -e "N_Samples\t${N_SAMPLES}" >> ${SUBMISSION_LOG}
echo -e "Job_Step1\t${JOB1}" >> ${SUBMISSION_LOG}
echo -e "Job_Step2\t${JOB2}" >> ${SUBMISSION_LOG}
echo -e "Job_Step3\t${JOB3}" >> ${SUBMISSION_LOG}
echo -e "STAR_Index\t${STAR_INDEX}" >> ${SUBMISSION_LOG}
echo -e "GTF\t${GTF}" >> ${SUBMISSION_LOG}
echo -e "Chemistry\t${CHEMISTRY_PRESET}" >> ${SUBMISSION_LOG}
echo -e "CB_START\t${CB_START}" >> ${SUBMISSION_LOG}
echo -e "CB_LEN\t${CB_LEN}" >> ${SUBMISSION_LOG}
echo -e "UMI_START\t${UMI_START}" >> ${SUBMISSION_LOG}
echo -e "UMI_LEN\t${UMI_LEN}" >> ${SUBMISSION_LOG}

echo ""
echo "Submission record saved to: ${SUBMISSION_LOG}"
