#!/bin/bash
#
# 03_repooling_report.sh
# Wrapper for automated repooling calculation (spike-in only)
#
# Called by submit_pipeline.sh after MultiQC completes.
# Calls the Python repooling script with appropriate arguments.
#
# Environment variables required (set by submit_pipeline.sh):
#   PROJECT_DIR, SAMPLES_FILE, TARGET_READS, EVENNESS_THRESHOLD
#

set -euo pipefail

################
# Activate conda environment
################

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV:-brb_seq}"

################
# Header
################

echo "======================================"
echo "BRB-Seq Repooling Report"
echo "======================================"
echo "Project:      ${PROJECT_NAME}"
echo "Target reads: ${TARGET_READS}M"
echo "Started:      $(date)"
echo "======================================"
echo ""

################
# Find the Python script
################

# Look in common locations
SCRIPT_LOCATIONS=(
    "$(dirname $0)/03_repooling_report.py"
    "${PROJECT_DIR}/scripts/03_repooling_report.py"
    "/ref/bmlab/software/BRB-Seq/03_repooling_report.py"
)

REPOOLING_SCRIPT=""
for loc in "${SCRIPT_LOCATIONS[@]}"; do
    if [ -f "${loc}" ]; then
        REPOOLING_SCRIPT="${loc}"
        break
    fi
done

if [ -z "${REPOOLING_SCRIPT}" ]; then
    echo "ERROR: Could not find 03_repooling_report.py"
    echo "Searched:"
    for loc in "${SCRIPT_LOCATIONS[@]}"; do
        echo "  ${loc}"
    done
    exit 1
fi

echo "Using repooling script: ${REPOOLING_SCRIPT}"
echo ""

################
# Run the repooling calculation
################

python3 ${REPOOLING_SCRIPT} \
    --project-dir ${PROJECT_DIR} \
    --samples-file ${SAMPLES_FILE} \
    --target-reads ${TARGET_READS} \
    --output-dir ${PROJECT_DIR}/repooling \
    --output-excel

################
# Complete
################

echo ""
echo "======================================"
echo "Repooling report complete"
echo "Outputs:"
echo "  ${PROJECT_DIR}/repooling/repooling_report.tsv"
echo "  ${PROJECT_DIR}/repooling/epmotion_export.tsv"
echo "  ${PROJECT_DIR}/repooling/repooling_summary.txt"
echo "  ${PROJECT_DIR}/repooling/repooling_template_filled.xlsx"
echo "Finished: $(date)"
echo "======================================"
