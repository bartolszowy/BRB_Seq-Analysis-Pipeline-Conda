#!/bin/bash
#
# 02_multiqc.sh
# BRB-Seq MultiQC aggregation pipeline
#
# Called by submit_pipeline.sh after all per-sample jobs complete.
# Reads configuration from environment variables.
#
# Tasks:
#   1. Combine FeatureCounts outputs into single matrix
#   2. Combine log files into single file
#   3. Generate MultiQC HTML report
#   4. Clean up intermediate QC files (if REMOVE_INTERMEDIATE=TRUE)
#
# Environment variables required (set by submit_pipeline.sh):
#   PROJECT_DIR, PROJECT_NAME, RUN_TYPE, REMOVE_INTERMEDIATE
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
echo "BRB-Seq MultiQC Aggregation"
echo "======================================"
echo "Project:  ${PROJECT_NAME}"
echo "Run type: ${RUN_TYPE}"
echo "Started:  $(date)"
echo "======================================"
echo ""

################
# Create MultiQC directory
################

mkdir -p ${PROJECT_DIR}/MultiQC/

################
# Record tool versions
################

multiqcVersion="$(multiqc --version 2>/dev/null)"
echo "MultiQC Version: ${multiqcVersion}"
echo ""

################
# Combine FeatureCounts outputs into single matrix
################

echo "Combining FeatureCounts outputs..."

FEATURE_LOCATION=${PROJECT_DIR}/${PROJECT_NAME}_${RUN_TYPE}_FeatureCounts.txt

# Paste combines all RMatrix files, awk keeps first column and every other column
paste -d '\t' ${PROJECT_DIR}/FeatureCounts/*RMatrix.txt | \
    awk -F"\t" '
        {OFS="\t"}
        {DL="";
         for (i=1;i<=NF;i+=(i<2 ? 1: 2))
             {printf "%s%s", DL, $i, DL="\t"};
         printf "\n"}
    ' > ${FEATURE_LOCATION}

echo "Combined FeatureCounts written to: ${FEATURE_LOCATION}"

################
# Combine log files
################

echo "Combining log files..."

LOG_LOCATION=${PROJECT_DIR}/${PROJECT_NAME}_${RUN_TYPE}_LogCombined.txt

paste -d '\t' ${PROJECT_DIR}/logfile/*Log.txt | \
    awk -F"\t" '
        {OFS="\t"}
        {DL="";
         for (i=1;i<=NF;i+=(i<2 ? 1: 2))
             {printf "%s%s", DL, $i, DL="\t"};
         printf "\n"}
    ' > ${LOG_LOCATION}

echo "Combined logs written to: ${LOG_LOCATION}"

################
# Generate MultiQC report
################

echo "Generating MultiQC report..."

# Determine which config to use based on run type
if [ "${RUN_TYPE}" = "spike_in" ]; then
    CONFIG_SRC="/ref/bmlab/software/BRB-Seq/multiqc_config_spikein_240116.yaml"
    CONFIG_NAME="multiqc_config_spike-in.yaml"
else
    CONFIG_SRC="/ref/bmlab/software/BRB-Seq/multiqc_config_fullrun_240116.yaml"
    CONFIG_NAME="multiqc_config_fullrun.yaml"
fi

CONFIG_TGT="${PROJECT_DIR}/MultiQC/${CONFIG_NAME}"

# Copy config and substitute project name
if [ -f "${CONFIG_SRC}" ]; then
    PROJECT_NAME=${PROJECT_NAME} envsubst < ${CONFIG_SRC} > ${CONFIG_TGT}
else
    echo "WARNING: MultiQC config template not found: ${CONFIG_SRC}"
    echo "         Proceeding without custom config."
    CONFIG_TGT=""
fi

# Build MultiQC command based on run type
MULTIQC_INPUTS="${PROJECT_DIR}/fastqc/*fastqc.zip \
    ${PROJECT_DIR}/STAR/*Log.final.out \
    ${PROJECT_DIR}/FeatureCounts/*featureCounts.txt.summary"

if [ "${RUN_TYPE}" = "full_run" ]; then
    MULTIQC_INPUTS="${MULTIQC_INPUTS} \
        ${PROJECT_DIR}/RSeQC/* \
        ${PROJECT_DIR}/Qualimap/*/rnaseq_qc_results.txt \
        ${PROJECT_DIR}/Qualimap/*/raw_data_qualimapReport/*.txt"
fi

# Run MultiQC
if [ -n "${CONFIG_TGT}" ]; then
    multiqc ${MULTIQC_INPUTS} \
        --filename "${PROJECT_DIR}/MultiQC/${PROJECT_NAME}_QCReport" \
        -v \
        -f \
        --config ${CONFIG_TGT} \
        2> ${PROJECT_DIR}/MultiQC/MultiQC_Report.log.txt
else
    multiqc ${MULTIQC_INPUTS} \
        --filename "${PROJECT_DIR}/MultiQC/${PROJECT_NAME}_QCReport" \
        -v \
        -f \
        2> ${PROJECT_DIR}/MultiQC/MultiQC_Report.log.txt
fi

echo "MultiQC report written to: ${PROJECT_DIR}/MultiQC/${PROJECT_NAME}_QCReport.html"

################
# Clean up QC intermediate files
################

if [ "${REMOVE_INTERMEDIATE}" = "TRUE" ]; then
    echo "Cleaning up intermediate QC files..."
    
    rm -f ${PROJECT_DIR}/FeatureCounts/*RMatrix.txt
    rm -f ${PROJECT_DIR}/FeatureCounts/*screen-output.log
    rm -f ${PROJECT_DIR}/STAR/*Log.out
    rm -f ${PROJECT_DIR}/STAR/*bam.bai
    rm -f ${PROJECT_DIR}/cutadapt/*fq.gz
    rm -rf ${PROJECT_DIR}/Demultiplexed_Fastq/*
    
    if [ "${RUN_TYPE}" = "full_run" ]; then
        rm -rf ${PROJECT_DIR}/RSeQC/*
        rm -rf ${PROJECT_DIR}/Qualimap/*
    fi
    
    # Keep fastqc zips — they're small and needed for MultiQC
    # rm -rf ${PROJECT_DIR}/fastqc/*
fi

################
# Complete
################

echo ""
echo "======================================"
echo "MultiQC aggregation complete"
echo "Report: ${PROJECT_DIR}/MultiQC/${PROJECT_NAME}_QCReport.html"
echo "Finished: $(date)"
echo "======================================"
