#!/bin/bash
#
# 03_multiplex_align.sh
# BRB-Seq multiplex alignment and counting (full run only)
#
# Called by submit_pipeline.sh for full_run mode.
# Aligns multiplexed reads using STARsolo and generates count matrices.
#
# Environment variables required (set by submit_pipeline.sh):
#   PROJECT_DIR, PROJECT_NAME, SAMPLES_FILE, READ1, READ2
#   CB_START, CB_LEN, UMI_START, UMI_LEN
#   STAR_INDEX, GTF, REMOVE_BAM
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
echo "BRB-Seq Multiplex Alignment"
echo "======================================"
echo "Project:  ${PROJECT_NAME}"
echo "Started:  $(date)"
echo "======================================"
echo ""

################
# Create output directories
################

mkdir -p ${PROJECT_DIR}/Star/
mkdir -p ${PROJECT_DIR}/cutadapt/
mkdir -p ${PROJECT_DIR}/Counts_Files/

################
# Create log file
################

LOGFILE=${PROJECT_DIR}/${PROJECT_NAME}_MultiplexAlign_Log.txt
echo -e "Parameter:\tValue" > ${LOGFILE}
echo -e "Project_Directory:\t${PROJECT_DIR}" >> ${LOGFILE}
echo -e "Project_Name:\t${PROJECT_NAME}" >> ${LOGFILE}
echo -e "Samples_File:\t${SAMPLES_FILE}" >> ${LOGFILE}
echo -e "Read1_File:\t${READ1}" >> ${LOGFILE}
echo -e "Read2_File:\t${READ2}" >> ${LOGFILE}
echo -e "CellBarcode_Start:\t${CB_START}" >> ${LOGFILE}
echo -e "CellBarcode_Length:\t${CB_LEN}" >> ${LOGFILE}
echo -e "UMI_Start:\t${UMI_START}" >> ${LOGFILE}
echo -e "UMI_Length:\t${UMI_LEN}" >> ${LOGFILE}
echo -e "STAR_Index:\t${STAR_INDEX}" >> ${LOGFILE}
echo -e "GTF:\t${GTF}" >> ${LOGFILE}
echo -e "Remove_BAM_Flag:\t${REMOVE_BAM}" >> ${LOGFILE}

################
# Create barcode whitelist from samples file
################

echo "Creating barcode whitelist..."

MY_WHITELIST=${PROJECT_DIR}/${PROJECT_NAME}_Whitelist.txt

# Skip header line, extract barcode column (column 2)
tail -n +2 ${SAMPLES_FILE} | cut -f 2 > ${MY_WHITELIST}

n_barcodes=$(wc -l < ${MY_WHITELIST})
echo "Whitelist created with ${n_barcodes} barcodes"

################
# Trim adapters with cutadapt (stage 1: Illumina adapter)
################

echo "Stage 1: Trimming Illumina adapters..."

cutadaptVersion="$(cutadapt --version 2>&1)"
echo -e "Cutadapt_version:\t${cutadaptVersion}" >> ${LOGFILE}
echo -e "Cutadapt_Adapter_Trim_Parameters:\t-A AGATCGGAAGAG --minimum-length=25 -j ${SLURM_CPUS_PER_TASK}" >> ${LOGFILE}

cutadapt \
    -A AGATCGGAAGAG \
    --minimum-length=25 \
    -o ${PROJECT_DIR}/cutadapt/${PROJECT_NAME}_temp_trimmed_read1.fq.gz \
    -p ${PROJECT_DIR}/cutadapt/${PROJECT_NAME}_temp_trimmed_read2.fq.gz \
    -j ${SLURM_CPUS_PER_TASK} \
    ${READ1} ${READ2} \
    > ${PROJECT_DIR}/cutadapt/${PROJECT_NAME}_AdapterTrim.log

# Extract metrics
numberDemultReads="$(grep "Total read pairs processed:" ${PROJECT_DIR}/cutadapt/${PROJECT_NAME}_AdapterTrim.log)"
numberDemultReads="$(echo ${numberDemultReads} | tr -dc '[:digit:]')"

numberReadsAdapter="$(grep "Read 2 with adapter:" ${PROJECT_DIR}/cutadapt/${PROJECT_NAME}_AdapterTrim.log)"
numberReadsAdapter="$(echo ${numberReadsAdapter} | sed 's/([0-9]\+\.[0-9]\+%)$//g')"
numberReadsAdapter="$(echo ${numberReadsAdapter} | tr -dc '[:digit:]')"

numberReadsWritten="$(grep "Pairs written (passing filters):" ${PROJECT_DIR}/cutadapt/${PROJECT_NAME}_AdapterTrim.log)"
numberReadsWritten="$(echo ${numberReadsWritten} | sed 's/([0-9]\+\.[0-9]\+%)$//g')"
numberReadsWritten="$(echo ${numberReadsWritten} | tr -dc '[:digit:]')"

echo -e "Number demultiplexed reads:\t${numberDemultReads}" >> ${LOGFILE}
echo -e "Number reads with adapter:\t${numberReadsAdapter}" >> ${LOGFILE}
echo -e "Number reads written:\t${numberReadsWritten}" >> ${LOGFILE}

################
# Trim polyA tails (stage 2)
################

echo "Stage 2: Trimming polyA tails..."

cutadapt \
    -A "A{30}" \
    --overlap=15 \
    --minimum-length=25 \
    -o ${PROJECT_DIR}/cutadapt/${PROJECT_NAME}_trimmed_read1.fq.gz \
    -p ${PROJECT_DIR}/cutadapt/${PROJECT_NAME}_trimmed_read2.fq.gz \
    -j ${SLURM_CPUS_PER_TASK} \
    ${PROJECT_DIR}/cutadapt/${PROJECT_NAME}_temp_trimmed_read1.fq.gz \
    ${PROJECT_DIR}/cutadapt/${PROJECT_NAME}_temp_trimmed_read2.fq.gz \
    > ${PROJECT_DIR}/cutadapt/${PROJECT_NAME}_polyA.log

numberReadsPolyA="$(grep "Read 2 with adapter:" ${PROJECT_DIR}/cutadapt/${PROJECT_NAME}_polyA.log)"
numberReadsPolyA="$(echo ${numberReadsPolyA} | sed 's/([0-9]\+\.[0-9]\+%)//g')"
numberReadsPolyA="$(echo ${numberReadsPolyA} | tr -dc '[:digit:]')"

echo -e "Number reads with polyA:\t${numberReadsPolyA}" >> ${LOGFILE}

# Clean up temp files
rm ${PROJECT_DIR}/cutadapt/${PROJECT_NAME}_temp_trimmed_read1.fq.gz
rm ${PROJECT_DIR}/cutadapt/${PROJECT_NAME}_temp_trimmed_read2.fq.gz

################
# Align with STARsolo
################

echo "Aligning with STARsolo..."

starVersion="$(STAR --version 2>&1)"
echo -e "STAR_Version:\t${starVersion}" >> ${LOGFILE}
echo -e "STAR_Database:\t${STAR_INDEX}" >> ${LOGFILE}

OUT_STEM=${PROJECT_DIR}/Star/${PROJECT_NAME}_

STAR \
    --runMode alignReads \
    --outSAMmapqUnique 60 \
    --soloType CB_UMI_Simple \
    --outSAMunmapped Within \
    --soloStrand Forward \
    --quantMode GeneCounts \
    --soloCBwhitelist ${MY_WHITELIST} \
    --soloFeatures Gene \
    --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
    --outFilterMultimapNmax 1 \
    --runThreadN ${SLURM_CPUS_PER_TASK} \
    --outBAMsortingThreadN ${SLURM_CPUS_PER_TASK} \
    --soloCBstart ${CB_START} \
    --soloCBlen ${CB_LEN} \
    --soloUMIstart ${UMI_START} \
    --soloUMIlen ${UMI_LEN} \
    --soloBarcodeReadLength 0 \
    --soloUMIdedup NoDedup Exact 1MM_All \
    --soloCellFilter None \
    --genomeDir ${STAR_INDEX} \
    --outFileNamePrefix ${OUT_STEM} \
    --readFilesIn ${PROJECT_DIR}/cutadapt/${PROJECT_NAME}_trimmed_read2.fq.gz \
                  ${PROJECT_DIR}/cutadapt/${PROJECT_NAME}_trimmed_read1.fq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate

################
# Convert MTX files to counts tables
################

echo "Converting MTX to count matrices..."

Rcode="/ref/bmlab/software/umi_dup_shared/Mtx_to_Counts.R"
echo -e "R_Code:\t${Rcode}" >> ${LOGFILE}

Rscript ${Rcode} \
    ${OUT_STEM}Solo.out/Gene/raw/ \
    ${SAMPLES_FILE} \
    ${PROJECT_DIR}/Counts_Files/ \
    >> ${LOGFILE}

################
# Clean up BAM and trimmed reads
################

if [ "${REMOVE_BAM}" = "TRUE" ]; then
    echo "Cleaning up BAM and intermediate files..."
    rm -f ${OUT_STEM}*.bam
    rm -f ${OUT_STEM}SJ.out.tab
    rm -f ${OUT_STEM}ReadsPerGene.out.tab
    rm -f ${OUT_STEM}Log.progress.out
    rm -f ${OUT_STEM}Log.out
    rm -f ${PROJECT_DIR}/cutadapt/${PROJECT_NAME}_trimmed_read1.fq.gz
    rm -f ${PROJECT_DIR}/cutadapt/${PROJECT_NAME}_trimmed_read2.fq.gz
fi

################
# Complete
################

echo ""
echo "======================================"
echo "Multiplex alignment complete"
echo "Counts: ${PROJECT_DIR}/Counts_Files/"
echo "Finished: $(date)"
echo "======================================"
