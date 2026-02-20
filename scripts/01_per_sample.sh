#!/bin/bash
#
# 01_per_sample.sh
# BRB-Seq per-sample processing pipeline
#
# Called by submit_pipeline.sh as a SLURM array job.
# Reads configuration from environment variables.
#
# For each sample:
#   1. Demultiplex by RT barcode (if DEMULTIPLEX=TRUE)
#   2. Run FastQC on demultiplexed reads
#   3. Trim adapters and polyA tails with cutadapt
#   4. Run FastQC on trimmed reads
#   5. Align with STAR
#   6. Count features with featureCounts
#   7. Post-alignment QC (RSeQC, Qualimap) - full run only
#
# Environment variables required (set by submit_pipeline.sh):
#   PROJECT_DIR, SAMPLES_FILE, READ1, READ2, DEMULTIPLEX
#   CB_START, CB_LEN, UMI_START, UMI_LEN
#   STAR_INDEX, GTF, BED, REMOVE_INTERMEDIATE, RUN_TYPE
#

set -euo pipefail

################
# Activate conda environment
################

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV:-brb_seq}"

################
# Read sample info from mapping file
################

# The mapping file now includes a header. Skip it with +1
MAPPING_LINE_READ=$(($SLURM_ARRAY_TASK_ID + 1))

read samplename barcode well < <( sed -n ${MAPPING_LINE_READ}p ${SAMPLES_FILE} )

echo "======================================"
echo "BRB-Seq Sample Processing"
echo "======================================"
echo "Sample:  ${samplename}"
echo "Barcode: ${barcode}"
echo "Well:    ${well}"
echo "Started: $(date)"
echo "======================================"
echo ""

################
# Create output directories
################

mkdir -p ${PROJECT_DIR}/logfile
mkdir -p ${PROJECT_DIR}/fastqc
mkdir -p ${PROJECT_DIR}/cutadapt
mkdir -p ${PROJECT_DIR}/STAR
mkdir -p ${PROJECT_DIR}/FeatureCounts
mkdir -p ${PROJECT_DIR}/Demultiplexed_Fastq

if [ "${RUN_TYPE}" = "full_run" ]; then
    mkdir -p ${PROJECT_DIR}/RSeQC
    mkdir -p ${PROJECT_DIR}/Qualimap
fi

################
# Initialize log file
################

LOGFILE=${PROJECT_DIR}/logfile/${samplename}_Log.txt
echo -e "Parameter:\tValue" > ${LOGFILE}
echo -e "Samplename:\t${samplename}" >> ${LOGFILE}
echo -e "Barcode:\t${barcode}" >> ${LOGFILE}
echo -e "Well:\t${well}" >> ${LOGFILE}
echo -e "Read1:\t${READ1}" >> ${LOGFILE}
echo -e "Read2:\t${READ2}" >> ${LOGFILE}
echo -e "CB_START:\t${CB_START}" >> ${LOGFILE}
echo -e "CB_LEN:\t${CB_LEN}" >> ${LOGFILE}
echo -e "UMI_START:\t${UMI_START}" >> ${LOGFILE}
echo -e "UMI_LEN:\t${UMI_LEN}" >> ${LOGFILE}
echo -e "Demultiplex:\t${DEMULTIPLEX}" >> ${LOGFILE}

################
# Demultiplex by RT Barcode (if enabled)
################

read2=""

if [ "${DEMULTIPLEX}" = "TRUE" ]; then
    echo "Demultiplexing with cutadapt..."
    
    cutadaptVersion="$(cutadapt --version 2>&1)"
    echo -e "Cutadapt_version:\t${cutadaptVersion}" >> ${LOGFILE}
    echo -e "Cutadapt_Barcode_Parameters:\t-g=^${barcode} --no-indels --minimum-length 1 --discard-untrimmed -e 0.15 -j ${SLURM_CPUS_PER_TASK}" >> ${LOGFILE}

    # Split comma-separated file lists into arrays
    IFS=',' read -r -a array1 <<< "${READ1}"
    IFS=',' read -r -a array2 <<< "${READ2}"

    for index in "${!array1[@]}"; do
        filename=$(basename -- ${array2[index]})
        filename="${filename%%.*}"
        
        echo "Processing ${filename}..."
        
        cutadapt \
            -g="^${barcode}" \
            -e 0.15 \
            --no-indels \
            --minimum-length 1 \
            -o /dev/null \
            -p ${PROJECT_DIR}/cutadapt/${samplename}_${filename}.fq.gz \
            --discard-untrimmed \
            --cores ${SLURM_CPUS_PER_TASK} \
            ${array1[$index]} ${array2[$index]} \
            > ${PROJECT_DIR}/cutadapt/${samplename}_${filename}.log
    done
    
    # Concatenate all demultiplexed files
    cat ${PROJECT_DIR}/cutadapt/${samplename}_*.fq.gz > ${PROJECT_DIR}/Demultiplexed_Fastq/${samplename}_read2.fq.gz
    rm ${PROJECT_DIR}/cutadapt/${samplename}_*.fq.gz
    
    read2=${PROJECT_DIR}/Demultiplexed_Fastq/${samplename}_read2.fq.gz
else
    read2=${READ2}
fi

# Generate MD5 checksum
md5sum ${read2} > ${PROJECT_DIR}/Demultiplexed_Fastq/${samplename}_read2.fq.gz.md5sum

################
# FastQC on input reads
################

echo "Running FastQC on input reads..."

fastqcVersion="$(fastqc --version 2>&1)"
echo -e "Fastqc_version:\t${fastqcVersion}" >> ${LOGFILE}
echo -e "Fastqc_parameters:\tDefault" >> ${LOGFILE}

fastqc -t ${SLURM_CPUS_PER_TASK} -o ${PROJECT_DIR}/fastqc/ ${read2} \
    2> ${PROJECT_DIR}/fastqc/${samplename}_fastq.stderr.txt

################
# Trim adapters with cutadapt
################

echo "Trimming adapters and short reads..."

cutadaptVersion="$(cutadapt --version 2>&1)"
echo -e "Cutadapt_Adapter_Trim_Parameters:\t--adapter=AGATCGGAAGAG --minimum-length=25 -j ${SLURM_CPUS_PER_TASK}" >> ${LOGFILE}

cutadapt \
    --adapter=AGATCGGAAGAG \
    --minimum-length=25 \
    -o ${PROJECT_DIR}/cutadapt/${samplename}_temp_trimmed.fq.gz \
    -j ${SLURM_CPUS_PER_TASK} \
    ${read2} \
    > ${PROJECT_DIR}/cutadapt/${samplename}_AdapterTrim.log

# Extract metrics from cutadapt log
numberDemultReads="$(grep "Total reads processed:" ${PROJECT_DIR}/cutadapt/${samplename}_AdapterTrim.log)"
numberDemultReads="$(echo ${numberDemultReads} | tr -dc '[:digit:]')"

numberReadsAdapter="$(grep "Reads with adapters:" ${PROJECT_DIR}/cutadapt/${samplename}_AdapterTrim.log)"
numberReadsAdapter="$(echo ${numberReadsAdapter} | sed 's/([0-9]\+\.[0-9]\+%)$//g')"
numberReadsAdapter="$(echo ${numberReadsAdapter} | tr -dc '[:digit:]')"

numberReadsWritten="$(grep "Reads written (passing filters):" ${PROJECT_DIR}/cutadapt/${samplename}_AdapterTrim.log)"
numberReadsWritten="$(echo ${numberReadsWritten} | sed 's/([0-9]\+\.[0-9]\+%)$//g')"
numberReadsWritten="$(echo ${numberReadsWritten} | tr -dc '[:digit:]')"

echo -e "Number demultiplexed reads:\t${numberDemultReads}" >> ${LOGFILE}
echo -e "Number reads with adapter:\t${numberReadsAdapter}" >> ${LOGFILE}
echo -e "Number reads written:\t${numberReadsWritten}" >> ${LOGFILE}

################
# Trim polyA tails
################

echo "Trimming polyA tails..."

cutadapt \
    --adapter="A{30}" \
    --overlap=15 \
    --minimum-length=25 \
    -o ${PROJECT_DIR}/cutadapt/${samplename}_trimmed.fq.gz \
    -j ${SLURM_CPUS_PER_TASK} \
    ${PROJECT_DIR}/cutadapt/${samplename}_temp_trimmed.fq.gz \
    > ${PROJECT_DIR}/cutadapt/${samplename}_polyA.log

# Fixed: removed stray closing parenthesis from grep
numberReadsPolyA="$(grep "Reads with adapters:" ${PROJECT_DIR}/cutadapt/${samplename}_polyA.log)"
numberReadsPolyA="$(echo ${numberReadsPolyA} | sed 's/([0-9]\+\.[0-9]\+%)//g')"
numberReadsPolyA="$(echo ${numberReadsPolyA} | tr -dc '[:digit:]')"

echo -e "Number reads with polyA:\t${numberReadsPolyA}" >> ${LOGFILE}

# Clean up temp file
rm ${PROJECT_DIR}/cutadapt/${samplename}_temp_trimmed.fq.gz

################
# FastQC on trimmed reads
################

echo "Running FastQC on trimmed reads..."

fastqc -t ${SLURM_CPUS_PER_TASK} -o ${PROJECT_DIR}/fastqc/ \
    ${PROJECT_DIR}/cutadapt/${samplename}_trimmed.fq.gz \
    2> ${PROJECT_DIR}/fastqc/${samplename}_adapterTrimmed_fastq.stderr.txt

################
# Align with STAR
################

echo "Aligning reads with STAR..."

starVersion="$(STAR --version 2>&1)"
echo -e "STAR_Version:\t${starVersion}" >> ${LOGFILE}
echo -e "STAR_Database:\t${STAR_INDEX}" >> ${LOGFILE}
echo -e "STAR_Parameters:\t--outFilterMultimapNmax 1" >> ${LOGFILE}

STAR \
    --runThreadN ${SLURM_CPUS_PER_TASK} \
    --genomeDir ${STAR_INDEX} \
    --readFilesCommand zcat \
    --readFilesIn ${PROJECT_DIR}/cutadapt/${samplename}_trimmed.fq.gz \
    --outSAMtype BAM SortedByCoordinate \
    --outFilterMultimapNmax 1 \
    --outFileNamePrefix ${PROJECT_DIR}/STAR/${samplename}_

# Extract STAR metrics
numberReadsSTAR="$(grep "Number of input reads" ${PROJECT_DIR}/STAR/${samplename}_Log.final.out)"
numberReadsSTAR="$(echo ${numberReadsSTAR} | tr -dc '[:digit:]')"

numberUniqueMap="$(grep "Uniquely mapped reads number" ${PROJECT_DIR}/STAR/${samplename}_Log.final.out)"
numberUniqueMap="$(echo ${numberUniqueMap} | tr -dc '[:digit:]')"

numberTooManyLoci="$(grep "Number of reads mapped to too many loci" ${PROJECT_DIR}/STAR/${samplename}_Log.final.out)"
numberTooManyLoci="$(echo ${numberTooManyLoci} | tr -dc '[:digit:]')"

numberTooShort="$(grep "Number of reads unmapped: too short" ${PROJECT_DIR}/STAR/${samplename}_Log.final.out)"
numberTooShort="$(echo ${numberTooShort} | tr -dc '[:digit:]')"

echo -e "Number reads input to STAR:\t${numberReadsSTAR}" >> ${LOGFILE}
echo -e "Number of reads uniquely mapped:\t${numberUniqueMap}" >> ${LOGFILE}
echo -e "Number of reads too many loci:\t${numberTooManyLoci}" >> ${LOGFILE}
echo -e "Number of reads too short:\t${numberTooShort}" >> ${LOGFILE}

################
# Count features with featureCounts
################

echo "Counting features with featureCounts..."

featureCountsVersion="$(featureCounts -v 2>&1)"
echo -e "FeatureCounts_version:\t${featureCountsVersion}" >> ${LOGFILE}
echo -e "FeatureCounts_Parameters:\t-a ${GTF}" >> ${LOGFILE}

featureCounts \
    -a ${GTF} \
    -o ${PROJECT_DIR}/FeatureCounts/${samplename}_featureCounts.txt \
    -T ${SLURM_CPUS_PER_TASK} \
    -R BAM \
    ${PROJECT_DIR}/STAR/${samplename}_Aligned.sortedByCoord.out.bam \
    2> ${PROJECT_DIR}/FeatureCounts/${samplename}_featurecounts.screen-output.log

# Create RMatrix format (gene ID + counts only)
cut -f 1,7 ${PROJECT_DIR}/FeatureCounts/${samplename}_featureCounts.txt \
    > ${PROJECT_DIR}/FeatureCounts/${samplename}_featureCounts.RMatrix.txt

# Remove first line (header)
sed -i '1d' ${PROJECT_DIR}/FeatureCounts/${samplename}_featureCounts.RMatrix.txt

# Replace second line (contains full path) with just sample name
sed -i "1c\Geneid\t${samplename}" ${PROJECT_DIR}/FeatureCounts/${samplename}_featureCounts.RMatrix.txt

################
# Post-alignment QC (full run only)
################

if [ "${RUN_TYPE}" = "full_run" ]; then
    
    echo "Running post-alignment QC..."
    
    # RSeQC
    # Fixed: typo in variable name (was rsqecVersion, now rseqcVersion)
    rseqcVersion="$(bam_stat.py --version 2>&1 || echo 'RSeQC')"
    samtoolsVersion="$(samtools --version 2>&1 | head -1)"
    
    echo -e "Rseqc_version:\t${rseqcVersion}" >> ${LOGFILE}
    echo -e "Samtools_version:\t${samtoolsVersion}" >> ${LOGFILE}
    
    # Index BAM
    samtools index ${PROJECT_DIR}/STAR/${samplename}_Aligned.sortedByCoord.out.bam \
        -@ ${SLURM_CPUS_PER_TASK} \
        ${PROJECT_DIR}/STAR/${samplename}_Aligned.sortedByCoord.out.bam.bai \
        2> /dev/null
    
    # BAM stats
    bam_stat.py -i ${PROJECT_DIR}/STAR/${samplename}_Aligned.sortedByCoord.out.bam \
        > ${PROJECT_DIR}/RSeQC/${samplename}_bam_stats_out.txt \
        2> /dev/null
    
    # Junction saturation
    junction_saturation.py -i ${PROJECT_DIR}/STAR/${samplename}_Aligned.sortedByCoord.out.bam \
        -r ${BED} \
        -o ${PROJECT_DIR}/RSeQC/${samplename}_junction_saturation \
        > /dev/null 2> /dev/null
    
    # Read distribution
    read_distribution.py -i ${PROJECT_DIR}/STAR/${samplename}_Aligned.sortedByCoord.out.bam \
        -r ${BED} \
        > ${PROJECT_DIR}/RSeQC/${samplename}_read_distribution.txt \
        2> /dev/null
    
    # Read duplication
    read_duplication.py -i ${PROJECT_DIR}/STAR/${samplename}_Aligned.sortedByCoord.out.bam \
        -o ${PROJECT_DIR}/RSeQC/${samplename}_read_duplication \
        > /dev/null 2>&1
    
    # Qualimap
    qualimapVersion="$(qualimap --version 2>&1 || echo 'QualiMap v.2.3')"
    echo -e "Qualimap_version:\t${qualimapVersion}" >> ${LOGFILE}
    
    qualimap rnaseq \
        --java-mem-size=70G \
        -outdir ${PROJECT_DIR}/Qualimap/${samplename} \
        -gtf ${GTF} \
        -bam ${PROJECT_DIR}/STAR/${samplename}_Aligned.sortedByCoord.out.bam \
        2> /dev/null
fi

################
# Clean up intermediate files
################

if [ "${REMOVE_INTERMEDIATE}" = "TRUE" ]; then
    rm -f ${PROJECT_DIR}/STAR/${samplename}_SJ.out.tab
    rm -f ${PROJECT_DIR}/STAR/${samplename}_Aligned.sortedByCoord.out.bam
    rm -f ${PROJECT_DIR}/STAR/${samplename}_Log.progress.out
    rm -f ${PROJECT_DIR}/FeatureCounts/${samplename}_Aligned.sortedByCoord.out.bam.featureCounts.bam
    rm -f ${PROJECT_DIR}/FeatureCounts/${samplename}_featureCounts.txt
fi

################
# Complete
################

echo ""
echo "======================================"
echo "Sample processing complete: ${samplename}"
echo "Finished: $(date)"
echo "======================================"
