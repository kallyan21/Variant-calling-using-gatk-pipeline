#!/bin/bash

# Usage: ./wes.sh <reference_genome> <sample_name> <read1.fastq.gz> <read2.fastq.gz> <known_sites.vcf>

# Ensure all required arguments are provided
if [ "$#" -ne 5 ]; then
    echo "Usage: $0 <reference_genome> <sample_name> <read1.fastq.gz> <read2.fastq.gz> <known_sites.vcf>"
    exit 1
fi

# Arguments
REFERENCE_GENOME=$1
SAMPLE_NAME=$2
READ1=$3
READ2=$4
KNOWN_SITES=$5

# Define directories
BASE_DIR=~/Research_project/WES
SAMPLE_DIR=$BASE_DIR/$SAMPLE_NAME
READS_DIR=$SAMPLE_DIR/reads
REFSEQ_DIR=$SAMPLE_DIR/refseq
ALIGNED_DIR=$READS_DIR/aligned_read
DUPLICATES_DIR=$READS_DIR/duplicates
DATA_DIR=$SAMPLE_DIR/data
RESULT_DIR=$SAMPLE_DIR/result
FILTER_DIR=$RESULT_DIR/filter
TMP_DIR=~/tmp

# Create directories if they do not exist
mkdir -p $ALIGNED_DIR $DUPLICATES_DIR $DATA_DIR $RESULT_DIR $FILTER_DIR $TMP_DIR

# --------------------------------------
# STEP 1: Map to reference using BWA-MEM
# --------------------------------------

echo "Step 1: Indexing the reference genome"
bwa index $REFERENCE_GENOME

echo "Step 1: Mapping reads to the reference genome with BWA MEM"
bwa mem -t 12 -R "@RG\tID:${SAMPLE_NAME}\tPL:ILLUMINA\tSM:${SAMPLE_NAME}" \
    $REFERENCE_GENOME \
    $READ1 \
    $READ2 \
    > $ALIGNED_DIR/${SAMPLE_NAME}.paired.sam

# -----------------------------------------
# STEP 2: Mark Duplicates and Sort - GATK4
# -----------------------------------------

echo "Step 2: Marking duplicates and sorting reads with GATK"
gatk MarkDuplicatesSpark \
    --input $ALIGNED_DIR/${SAMPLE_NAME}.paired.sam \
    --output $DUPLICATES_DIR/${SAMPLE_NAME}_sorted_dedup_reads.bam \
    --duplicate-tagging-policy OpticalOnly \
    --tmp-dir $TMP_DIR

# ----------------------------------
# STEP 3: Base quality recalibration
# ----------------------------------

echo "Step 3: Building base quality recalibration model"
gatk BaseRecalibrator \
    -I $DUPLICATES_DIR/${SAMPLE_NAME}_sorted_dedup_reads.bam \
    -R $REFERENCE_GENOME \
    --known-sites $KNOWN_SITES \
    -O $DATA_DIR/recal_data.table

echo "Step 3: Applying base quality recalibration model"
gatk ApplyBQSR \
    -I $DUPLICATES_DIR/${SAMPLE_NAME}_sorted_dedup_reads.bam \
    -R $REFERENCE_GENOME \
    --bqsr-recal-file $DATA_DIR/recal_data.table \
    -O $DUPLICATES_DIR/${SAMPLE_NAME}_sorted_dedup_bqsr_reads.bam

# -----------------------------------------------
# STEP 4: Collect Alignment & Insert Size Metrics
# -----------------------------------------------

echo "Step 4: Collecting alignment summary metrics"
gatk CollectAlignmentSummaryMetrics \
    -R $REFERENCE_GENOME \
    -I $DUPLICATES_DIR/${SAMPLE_NAME}_sorted_dedup_bqsr_reads.bam \
    -O $DUPLICATES_DIR/alignment_metrics.txt

echo "Step 4: Collecting insert size metrics"
gatk CollectInsertSizeMetrics \
    -I $DUPLICATES_DIR/${SAMPLE_NAME}_sorted_dedup_bqsr_reads.bam \
    -O $DUPLICATES_DIR/insert_size_metrics.txt \
    --HISTOGRAM_FILE $DUPLICATES_DIR/insert_size_histogram.pdf

# ----------------------------------------------
# STEP 5: Call Variants - GATK HaplotypeCaller
# ----------------------------------------------

echo "Step 5: Calling variants with GATK HaplotypeCaller"
gatk HaplotypeCaller \
    -R $REFERENCE_GENOME \
    -I $DUPLICATES_DIR/${SAMPLE_NAME}_sorted_dedup_bqsr_reads.bam \
    -O $RESULT_DIR/${SAMPLE_NAME}_raw_variants.vcf

# ----------------------------------------------
# STEP 6: Extract SNPs & INDELS
# ----------------------------------------------

echo "Step 6: Extracting SNPs"
gatk SelectVariants \
    -R $REFERENCE_GENOME \
    -V $RESULT_DIR/${SAMPLE_NAME}_raw_variants.vcf \
    --select-type SNP \
    -O $RESULT_DIR/${SAMPLE_NAME}_raw_snps.vcf

echo "Step 6: Extracting INDELs"
gatk SelectVariants \
    -R $REFERENCE_GENOME \
    -V $RESULT_DIR/${SAMPLE_NAME}_raw_variants.vcf \
    --select-type INDEL \
    -O $RESULT_DIR/${SAMPLE_NAME}_raw_indels.vcf

# --------------------------
# STEP 7: Variant filtration 
# --------------------------

echo "Step 7: Filtering SNPs"
gatk VariantFiltration \
    -R $REFERENCE_GENOME \
    -V $RESULT_DIR/${SAMPLE_NAME}_raw_snps.vcf \
    -O $FILTER_DIR/${SAMPLE_NAME}_filtered_snps.vcf \
    -filter-name "QD_filter" -filter "QD < 2.0" \
    -filter-name "FS_filter" -filter "FS > 60.0" \
    -filter-name "MQ_filter" -filter "MQ < 40.0" \
    -filter-name "SOR_filter" -filter "SOR > 4.0" \
    -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
    -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0" \
    -genotype-filter-expression "DP < 10" \
    -genotype-filter-name "DP_filter" \
    -genotype-filter-expression "GQ < 10" \
    -genotype-filter-name "GQ_filter"

echo "Step 7: Filtering INDELs"
gatk VariantFiltration \
    -R $REFERENCE_GENOME \
    -V $RESULT_DIR/${SAMPLE_NAME}_raw_indels.vcf \
    -O $FILTER_DIR/${SAMPLE_NAME}_filtered_indels.vcf \
    -filter-name "QD_filter" -filter "QD < 2.0" \
    -filter-name "FS_filter" -filter "FS > 200.0" \
    -filter-name "SOR_filter" -filter "SOR > 10.0" \
    -genotype-filter-expression "DP < 10" \
    -genotype-filter-name "DP_filter" \
    -genotype-filter-expression "GQ < 10" \
    -genotype-filter-name "GQ_filter"

echo "Step 7: Selecting variants that pass filters for SNPs"
gatk SelectVariants \
    --exclude-filtered \
    -V $FILTER_DIR/${SAMPLE_NAME}_filtered_snps.vcf \
    -O $FILTER_DIR/${SAMPLE_NAME}_analysis-ready-snps.vcf

echo "Step 7: Selecting variants that pass filters for INDELs"
gatk SelectVariants \
    --exclude-filtered \
    -V $FILTER_DIR/${SAMPLE_NAME}_filtered_indels.vcf \
    -O $FILTER_DIR/${SAMPLE_NAME}_analysis-ready-indels.vcf

echo "Step 7: Excluding variants that failed genotype filters for SNPs"
grep -v -E "DP_filter|GQ_filter" $FILTER_DIR/${SAMPLE_NAME}_analysis-ready-snps.vcf > $FILTER_DIR/${SAMPLE_NAME}_analysis-ready-snps-filteredGT.vcf

echo "Step 7: Excluding variants that failed genotype filters for INDELs"
grep -v -E "DP_filter|GQ_filter" $FILTER_DIR/${SAMPLE_NAME}_analysis-ready-indels.vcf > $FILTER_DIR/${SAMPLE_NAME}_analysis-ready-indels-filteredGT.vcf

echo "Pipeline execution completed successfully!"


