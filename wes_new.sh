#!/bin/bash

# Usage: ./wes_new.sh <sample_name> <read1.fastq.gz> <read2.fastq.gz>

# Ensure all required arguments are provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <sample_name> <read1.fastq.gz> <read2.fastq.gz>"
    exit 1
fi

# Arguments
SAMPLE_NAME=$1
READ1=$2
READ2=$3

# Define directories
REFERENCE_GENOME="/home/dbme/Research_project/WES/refseq/hg38.fa"
dbsnp="/home/dbme/Research_project/WES/refseq/Homo_sapiens_assembly38.dbsnp138.vcf"
thousand_genomes="/home/dbme/Research_project/WES/database/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf.gz"
mills="/home/dbme/Research_project/WES/database/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"

BASE_DIR=~/Research_project/WES
SAMPLE_DIR=$BASE_DIR/$SAMPLE_NAME
READS_DIR=$SAMPLE_DIR/reads
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

echo "Step 3.1: Building base quality recalibration model"
gatk BaseRecalibrator \
    -I $DUPLICATES_DIR/${SAMPLE_NAME}_sorted_dedup_reads.bam \
    -R $REFERENCE_GENOME \
    --known-sites $dbsnp \
    --known-sites $thousand_genomes \
    --known-sites $mills \
    -O $DATA_DIR/recal_data.table

echo "Step 3.2: Applying base quality recalibration model"
gatk ApplyBQSR \
    -I $DUPLICATES_DIR/${SAMPLE_NAME}_sorted_dedup_reads.bam \
    -R $REFERENCE_GENOME \
    --bqsr-recal-file $DATA_DIR/recal_data.table \
    -O $DUPLICATES_DIR/${SAMPLE_NAME}_sorted_dedup_bqsr_reads.bam

echo "Step 3.3: Building base quality recalibration model again"  
gatk BaseRecalibrator \
    -I $DUPLICATES_DIR/${SAMPLE_NAME}_sorted_dedup_bqsr_reads.bam \
    -R $REFERENCE_GENOME \
    --known-sites $dbsnp \
    --known-sites $thousand_genomes \
    --known-sites $mills \
    -O $DATA_DIR/Post_recal_data.table

# -----------------------------------------------
# STEP 4: Collect Alignment & Insert Size Metrics
# -----------------------------------------------

echo "Step 4.1: Collecting alignment summary metrics"
gatk CollectAlignmentSummaryMetrics \
    -R $REFERENCE_GENOME \
    -I $DUPLICATES_DIR/${SAMPLE_NAME}_sorted_dedup_bqsr_reads.bam \
    -O $DUPLICATES_DIR/alignment_metrics.txt

echo "Step 4.2: Collecting insert size metrics"
gatk CollectInsertSizeMetrics \
    -I $DUPLICATES_DIR/${SAMPLE_NAME}_sorted_dedup_bqsr_reads.bam \
    -O $DUPLICATES_DIR/insert_size_metrics.txt -H $DUPLICATES_DIR/insert_size_histogram.pdf

# Analyze the recalibration
echo "Step 4.3: Analyzing recalibration"
gatk AnalyzeCovariates \
    -before $DATA_DIR/recal_data.table \
    -after $DATA_DIR/Post_recal_data.table \
    -plots $DATA_DIR/recalibration_plots.pdf

# ----------------------------------------------
# STEP 5: Call Variants - GATK HaplotypeCaller
# ----------------------------------------------

echo "Step 5: Calling variants with GATK HaplotypeCaller"
gatk HaplotypeCaller \
    -R $REFERENCE_GENOME \
    -I $DUPLICATES_DIR/${SAMPLE_NAME}_sorted_dedup_bqsr_reads.bam \
    -O $RESULT_DIR/${SAMPLE_NAME}_raw_variants.g.vcf \
    -ERC GVCF
    
# ----------------------------------------------
# STEP 6: Copy raw variants to combined VCF directory
# ----------------------------------------------

echo "Step 6: Copying raw variants to the combined gVCF directory"

# Create directories if they do not exist
mkdir -p /home/dbme/Research_project/WES/combined_vcf
cp $RESULT_DIR/${SAMPLE_NAME}_raw_variants.g.vcf /home/dbme/Research_project/WES/combined_vcf/

echo "Pipeline execution completed successfully!"
