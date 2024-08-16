
# Define directories

REFERENCE_GENOME="/home/dbme/Research_project/WES/refseq/hg38.fa"
dbsnp="/home/dbme/Research_project/WES/refseq/Homo_sapiens_assembly38.dbsnp138.vcf"
thousand_genomes="/home/dbme/Research_project/WES/database/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf.gz"
indels="/home/dbme/Research_project/WES/database/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
omni="/home/dbme/Research_project/WES/database/resources_broad_hg38_v0_1000G_omni2.5.hg38.vcf.gz"
hapmap="/home/dbme/Research_project/WES/database/resources_broad_hg38_v0_hapmap_3.3.hg38.vcf.gz"
primary_targets="/home/dbme/Research_project/WES/database/exome_regions.bed.interval_list"
TMP_DIR=~/tmp
# Make directories

mkdir -p gatk_output
mkdir -p gatk_output/samples

# --------------------------------------
# STEP 1: Combine Sample 
# --------------------------------------

echo "Step 1: Combining Sample "
gatk GenomicsDBImport \
    -V /home/dbme/Research_project/WES/combined_vcf/AMN_1_raw_variants.g.vcf \
    -V /home/dbme/Research_project/WES/combined_vcf/AMN_2_raw_variants.g.vcf \
    -V /home/dbme/Research_project/WES/combined_vcf/AMN_4_raw_variants.g.vcf \
    -V /home/dbme/Research_project/WES/combined_vcf/AMN_6_raw_variants.g.vcf \
    -V /home/dbme/Research_project/WES/combined_vcf/AMN_10_raw_variants.g.vcf \
    -V /home/dbme/Research_project/WES/combined_vcf/AMN_13_raw_variants.g.vcf \
    -V /home/dbme/Research_project/WES/combined_vcf/AMN_14_raw_variants.g.vcf \
    -V /home/dbme/Research_project/WES/combined_vcf/AMN_15_raw_variants.g.vcf \
    -V /home/dbme/Research_project/WES/combined_vcf/AMN_16_raw_variants.g.vcf \
    -V /home/dbme/Research_project/WES/combined_vcf/AMN_17_raw_variants.g.vcf \
    -V /home/dbme/Research_project/WES/combined_vcf/AMN_18_raw_variants.g.vcf \
    -V /home/dbme/Research_project/WES/combined_vcf/AMN_19_raw_variants.g.vcf \
    -V /home/dbme/Research_project/WES/combined_vcf/AMN_20_raw_variants.g.vcf \
    -V /home/dbme/Research_project/WES/combined_vcf/AMN_A03_raw_variants.g.vcf \
    -V /home/dbme/Research_project/WES/combined_vcf/AMN_A06_raw_variants.g.vcf \
    -V /home/dbme/Research_project/WES/combined_vcf/AMN_A17_raw_variants.g.vcf \
    -V /home/dbme/Research_project/WES/combined_vcf/AMN_A35_raw_variants.g.vcf \
    -V /home/dbme/Research_project/WES/combined_vcf/AMN_A43_raw_variants.g.vcf \
    -V /home/dbme/Research_project/WES/combined_vcf/AMN_A50_raw_variants.g.vcf \
    -V /home/dbme/Research_project/WES/combined_vcf/AMN_E20_raw_variants.g.vcf \
    --genomicsdb-workspace-path gatk_output/sample_database \
    --intervals /home/dbme/Research_project/WES/database/exome_regions.bed.interval_list \
    --reader-threads 16 \
    --tmp-dir $TMP_DIR
# --------------------------------------
# STEP 1.1: Combine Sample 
# --------------------------------------
echo "Step 1.1: Combining Sample"
gatk GenotypeGVCFs \
    -R $REFERENCE_GENOME \
    -V gendb://gatk_output/sample_database \
    -O gatk_output/combined_samples.vcf \
   --tmp-dir $TMP_DIR
# --------------------------------------
# STEP 2: Select Varints
# --------------------------------------

echo " Step 2 : Select Variants " 

gatk SelectVariants \
    -R $REFERENCE_GENOME \
    -V gatk_output/combined_samples.vcf  \
    --select-type-to-include SNP \
    -O gatk_output/combined_samples_snps.vcf 


gatk SelectVariants \
    -R $REFERENCE_GENOME \
    -V gatk_output/combined_samples.vcf  \
    --select-type-to-include INDEL \
    -O gatk_output/combined_samples_indels.vcf 

# --------------------------------------
# STEP 3: Recalibaration 
# --------------------------------------

echo " Step 3.1  : Recalibration for SNPS" 


gatk VariantRecalibrator \
   -R $REFERENCE_GENOME \
   -V gatk_output/combined_samples_snps.vcf \
   -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap \
   -resource:omni,known=false,training=true,truth=false,prior=12.0 $omni \
   -resource:1000G,known=false,training=true,truth=false,prior=10.0 $thousand_genomes \
   -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnp \
   -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
   -mode SNP \
   -O gatk_output/combined_samples_snps.recal \
   --tranches-file gatk_output/combined_samples_snps.tranches \
   --rscript-file gatk_output/output.plots.R
   
   
 echo " Step 3.2  : Recalibration for INDEL" 
   
gatk VariantRecalibrator \
    -R $REFERENCE_GENOME \
    -V gatk_output/combined_samples_indels.vcf \
    -O gatk_output/combined_samples_indels.recal \
    --tranches-file gatk_output/combined_samples_indels.tranches \
    --resource:mills,known=false,training=true,truth=true,prior=12.0 $indels \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnp \
    -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
    -mode INDEL


# --------------------------------------
# STEP 4: Apply Recalibaration 
# --------------------------------------
echo " Step 4   :  Apply Recalibration" 



 gatk ApplyVQSR \
   -R $REFERENCE_GENOME \
   -V gatk_output/combined_samples_snps.vcf  \
   -O gatk_output/combined_samples_snps_recal.vcf \
   --truth-sensitivity-filter-level 99.0 \
   --tranches-file gatk_output/combined_samples_snps.tranches \
   --recal-file gatk_output/combined_samples_snps.recal \
   -mode SNP
   
 gatk ApplyVQSR \
   -R $REFERENCE_GENOME \
   -V gatk_output/combined_samples_indels.vcf  \
   -O gatk_output/combined_samples_indels_recal.vcf \
   --truth-sensitivity-filter-level 99.0 \
   --tranches-file gatk_output/combined_samples_indels.tranches \
   --recal-file gatk_output/combined_samples_indels.recal \
   -mode INDEL
   
# --------------------------------------
# STEP 5 : Evaluate Variants 
# --------------------------------------

echo " Step 5   :  Evaluate Variants "
gatk VariantEval \
    -R $REFERENCE_GENOME \
    -L $primary_targets \
    -eval gatk_output/combined_samples_snps_recal.vcf \
    -D $dbsnp \
    -O gatk_output/combined_samples_snps.eval 

gatk VariantEval \
    -R $REFERENCE_GENOME \
    -L $primary_targets \
    -eval gatk_output/combined_samples_indels_recal.vcf \
    -D $dbsnp \
    -O gatk_output/combined_samples_indels.eval 
    
# --------------------------------------
# STEP 6 : Select Final  variants
# --------------------------------------

echo " Select final variants " 

gatk SelectVariants \
    -R $REFERENCE_GENOME \
    -V gatk_output/combined_samples_snps_recal.vcf  \
    --select-type-to-include SNP \
    -O gatk_output/final_snps.vcf 

gatk SelectVariants \
    -R $REFERENCE_GENOME \
    -V gatk_output/combined_samples_indels_recal.vcf  \
    --select-type-to-include INDEL \
    -O gatk_output/final_indels.vcf 

