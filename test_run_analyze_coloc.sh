#!/bin/bash

# Configuration 
OUTPUT_DIR="/mnt/output/output/coloc"
GWAS_FP="GP2_et_al_2025_PD_case_control_EUR_ALL_hg38_rsID.standardized.tsv"
QTL_FP="CER_Mayo_cis_eQTL_release.standardized.tsv"
QTL_STRATA_KEY="TISSUE"
GWAS_SAMPLE_SIZE=226196
GWAS_CASE_FRACTION=0.207306937
QTL_SAMPLE_SIZE=377
MIN_OVERLAP=10
MIN_SNPS_SUSIE=10

# Run the analyze_colocalization.R script with the specified parameters
Rscript src/analyze_colocalization.R \
    --gwas $GWAS_FP \
    --qtl $QTL_FP \
    --gwas_sample_size $GWAS_SAMPLE_SIZE \
    --gwas_case_fraction $GWAS_CASE_FRACTION \
    --qtl_sample_size $QTL_SAMPLE_SIZE \
    --output_dir $OUTPUT_DIR \
    --min_overlap $MIN_OVERLAP \
    --susie_min_snps $MIN_SNPS_SUSIE \
    --qtl_strata_key $QTL_STRATA_KEY \
    2>&1 | tee /mnt/output/output/coloc/logs/analyze_coloc.log
