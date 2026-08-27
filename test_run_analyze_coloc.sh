#!/bin/bash

# Configuration 
OUTPUT_DIR="/mnt/output/output/coloc"
GWAS_FP="/mnt/output/output/coloc/GP2_et_al_2025_PD_case_control_EUR_ALL_hg38_rsID.standardized.tsv"
QTL_FP="/mnt/output/output/coloc/CER_Mayo_cis_eQTL_release.standardized.tsv"
GWAS_SAMPLE_SIZE=226196
GWAS_CASE_FRACTION=0.207306937
QTL_SAMPLE_SIZE=377

# Run the analyze_colocalization.R script with the specified parameters
Rscript src/analyze_colocalization.R \
    --gwas_fp "$GWAS_FP" \
    --qtl_fp "$QTL_FP" \
    --gwas_sample_size "$GWAS_SAMPLE_SIZE" \
    --gwas_case_fraction "$GWAS_CASE_FRACTION" \
    --qtl_sample_size "$QTL_SAMPLE_SIZE" \
    --out_dir "$OUTPUT_DIR" \
    2>&1 | tee /mnt/output/output/coloc/logs/analyze_coloc.log
