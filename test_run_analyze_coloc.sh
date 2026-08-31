#!/bin/bash

# Configuration 
OUTPUT_DIR="/mnt/output/output/coloc/amp_ad/tcx_mayo"
LD_DIR="/mnt/output/output/coloc/ld"
QC_DIR="/mnt/output/output/coloc/amp_ad/tcx_mayo/qc"
GWAS_FP="/mnt/output/output/coloc/GP2_et_al_2025_PD_case_control_EUR_ALL_hg38_rsID.standardized.tsv"
QTL_FP="/mnt/output/output/coloc/amp_ad/tcx_mayo/TCX_Mayo_cis_eQTL_release.standardized.tsv"
QTL_STRATA_KEY="REGION"
GWAS_SAMPLE_SIZE=226196
GWAS_CASE_FRACTION=0.207306937
QTL_SAMPLE_SIZE=940
MIN_OVERLAP=10
MIN_SNPS_SUSIE=10

# Run the analyze_colocalization.R script with the specified parameters
Rscript src/analyze_colocalization.R \
    --gwas_fp $GWAS_FP \
    --qtl_fp $QTL_FP \
    --gwas_sample_size $GWAS_SAMPLE_SIZE \
    --gwas_case_fraction $GWAS_CASE_FRACTION \
    --qtl_sample_size $QTL_SAMPLE_SIZE \
    --output_dir $OUTPUT_DIR \
    --ld_dir $LD_DIR \
    --qc_dir $QC_DIR \
    --min_overlap $MIN_OVERLAP \
    --susie_min_snps $MIN_SNPS_SUSIE \
    --qtl_strata_key $QTL_STRATA_KEY \
    2>&1 | tee /mnt/output/output/coloc/amp_ad/tcx_mayo/logs/analyze_coloc.log
