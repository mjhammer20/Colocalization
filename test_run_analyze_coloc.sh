#!/bin/bash

# Configuration 
OUTPUT_DIR="/mnt/output/output/coloc/amp_ad/cortex_meta/1kg_full_EUR"
LD_DIR="/mnt/output/output/coloc/ld_1kg_full_EUR"
QC_DIR="/mnt/output/output/coloc/amp_ad/cortex_meta/qc/1kg_full_EUR"
GWAS_FP="/mnt/output/output/coloc/GP2_et_al_2025_PD_case_control_EUR_ALL_hg38_rsID.standardized.tsv"
QTL_FP="/mnt/output/output/coloc/amp_ad/cortex_meta/Cortex_MetaAnalysis_ROSMAP_CMC_HBCC_Mayo_cis_eQTL_release.standardized.tsv"
QTL_STRATA_KEY="REGION"
GWAS_SAMPLE_SIZE=226196
GWAS_CASE_FRACTION=0.207306937
QTL_SAMPLE_SIZE=1694
MIN_OVERLAP=10
MIN_SNPS_SUSIE=10
REPEAT_UNTIL_CONVERGED=FALSE
MAX_ITER=1000
LD_MANIFEST="ld_manifest_MHC.tsv"

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
    --susie_max_iter 1000 \
    --qtl_strata_key $QTL_STRATA_KEY \
    --susie_repeat_until_converged $REPEAT_UNTIL_CONVERGED \
    --susie_max_iter $MAX_ITER \
    --ld_manifest  $LD_MANIFEST \
    2>&1 | tee /mnt/output/output/coloc/amp_ad/cortex_meta/logs/analyze_coloc_1kg_full_EUR_MHC.log
