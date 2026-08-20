#!/bin/bash

# Configuration
LD_OUTPUT_DIR="/mnt/output/output/coloc/ld"
LD_MANIFEST="ld_manifest.tsv"
MANIFEST_LOC_KEY="LOCUS_ID"
MANIFEST_CHR_KEY="CHR"
MANIFEST_LEFT_KEY="LEFT_500KB"
MANIFEST_RIGHT_KEY="RIGHT_500KB"
MANIFEST_BIM_KEY="BIM"
GWAS_FP="/mnt/working/summary_analysis_data/gwas_stats/GP2_et_al_2025_PD_case_control_EUR_ALL_hg38_rsID.txt.gz"
GWAS_STRATA_KEY=""
QTL_FP="/mnt/output/output/coloc/CER_Mayo_cis_eQTL_release.standardized.tsv"
QTL_STRATA_KEY=""
STANDARD_SNP_KEY="SNP"
STANDARD_CHR_KEY="CHR"
STANDARD_POS_KEY="BP"
STANDARD_VAR_KEY="VAR"
STANDARD_P_COL="P"
HIGH_OVERLAP_MIN=0.90
MEDIUM_OVERLAP_MIN=0.70
QC_OUTPUT_DIR="/mnt/output/output/coloc/qc"