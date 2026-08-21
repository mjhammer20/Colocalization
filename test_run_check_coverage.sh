#!/bin/bash

# Configuration
LD_OUTPUT_DIR="/mnt/output/output/coloc/ld"
LD_MANIFEST="ld_manifest.tsv"
STANDARD_LOC_KEY="LOCUS_ID"
STANDARD_LEFT_KEY="LEFT_500KB"
STANDARD_RIGHT_KEY="RIGHT_500KB"
STANDARD_SNP_KEY="SNP"
STANDARD_CHR_KEY="CHR"
STANDARD_POS_KEY="BP"
STANDARD_VAR_KEY="VAR"
STANDARD_P_KEY="P"
STANDARD_EFFECT_ALLELE_KEY="EFFECT_ALLELE"
STANDARD_NON_EFFECT_ALLELE_KEY="NON_EFFECT_ALLELE"
MANIFEST_BIM_KEY="BIM"
GWAS_FP="/mnt/output/output/coloc/GP2_et_al_2025_PD_case_control_EUR_ALL_hg38_rsID.standardized.tsv"
GWAS_STRATA_KEY=""
QTL_FP="/mnt/output/output/coloc/CER_Mayo_cis_eQTL_release.standardized.tsv"
QTL_STRATA_KEY=""
HIGH_OVERLAP_MIN=0.90
MEDIUM_OVERLAP_MIN=0.70
QC_OUTPUT_DIR="/mnt/output/output/coloc/qc"

# Run the check_coverage.py script with the specified parameters
python3 src/check_coverage.py \
    --out_ld_dir "$LD_OUTPUT_DIR" \
    --ld_manifest "$LD_MANIFEST" \
    --standardized_locus_id_key "$STANDARD_LOC_KEY" \
    --standardized_chr_key "$STANDARD_CHR_KEY" \
    --standardized_left_bound_key "$STANDARD_LEFT_KEY" \
    --standardized_right_bound_key "$STANDARD_RIGHT_KEY" \
    --standardized_snp_key "$STANDARD_SNP_KEY" \
    --standardized_pos_key "$STANDARD_POS_KEY" \
    --standardized_var_key "$STANDARD_VAR_KEY" \
    --standardized_p_key "$STANDARD_P_KEY" \
    --standardized_effect_allele_key "$STANDARD_EFFECT_ALLELE_KEY" \
    --standardized_non_effect_allele_key "$STANDARD_NON_EFFECT_ALLELE_KEY" \
    --manifest_bim_key "$MANIFEST_BIM_KEY" \
    --gwas_fp "$GWAS_FP" \
    --qtl_fp "$QTL_FP" \
    --high_overlap_min "$HIGH_OVERLAP_MIN" \
    --medium_overlap_min "$MEDIUM_OVERLAP_MIN" \
    --out_qc_dir "$QC_OUTPUT_DIR" \
    2>&1 | tee /mnt/output/output/coloc/logs/check_coverage.log