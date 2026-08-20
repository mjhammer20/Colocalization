#!/bin/bash

# Configuration
ENTREZ_EMAIL="matt@datatecnica.com"

# Input Summary Statistics
SUM_STATS_FILE="/mnt/working/summary_analysis_data/gwas_stats/GP2_et_al_2025_PD_case_control_EUR_ALL_hg38_rsID.txt.gz"
HEADER_LINES=0
SS_CHR_KEY="chromosome"
SS_POS_KEY="base_pair_position"
SS_RSID_KEY="rsID"
SS_GENE_ID_KEY=""
SS_NONEFFECT_KEY="other_allele"
SS_EFFECT_KEY="effect_allele"
SS_P_KEY="p_value"
SS_BETA_KEY="beta"
SS_SE_KEY="standard_error"
SS_N_KEY=""
SS_STATISTIC_KEY=""
SS_MAF_KEY="effect_allele_frequency"
SS_MAC_KEY=""
SS_VAR_BETA_KEY=""
SS_SDY_KEY=""
SS_TISSUE_KEY=""

# Input Loci File
OUTPUT_DIR="/mnt/output/output/coloc"
LOCI_FILE="meta6_gwas_loci.merged.tsv"
LOCI_LEFT_BOUND_KEY="LEFT_500KB"
LOCI_RIGHT_BOUND_KEY="RIGHT_500KB"

# Standardized Column Names
STANDARDIZED_CHR_KEY="CHR"
STANDARDIZED_POS_KEY="BP"
STANDARDIZED_RSID_KEY="SNP"
STANDARDIZED_VAR_ID_KEY="VAR"
STANDARDIZED_GENE_ID_KEY="GENE"
STANDARDIZED_NONEFFECT_KEY="NON_EFFECT"
STANDARDIZED_EFFECT_KEY="EFFECT"
STANDARDIZED_P_KEY="P"
STANDARDIZED_BETA_KEY="BETA"
STANDARDIZED_SE_KEY="SE"
STANDARDIZED_N_KEY="N"
STANDARDIZED_STATISTIC_KEY="STAT"
STANDARDIZED_MAF_KEY="MAF"
STANDARDIZED_VAR_BETA_KEY="VARBETA"
STANDARDIZED_SDY_KEY="SDY"
STANDARDIZED_TISSUE_KEY="TISSUE"

# Run sum_stats_transform.py
python3 -u src/sum_stats_transform.py \
    --sum_stats_file $SUM_STATS_FILE \
    --loci_file $LOCI_FILE \
    --output_dir $OUTPUT_DIR \
    --header_lines $HEADER_LINES \
    --entrez_email $ENTREZ_EMAIL \
    --ss_chr_key $SS_CHR_KEY \
    --ss_pos_key $SS_POS_KEY \
    --ss_rsid_key $SS_RSID_KEY \
    --ss_non_effect_allele_key $SS_NONEFFECT_KEY \
    --ss_effect_allele_key $SS_EFFECT_KEY \
    --ss_p_key $SS_P_KEY \
    --ss_beta_key $SS_BETA_KEY \
    --ss_se_key $SS_SE_KEY \
    --ss_maf_key $SS_MAF_KEY \
    --standardized_chr_key $STANDARDIZED_CHR_KEY \
    --standardized_pos_key $STANDARDIZED_POS_KEY \
    --standardized_rsid_key $STANDARDIZED_RSID_KEY \
    --standardized_variant_id_key $STANDARDIZED_VAR_ID_KEY \
    --standardized_gene_id_key $STANDARDIZED_GENE_ID_KEY \
    --standardized_non_effect_allele_key $STANDARDIZED_NONEFFECT_KEY \
    --standardized_effect_allele_key $STANDARDIZED_EFFECT_KEY \
    --standardized_p_key $STANDARDIZED_P_KEY \
    --standardized_beta_key $STANDARDIZED_BETA_KEY \
    --standardized_se_key $STANDARDIZED_SE_KEY \
    --standardized_n_key $STANDARDIZED_N_KEY \
    --standardized_statistic_key $STANDARDIZED_STATISTIC_KEY \
    --standardized_maf_key $STANDARDIZED_MAF_KEY \
    --standardized_var_beta_key $STANDARDIZED_VAR_BETA_KEY \
    --standardized_sdy_key $STANDARDIZED_SDY_KEY \
    --standardized_tissue_key $STANDARDIZED_TISSUE_KEY \
    --loci_left_bound_key $LOCI_LEFT_BOUND_KEY \
    --loci_right_bound_key $LOCI_RIGHT_BOUND_KEY \
    2>&1 | tee /mnt/output/output/coloc/test_run_sum_stats_transform_gwas.log