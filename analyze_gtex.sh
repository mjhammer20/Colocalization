#!/bin/bash

# Conda Initialization
eval "$(conda shell.bash hook)"

# Standardization Parameters
ENTREZ_EMAIL="matt@datatecnica.com"
QTL_SUM_STATS_GENOME_BUILD="hg38"
HEADER_LINES=0
SS_CHR_KEY="chr"
SS_POS_KEY="variant_pos"
SS_RSID_KEY="rs_id_dbSNP155_GRCh38p13"
SS_GENE_ID_KEY="gene_name"
SS_NONEFFECT_KEY="ref"
SS_EFFECT_KEY="alt"
SS_P_KEY="pval_nominal"
SS_BETA_KEY="slope"
SS_SE_KEY="slope_se"
SS_STATISTIC_KEY=""
SS_MAF_KEY="af"
SS_MAC_KEY=""
SS_VAR_BETA_KEY=""
SS_SDY_KEY=""
SS_STRATA_KEY="tissue"

# Input Loci File Parameters
LOCI_FILE="/mnt/output/output/coloc/meta6_gwas_loci.merged.tsv"
LOCI_LEFT_BOUND_KEY="LEFT_500KB"
LOCI_RIGHT_BOUND_KEY="RIGHT_500KB"

# GWAS Summary Statistics File Parameters
STANDARDIZED_GWAS_SUM_STATS_FP="/mnt/output/output/coloc/GP2_et_al_2025_PD_case_control_EUR_ALL_hg38_rsID.standardized.tsv"
GWAS_SAMPLE_SIZE=226196
GWAS_CASE_FRACTION=0.207306937

# Standardized Summary Statistics Parameters
STANDARDIZED_CHR_KEY="CHR"
STANDARDIZED_POS_KEY="BP"
STANDARDIZED_RSID_KEY="SNP"
STANDARDIZED_VAR_ID_KEY="VAR"
STANDARDIZED_GENE_ID_KEY="GENE"
STANDARDIZED_NON_EFFECT_KEY="NON_EFFECT"
STANDARDIZED_EFFECT_KEY="EFFECT"
STANDARDIZED_P_KEY="P"
STANDARDIZED_BETA_KEY="BETA"
STANDARDIZED_SE_KEY="SE"
STANDARDIZED_STATISTIC_KEY="STAT"
STANDARDIZED_MAF_KEY="MAF"
STANDARDIZED_VAR_BETA_KEY="VARBETA"
STANDARDIZED_STRATA_KEY="TISSUE"

# LD Manifest File Parameters
LD_OUTPUT_DIR="/mnt/output/output/coloc/ld"
LD_MANIFEST="ld_manifest.tsv"
MANIFEST_LOC_KEY="LOCUS_ID"
MANIFEST_BIM_KEY="BIM"
HIGH_OVERLAP_MIN=0.90
MEDIUM_OVERLAP_MIN=0.70

# Susie Parameters
MIN_OVERLAP=10
MIN_SNPS_SUSIE=10

# GTex Tissue Manifest File Parameters
GTEX_MANIFEST_FILE="/mnt/working/locus_reports/qtl/eqtl/gtex/significant/GTEx_Analysis_v10_eQTL_updated/GTEx_tissue_manifest.tsv"

# Loop through each tissue in the GTEx manifest and run workflow
while IFS=$'\t' read -r index tissue_name sample_size tissue_file; do

    # Activate conda environment
    conda activate omics-base
    
    # Extract tissue information
    TISSUE_NAME=${tissue_name}
    QTL_SUM_STATS_FILE=${tissue_file}
    QTL_SS_N=${sample_size}

    # Log Status
    echo "Processing tissue: $TISSUE_NAME"

    # Output Directories
    OUTPUT_DIR="/mnt/output/output/coloc/gtex/${TISSUE_NAME}"
    QC_OUTPUT_DIR="/mnt/output/output/coloc/gtex/${TISSUE_NAME}/qc"
    LOGS_OUTPUT_DIR="/mnt/output/output/coloc/gtex/${TISSUE_NAME}/logs"

    # Create output directories if they don't exist
    mkdir -p $OUTPUT_DIR
    mkdir -p $QC_OUTPUT_DIR
    mkdir -p $LOGS_OUTPUT_DIR

    # Run sum_stats_transform.py
    # python3 -u src/sum_stats_transform.py \
    #     --sum_stats_file $QTL_SUM_STATS_FILE \
    #     --sum_stats_genome_build $QTL_SUM_STATS_GENOME_BUILD \
    #     --loci_file $LOCI_FILE \
    #     --output_dir $OUTPUT_DIR \
    #     --qc_output_dir $QC_OUTPUT_DIR \
    #     --entrez_email $ENTREZ_EMAIL \
    #     --ss_chr_key $SS_CHR_KEY \
    #     --ss_pos_key $SS_POS_KEY \
    #     --ss_rsid_key $SS_RSID_KEY \
    #     --ss_gene_id_key $SS_GENE_ID_KEY \
    #     --ss_non_effect_allele_key $SS_NONEFFECT_KEY \
    #     --ss_effect_allele_key $SS_EFFECT_KEY \
    #     --ss_p_key $SS_P_KEY \
    #     --ss_beta_key $SS_BETA_KEY \
    #     --ss_se_key $SS_SE_KEY \
    #     --ss_maf_key $SS_MAF_KEY \
    #     --ss_n $QTL_SS_N \
    #     --ss_strata_key $SS_STRATA_KEY \
    #     --standardized_chr_key $STANDARDIZED_CHR_KEY \
    #     --standardized_pos_key $STANDARDIZED_POS_KEY \
    #     --standardized_rsid_key $STANDARDIZED_RSID_KEY \
    #     --standardized_variant_id_key $STANDARDIZED_VAR_ID_KEY \
    #     --standardized_gene_id_key $STANDARDIZED_GENE_ID_KEY \
    #     --standardized_non_effect_allele_key $STANDARDIZED_NON_EFFECT_KEY \
    #     --standardized_effect_allele_key $STANDARDIZED_EFFECT_KEY \
    #     --standardized_p_key $STANDARDIZED_P_KEY \
    #     --standardized_beta_key $STANDARDIZED_BETA_KEY \
    #     --standardized_se_key $STANDARDIZED_SE_KEY \
    #     --standardized_statistic_key $STANDARDIZED_STATISTIC_KEY \
    #     --standardized_maf_key $STANDARDIZED_MAF_KEY \
    #     --standardized_var_beta_key $STANDARDIZED_VAR_BETA_KEY \
    #     --standardized_strata_key $STANDARDIZED_STRATA_KEY \
    #     --loci_left_bound_key $LOCI_LEFT_BOUND_KEY \
    #     --loci_right_bound_key $LOCI_RIGHT_BOUND_KEY \
    #     2>&1 | tee $LOGS_OUTPUT_DIR/sum_stats_transform.log

    # Standardized QTL Summary Statistics File
    filename=$(basename "$QTL_SUM_STATS_FILE" .tsv)
    STANDARDIZED_QTL_SUM_STATS_FP="${OUTPUT_DIR}/$filename.standardized.tsv"

    # Run check_coverage.py
    # python3 -u src/check_coverage.py \
    #     --out_ld_dir $LD_OUTPUT_DIR \
    #     --ld_manifest $LD_MANIFEST \
    #     --standardized_locus_id_key $MANIFEST_LOC_KEY \
    #     --standardized_chr_key $STANDARDIZED_CHR_KEY \
    #     --standardized_left_bound_key $LOCI_LEFT_BOUND_KEY \
    #     --standardized_right_bound_key $LOCI_RIGHT_BOUND_KEY \
    #     --standardized_snp_key $STANDARDIZED_VAR_ID_KEY \
    #     --standardized_pos_key $STANDARDIZED_POS_KEY \
    #     --standardized_var_key $STANDARDIZED_VAR_ID_KEY \
    #     --standardized_p_key $STANDARDIZED_P_KEY \
    #     --standardized_effect_allele_key $STANDARDIZED_EFFECT_KEY \
    #     --standardized_non_effect_allele_key $STANDARDIZED_NON_EFFECT_KEY \
    #     --manifest_bim_key $MANIFEST_BIM_KEY \
    #     --gwas_fp $STANDARDIZED_GWAS_SUM_STATS_FP \
    #     --qtl_fp $STANDARDIZED_QTL_SUM_STATS_FP \
    #     --qtl_strata_key $STANDARDIZED_STRATA_KEY \
    #     --high_overlap_min $HIGH_OVERLAP_MIN \
    #     --medium_overlap_min $MEDIUM_OVERLAP_MIN \
    #     --out_qc_dir $QC_OUTPUT_DIR \
    #     2>&1 | tee $LOGS_OUTPUT_DIR/check_coverage.log

    # Activate conda environment for coloc
    conda activate R

    # Run analyze_colocalization.R
    Rscript src/analyze_colocalization.R \
        --gwas_fp $STANDARDIZED_GWAS_SUM_STATS_FP \
        --qtl_fp $STANDARDIZED_QTL_SUM_STATS_FP \
        --gwas_sample_size $GWAS_SAMPLE_SIZE \
        --gwas_case_fraction $GWAS_CASE_FRACTION \
        --qtl_sample_size $QTL_SS_N \
        --output_dir $OUTPUT_DIR \
        --ld_dir $LD_OUTPUT_DIR \
        --qc_dir $QC_OUTPUT_DIR \
        --min_overlap $MIN_OVERLAP \
        --susie_min_snps $MIN_SNPS_SUSIE \
        --qtl_strata_key $STANDARDIZED_STRATA_KEY \
        2>&1 | tee $LOGS_OUTPUT_DIR/analyze_coloc.log

    # Log Status
    echo "Finished processing tissue: $TISSUE_NAME"

done < <(tail -n +2 $GTEX_MANIFEST_FILE)