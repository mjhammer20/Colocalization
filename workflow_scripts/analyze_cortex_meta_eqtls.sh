#!/bin/bash

# Entrez Email for NCBI API
ENTREZ_EMAIL="matt@datatecnica.com"

# Input Loci File Parameters
LOCI_FILE="/mnt/disks/output/output/coloc/meta6_gwas_loci.merged.tsv"
LOCI_LEFT_BOUND_KEY="LEFT_500KB"
LOCI_RIGHT_BOUND_KEY="RIGHT_500KB"

# Input Summary Statistics GWAS
GWAS_SUM_STATS_FILE="/mnt/disks/working/summary_analysis_data/gwas_stats/GP2_et_al_2025_PD_case_control_EUR_ALL_hg38_rsID.txt.gz"
GWAS_SUM_STATS_GENOME_BUILD="hg38"
HEADER_LINES=0
GWAS_SS_CHR_KEY="chromosome"
GWAS_SS_POS_KEY="base_pair_position"
GWAS_SS_RSID_KEY="rsID"
GWAS_SS_GENE_ID_KEY=""
GWAS_SS_NONEFFECT_KEY="other_allele"
GWAS_SS_EFFECT_KEY="effect_allele"
GWAS_SS_P_KEY="p_value"
GWAS_SS_BETA_KEY="beta"
GWAS_SS_SE_KEY="standard_error"
GWAS_SS_STATISTIC_KEY=""
GWAS_SS_MAF_KEY="effect_allele_frequency"
GWAS_SS_MAC_KEY=""
GWAS_SS_VAR_BETA_KEY=""
GWAS_SS_SDY_KEY=""
GWAS_SS_TISSUE_KEY=""
GWAS_SAMPLE_SIZE=226196
GWAS_CASE_FRACTION=0.207306937

# Input Summary Statistics QTL
QTL_SUM_STATS_FILE="/mnt/disks/working/locus_reports/qtl/eqtl/amp_ad/raw/Cortex_MetaAnalysis_ROSMAP_CMC_HBCC_Mayo_cis_eQTL_release.csv"
QTL_SUM_STATS_GENOME_BUILD="hg19"
HEADER_LINES=0
QTL_SS_CHR_KEY="chromosome"
QTL_SS_POS_KEY="snpLocation"
QTL_SS_RSID_KEY="snpid"
QTL_SS_GENE_SYMBOL_KEY="geneSymbol"
QTL_SS_GENE_ID_KEY="gene"
QTL_SS_NONEFFECT_KEY="A1"
QTL_SS_EFFECT_KEY="A2"
QTL_SS_P_KEY="pvalue"
QTL_SS_BETA_KEY="beta"
QTL_SS_SE_KEY=""
QTL_SS_N=1694
QTL_SS_STATISTIC_KEY="statistic"
QTL_SS_MAF_KEY="A2freq"
QTL_SS_MAC_KEY=""
QTL_SS_VAR_BETA_KEY=""
QTL_SS_SDY_KEY=""
QTL_SS_STRATA_KEY="region"
QTL_SAMPLE_SIZE=1694

# Standardized Summary Statistics File Paths
STANDARDIZED_GWAS_SUM_STATS_FP="/mnt/disks/output/output/coloc/not_validated/GP2_et_al_2025_PD_case_control_EUR_ALL_hg38_rsID.standardized.tsv"
STANDARDIZED_QTL_SUM_STATS_FP="/mnt/disks/output/output/coloc/amp_ad/cortex_meta/updated_retry/Cortex_MetaAnalysis_ROSMAP_CMC_HBCC_Mayo_cis_eQTL_release.standardized.tsv"

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
STANDARDIZED_STRATA_KEY="REGION"

# LD Manifest File Parameters
LD_OUTPUT_DIR="/mnt/disks/output/output/coloc/ld_1kg_full_EUR"
LD_MANIFEST="ld_manifest.tsv"
MANIFEST_LOC_KEY="LOCUS_ID"
MANIFEST_BIM_KEY="BIM"
HIGH_OVERLAP_MIN=0.90
MEDIUM_OVERLAP_MIN=0.70

# Susie Parameters
MIN_OVERLAP=10
MIN_SNPS_SUSIE=10

# Directories
GWAS_OUTPUT_DIR="/mnt/disks/output/output/coloc/not_validated"
GWAS_QC_OUTPUT_DIR="/mnt/disks/output/output/coloc/not_validated/qc"
OUTPUT_DIR="/mnt/disks/output/output/coloc/amp_ad/cortex_meta/updated_retry"
QC_OUTPUT_DIR="/mnt/disks/output/output/coloc/amp_ad/cortex_meta/updated_retry/qc"
LOGS_OUTPUT_DIR="/mnt/disks/output/output/coloc/amp_ad/cortex_meta/updated_retry/logs"
SUSIE_RESULTS_DIR="/mnt/disks/output/output/coloc/susie_results"
SUSIE_QC_DIR="/mnt/disks/output/output/coloc/susie_qc"

# Study Labels
GWAS_LABEL="META6_PD_Corrected"
QTL_LABEL="Cortex_Meta"

# Create necessary directories if they don't exist
mkdir -p "$OUTPUT_DIR"
mkdir -p "$QC_OUTPUT_DIR"
mkdir -p "$LOGS_OUTPUT_DIR"
mkdir -p "$SUSIE_RESULTS_DIR"
mkdir -p "$SUSIE_QC_DIR"
echo "Output Directory: $OUTPUT_DIR"
echo "QC Output Directory: $QC_OUTPUT_DIR"
echo "Logs Output Directory: $LOGS_OUTPUT_DIR"
echo "Susie Results Directory: $SUSIE_RESULTS_DIR"
echo "Susie QC Directory: $SUSIE_QC_DIR"

# Run standardize_sum_stats.py
# echo "Running standardize_sum_stats.py for GWAS..."
# python3 -u src/standardize_sum_stats.py \
#     --sum_stats_file $GWAS_SUM_STATS_FILE \
#     --sum_stats_genome_build $GWAS_SUM_STATS_GENOME_BUILD \
#     --loci_file $LOCI_FILE \
#     --output_dir $GWAS_OUTPUT_DIR \
#     --qc_output_dir $GWAS_QC_OUTPUT_DIR \
#     --header_lines $HEADER_LINES \
#     --entrez_email $ENTREZ_EMAIL \
#     --ss_chr_key $GWAS_SS_CHR_KEY \
#     --ss_pos_key $GWAS_SS_POS_KEY \
#     --ss_rsid_key $GWAS_SS_RSID_KEY \
#     --ss_non_effect_allele_key $GWAS_SS_NONEFFECT_KEY \
#     --ss_effect_allele_key $GWAS_SS_EFFECT_KEY \
#     --ss_p_key $GWAS_SS_P_KEY \
#     --ss_se_key $GWAS_SS_SE_KEY \
#     --ss_beta_key $GWAS_SS_BETA_KEY \
#     --ss_maf_key $GWAS_SS_MAF_KEY \
#     --ss_n $GWAS_SAMPLE_SIZE \
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
#     --loci_left_bound_key $LOCI_LEFT_BOUND_KEY \
#     --loci_right_bound_key $LOCI_RIGHT_BOUND_KEY \
#     > "$LOGS_OUTPUT_DIR/sum_stats_standardization_gwas.log" 2>&1

# Run standardize_sum_stats.py for QTL
# echo "Running standardize_sum_stats.py for QTL..."
# python3 -u src/standardize_sum_stats.py \
#     --sum_stats_file $QTL_SUM_STATS_FILE \
#     --sum_stats_genome_build $QTL_SUM_STATS_GENOME_BUILD \
#     --loci_file $LOCI_FILE \
#     --output_dir $OUTPUT_DIR \
#     --qc_output_dir $QC_OUTPUT_DIR \
#     --header_lines $HEADER_LINES \
#     --entrez_email $ENTREZ_EMAIL \
#     --ss_chr_key $QTL_SS_CHR_KEY \
#     --ss_pos_key $QTL_SS_POS_KEY \
#     --ss_rsid_key $QTL_SS_RSID_KEY \
#     --ss_gene_id_key $QTL_SS_GENE_ID_KEY \
#     --ss_gene_symbol_key $QTL_SS_GENE_SYMBOL_KEY \
#     --ss_non_effect_allele_key $QTL_SS_NONEFFECT_KEY \
#     --ss_effect_allele_key $QTL_SS_EFFECT_KEY \
#     --ss_statistic_key $QTL_SS_STATISTIC_KEY \
#     --ss_p_key $QTL_SS_P_KEY \
#     --ss_beta_key $QTL_SS_BETA_KEY \
#     --ss_maf_key $QTL_SS_MAF_KEY \
#     --ss_n $QTL_SAMPLE_SIZE \
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
#     > "$LOGS_OUTPUT_DIR/sum_stats_standardization_qtl.log" 2>&1

# Run the check_coverage.py script with the specified parameters
echo "Running check_coverage.py..."
python3 -u src/check_coverage.py \
    --out_ld_dir $LD_OUTPUT_DIR \
    --ld_manifest $LD_MANIFEST \
    --standardized_locus_id_key $MANIFEST_LOC_KEY \
    --standardized_chr_key $STANDARDIZED_CHR_KEY \
    --standardized_left_bound_key $LOCI_LEFT_BOUND_KEY \
    --standardized_right_bound_key $LOCI_RIGHT_BOUND_KEY \
    --standardized_snp_key $STANDARDIZED_VAR_ID_KEY \
    --standardized_pos_key $STANDARDIZED_POS_KEY \
    --standardized_var_key $STANDARDIZED_VAR_ID_KEY \
    --standardized_p_key $STANDARDIZED_P_KEY \
    --standardized_effect_allele_key $STANDARDIZED_EFFECT_KEY \
    --standardized_non_effect_allele_key $STANDARDIZED_NON_EFFECT_KEY \
    --manifest_bim_key $MANIFEST_BIM_KEY \
    --gwas_fp $STANDARDIZED_GWAS_SUM_STATS_FP \
    --qtl_fp $STANDARDIZED_QTL_SUM_STATS_FP \
    --qtl_strata_key $STANDARDIZED_STRATA_KEY \
    --high_overlap_min $HIGH_OVERLAP_MIN \
    --medium_overlap_min $MEDIUM_OVERLAP_MIN \
    --out_qc_dir $QC_OUTPUT_DIR \
    --gwas_label $GWAS_LABEL \
    --qtl_label $QTL_LABEL \
    > "$LOGS_OUTPUT_DIR/check_coverage.log" 2>&1

# echo "check_coverage.py completed. Log available at $LOGS_OUTPUT_DIR/check_coverage.log"

# Run analyze_colocalization.R
echo "Running analyze_colocalization.R..."
Rscript src/analyze_colocalization.R \
    --gwas_fp $STANDARDIZED_GWAS_SUM_STATS_FP \
    --qtl_fp $STANDARDIZED_QTL_SUM_STATS_FP \
    --gwas_sample_size $GWAS_SAMPLE_SIZE \
    --gwas_case_fraction $GWAS_CASE_FRACTION \
    --qtl_sample_size $QTL_SAMPLE_SIZE \
    --output_dir $OUTPUT_DIR \
    --ld_dir $LD_OUTPUT_DIR \
    --qc_dir $QC_OUTPUT_DIR \
    --min_overlap $MIN_OVERLAP \
    --susie_min_snps $MIN_SNPS_SUSIE \
    --qtl_strata_key $STANDARDIZED_STRATA_KEY \
    --susie_results_dir $SUSIE_RESULTS_DIR \
    --susie_qc_dir $SUSIE_QC_DIR \
    --gwas_label $GWAS_LABEL \
    --qtl_label $QTL_LABEL \
    2>&1 | tee "$LOGS_OUTPUT_DIR/analyze_coloc.log"

echo "analyze_colocalization.R completed. Log available at $LOGS_OUTPUT_DIR/analyze_coloc.log"
echo "Colocalization analysis completed. Results are in $OUTPUT_DIR"