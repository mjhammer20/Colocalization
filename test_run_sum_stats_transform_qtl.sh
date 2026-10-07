# Entrez Email for NCBI API
ENTREZ_EMAIL="matt@datatecnica.com"

# Input Loci File Parameters
LOCI_FILE="/mnt/output/output/coloc/meta6_gwas_loci.merged.tsv"
LOCI_LEFT_BOUND_KEY="LEFT_500KB"
LOCI_RIGHT_BOUND_KEY="RIGHT_500KB"

# Liftover Chain File Path
LIFTOVER_CHAIN_FILE="/mnt/working/chain_files/hg19ToHg38.chain.gz"

# Input Summary Statistics QTL
QTL_SUM_STATS_FILE="/mnt/working/locus_reports/qtl/eqtl/amp_ad/raw/Cortex_MetaAnalysis_ROSMAP_CMC_HBCC_Mayo_cis_eQTL_release.csv"
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

# Output Directories
OUTPUT_DIR="/mnt/output/output/coloc/amp_ad/cortex_meta/updated_retry"
LOGS_OUTPUT_DIR="/mnt/output/output/coloc/amp_ad/cortex_meta/updated_retry/logs"
QC_OUTPUT_DIR="/mnt/output/output/coloc/amp_ad/cortex_meta/updated_retry/qc"

# Study Labels
GWAS_LABEL="META6_PD"
QTL_LABEL="Cortex_Meta"

echo "Running standardize_sum_stats.py for QTL..."
python3 -u src/standardize_sum_stats.py \
    --sum_stats_file $QTL_SUM_STATS_FILE \
    --sum_stats_genome_build $QTL_SUM_STATS_GENOME_BUILD \
    --liftover_chain_file $LIFTOVER_CHAIN_FILE \
    --loci_file $LOCI_FILE \
    --output_dir $OUTPUT_DIR \
    --qc_output_dir $QC_OUTPUT_DIR \
    --header_lines $HEADER_LINES \
    --entrez_email $ENTREZ_EMAIL \
    --ss_chr_key $QTL_SS_CHR_KEY \
    --ss_pos_key $QTL_SS_POS_KEY \
    --ss_rsid_key $QTL_SS_RSID_KEY \
    --ss_gene_id_key $QTL_SS_GENE_ID_KEY \
    --ss_gene_symbol_key $QTL_SS_GENE_SYMBOL_KEY \
    --ss_non_effect_allele_key $QTL_SS_NONEFFECT_KEY \
    --ss_effect_allele_key $QTL_SS_EFFECT_KEY \
    --ss_statistic_key $QTL_SS_STATISTIC_KEY \
    --ss_p_key $QTL_SS_P_KEY \
    --ss_beta_key $QTL_SS_BETA_KEY \
    --ss_maf_key $QTL_SS_MAF_KEY \
    --ss_n $QTL_SAMPLE_SIZE \
    --standardized_chr_key $STANDARDIZED_CHR_KEY \
    --standardized_pos_key $STANDARDIZED_POS_KEY \
    --standardized_rsid_key $STANDARDIZED_RSID_KEY \
    --standardized_variant_id_key $STANDARDIZED_VAR_ID_KEY \
    --standardized_gene_id_key $STANDARDIZED_GENE_ID_KEY \
    --standardized_non_effect_allele_key $STANDARDIZED_NON_EFFECT_KEY \
    --standardized_effect_allele_key $STANDARDIZED_EFFECT_KEY \
    --standardized_p_key $STANDARDIZED_P_KEY \
    --standardized_beta_key $STANDARDIZED_BETA_KEY \
    --standardized_se_key $STANDARDIZED_SE_KEY \
    --standardized_statistic_key $STANDARDIZED_STATISTIC_KEY \
    --standardized_maf_key $STANDARDIZED_MAF_KEY \
    --standardized_var_beta_key $STANDARDIZED_VAR_BETA_KEY \
    --standardized_strata_key $STANDARDIZED_STRATA_KEY \
    --loci_left_bound_key $LOCI_LEFT_BOUND_KEY \
    --loci_right_bound_key $LOCI_RIGHT_BOUND_KEY \
    > "$LOGS_OUTPUT_DIR/sum_stats_standardization_qtl.log" 2>&1
