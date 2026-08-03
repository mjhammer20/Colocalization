#!/usr/bin/env bash

INPUT_FILES=(
  "/mnt/working/locus_reports/qtl/eqtl/coloc_ready/need_retry/Cortex_MetaAnalysis_ROSMAP_CMC_HBCC_Mayo_cis_eQTL_release.tsv"
  "/mnt/working/locus_reports/qtl/eqtl/coloc_ready/need_retry/CER_Mayo_cis_eQTL_release_no_rsid.tsv"
  "/mnt/working/locus_reports/qtl/eqtl/coloc_ready/need_retry/DLPFC_ROSMAP_cis_eQTL_release_no_rsid.tsv"
  "/mnt/working/locus_reports/qtl/eqtl/coloc_ready/need_retry/TCX_Mayo_cis_eQTL_release_no_rsid.tsv"
)
EXPANDED_RANGES_FILE="/mnt/working/locus_reports/expanded_locus_ranges_by_chromosome.json"
OUTPUT_DIR="/mnt/working/locus_reports/qtl/eqtl/coloc_ready/need_retry"
HEADER_ROWS=0
CHR_COL="chromosome"
POS_COL="snpLocation"
RSID_COL="snpid"
GENE_ID_COL="gene"
EFFECT_ALLELE_COL="A2"
NON_EFFECT_ALLELE_COL="A1"
P_COL="pvalue"
SE_COL=""
BETA_COL="beta"
MAF_COL="A2freq"
N_COL=""
STATISTIC_COL="statistic"
MAC_COL=""
VAR_BETA_COL=""
SDY_COL=""
REGION_COL="Region"


for INPUT_FILE in "${INPUT_FILES[@]}"; do
    python3 -u /mnt/working/scripts/coloc_transform.py \
        --input_file "$INPUT_FILE" \
        --expanded_ranges_file "$EXPANDED_RANGES_FILE" \
        --output_dir "$OUTPUT_DIR" \
        --chr_col "$CHR_COL" \
        --pos_col "$POS_COL" \
        --rsid_col "$RSID_COL" \
        --gene_id_col "$GENE_ID_COL" \
        --effect_allele_col "$EFFECT_ALLELE_COL" \
        --non_effect_allele_col "$NON_EFFECT_ALLELE_COL" \
        --p_col "$P_COL" \
        --beta_col "$BETA_COL" \
        --statistic_col "$STATISTIC_COL" \
        --maf_col "$MAF_COL" 
done