#!/bin/bash

# ===== COMMAND LINE ARGUMENT PARSING =====
usage() {
    cat <<EOF
Usage: $0 [OPTIONS]

OPTIONS:
    -s, --start START_INDEX      Start index for batch processing (0-based, default: 0)
    -e, --end END_INDEX          End index for batch processing (exclusive, default: all rows)
    -h, --help                   Show this help message

EXAMPLES:
    # Process all tissues (default)
    bash analyze_gtex_sqtls_batched_local.sh

    # Process tissues 0-4 (first 5 tissues)
    bash analyze_gtex_sqtls_batched_local.sh -s 0 -e 5

    # Process tissues 5-9 (second batch of 5)
    bash analyze_gtex_sqtls_batched_local.sh --start 5 --end 10

    # Process tissues 10 onwards (all remaining)
    bash analyze_gtex_sqtls_batched_local.sh -s 10
EOF
    exit 1
}

# Default values
BATCH_START_INDEX=0
BATCH_END_INDEX=""

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -s|--start)
            BATCH_START_INDEX="$2"
            shift 2
            ;;
        -e|--end)
            BATCH_END_INDEX="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "ERROR: Unknown option: $1"
            usage
            ;;
    esac
done

# Validate batch indices
if ! [[ "$BATCH_START_INDEX" =~ ^[0-9]+$ ]]; then
    echo "ERROR: START_INDEX must be a non-negative integer, got: $BATCH_START_INDEX"
    exit 1
fi

if [[ -n "$BATCH_END_INDEX" ]] && ! [[ "$BATCH_END_INDEX" =~ ^[0-9]+$ ]]; then
    echo "ERROR: END_INDEX must be a non-negative integer, got: $BATCH_END_INDEX"
    exit 1
fi

echo "Batch parameters: START=$BATCH_START_INDEX, END=${BATCH_END_INDEX:-all}"

# ======== Configuration Parameters ==========

# Standardization Parameters
ENTREZ_EMAIL="matt@datatecnica.com"
QTL_SUM_STATS_GENOME_BUILD="hg38"
HEADER_LINES=0
SS_CHR_KEY="chr"
SS_POS_KEY="pos"
SS_RSID_KEY="rsid"
SS_GENE_ID_KEY="gene_id"
SS_GENE_SYMBOL_KEY="gene_symbol"
SS_PHENOTYPE_ID_KEY="phenotype_id"
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
SS_STRATA_KEY=""

# Input Loci File Parameters
LOCI_FILE="/mnt/disks/output/output/coloc/meta6_gwas_loci.merged.tsv"
LOCI_LEFT_BOUND_KEY="LEFT_500KB"
LOCI_RIGHT_BOUND_KEY="RIGHT_500KB"

# GWAS Summary Statistics File Parameters
STANDARDIZED_GWAS_SUM_STATS_FP="/mnt/disks/output/output/coloc/not_validated/GP2_et_al_2025_PD_case_control_EUR_ALL_hg38_rsID.standardized.tsv"
GWAS_SAMPLE_SIZE=226196
GWAS_CASE_FRACTION=0.207306937

# Standardized Summary Statistics Parameters
STANDARDIZED_CHR_KEY="CHR"
STANDARDIZED_POS_KEY="BP"
STANDARDIZED_RSID_KEY="SNP"
STANDARDIZED_VAR_ID_KEY="VAR"
STANDARDIZED_GENE_ID_KEY="GENE"
STANDARDIZED_PHENOTYPE_ID_KEY="PHENOTYPE"
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
LD_OUTPUT_DIR="/mnt/disks/output/output/coloc/ld_1kg_full_EUR"
LD_MANIFEST="ld_manifest.tsv"
MANIFEST_LOC_KEY="LOCUS_ID"
MANIFEST_BIM_KEY="BIM"
HIGH_OVERLAP_MIN=0.90
MEDIUM_OVERLAP_MIN=0.70

# Susie Parameters
MIN_OVERLAP=10
MIN_SNPS_SUSIE=10

# GTex Tissue Manifest File Parameters
GTEX_MANIFEST_FILE="/mnt/disks/working/locus_reports/qtl/sqtl/gtex/all_associations/GTEx_Analysis_v11_sQTL_all_associations_prepared/GTEx_tissue_manifest.tsv"

# SuSiE Results Directory (optional, defaults to output directory if not provided)
SUSIE_RESULTS_DIR="/mnt/disks/output/output/coloc/susie_results"

# Study Labels
GWAS_LABEL="META6_PD_Corrected"
QTL_LABEL="GTEx_v11_sQTLs"

# ===== Main Workflow Loop =====

# Count total tissues in manifest
TOTAL_TISSUES=$(tail -n +2 "$GTEX_MANIFEST_FILE" | wc -l)
echo "Total tissues in manifest: $TOTAL_TISSUES"

# Validate and set batch end index
if [[ -z "$BATCH_END_INDEX" ]]; then
    BATCH_END_INDEX=$TOTAL_TISSUES
fi

# Validate batch indices
if [[ $BATCH_START_INDEX -ge $TOTAL_TISSUES ]]; then
    echo "ERROR: BATCH_START_INDEX ($BATCH_START_INDEX) >= total tissues ($TOTAL_TISSUES)"
    exit 1
fi

if [[ $BATCH_END_INDEX -gt $TOTAL_TISSUES ]]; then
    echo "WARNING: BATCH_END_INDEX ($BATCH_END_INDEX) > total tissues ($TOTAL_TISSUES), capping to $TOTAL_TISSUES"
    BATCH_END_INDEX=$TOTAL_TISSUES
fi

if [[ $BATCH_START_INDEX -gt $BATCH_END_INDEX ]]; then
    echo "ERROR: BATCH_START_INDEX ($BATCH_START_INDEX) > BATCH_END_INDEX ($BATCH_END_INDEX)"
    exit 1
fi

echo "Processing tissues $BATCH_START_INDEX to $BATCH_END_INDEX (batch range: $BATCH_START_INDEX-$BATCH_END_INDEX)"
echo "Total tissues in batch: $((BATCH_END_INDEX - BATCH_START_INDEX))"

# Loop through each tissue in the GTEx manifest and run workflow
TISSUE_INDEX=0
while IFS=$'\t' read -r tissue_name sample_size tissue_dir; do
    
    # Skip tissues outside the batch range
    if [[ $TISSUE_INDEX -lt $BATCH_START_INDEX ]]; then
        ((TISSUE_INDEX++))
        continue
    fi
    
    if [[ $TISSUE_INDEX -gt $BATCH_END_INDEX ]]; then
        break
    fi
    
    # Extract tissue information
    TISSUE_NAME="${tissue_name}"
    BATCH_POSITION=$((TISSUE_INDEX - BATCH_START_INDEX + 1))
    TISSUE_NUMBER=$((TISSUE_INDEX + 1))
    echo ""
    echo "[$BATCH_POSITION/$((BATCH_END_INDEX - BATCH_START_INDEX))] Processing tissue #$TISSUE_NUMBER: $TISSUE_NAME"
    echo "Tissue: $TISSUE_NAME"
    QTL_SS_N="${sample_size}"
    echo "Sample Size: $QTL_SS_N"
    TISSUE_DIR="${tissue_dir}"
    # Replace /mnt/working with /mnt/disks/working (correct path for worker VMs)
    TISSUE_DIR="${TISSUE_DIR//\/mnt\/working/\/mnt\/disks\/working}"
    echo "Pre-processed Tissue Summary Statistics Directory: $TISSUE_DIR"

    # Output Directories
    OUTPUT_DIR="/mnt/disks/output/output/coloc/gtex_sqtl/${TISSUE_NAME}"
    TEMP_OUTPUT_DIR="${OUTPUT_DIR}/temp"
    QC_OUTPUT_DIR="/mnt/disks/output/output/coloc/gtex_sqtl/${TISSUE_NAME}/qc"
    LOGS_OUTPUT_DIR="/mnt/disks/output/output/coloc/gtex_sqtl/${TISSUE_NAME}/logs"
    echo "Output Directory: $OUTPUT_DIR"
    echo "Temporary Output Directory: $TEMP_OUTPUT_DIR"
    echo "QC Output Directory: $QC_OUTPUT_DIR"
    echo "Logs Output Directory: $LOGS_OUTPUT_DIR"

    # Create output directories if they don't exist
    mkdir -p "${OUTPUT_DIR}"
    mkdir -p "${TEMP_OUTPUT_DIR}"
    mkdir -p "${QC_OUTPUT_DIR}"
    mkdir -p "${LOGS_OUTPUT_DIR}"
    mkdir -p "${SUSIE_RESULTS_DIR}"

    # Standardized QTL Summary Statistics File
    STANDARDIZED_QTL_SUM_STATS_FP="${OUTPUT_DIR}/GTEx_v11_cis_sQTLs_${TISSUE_NAME}.standardized.tsv"

    # Check if the standardized file already exists
    if [[ -f "$STANDARDIZED_QTL_SUM_STATS_FP" ]]; then
        echo "Standardized QTL summary statistics file already exists."
        echo "Skipping standardize_sum_stats.py step..."
    else
        echo "Standardized QTL summary statistics file not found. Running standardize_sum_stats.py..."

        # Get all parquet files
        files=("$TISSUE_DIR"/*.parquet)
        TOTAL_FILES=${#files[@]}
        echo "Total parquet files found: $TOTAL_FILES"

        # Iterate over all files for this tissue
        for ((file_idx = 0; file_idx < TOTAL_FILES; file_idx++)); do
            f="${files[$file_idx]}"
            if [[ -f "$f" ]]; then
                QTL_SUM_STATS_FILE="$f"
                echo "  [$((file_idx + 1))/$TOTAL_FILES] Processing: $(basename "$QTL_SUM_STATS_FILE")"

                # Run standardize_sum_stats.py
                python3 -u src/standardize_sum_stats.py \
                    --sum_stats_file "$QTL_SUM_STATS_FILE" \
                    --sum_stats_genome_build $QTL_SUM_STATS_GENOME_BUILD \
                    --loci_file $LOCI_FILE \
                    --output_dir $TEMP_OUTPUT_DIR \
                    --qc_output_dir $QC_OUTPUT_DIR \
                    --entrez_email $ENTREZ_EMAIL \
                    --ss_rsid_key $SS_RSID_KEY \
                    --ss_chr_key $SS_CHR_KEY \
                    --ss_pos_key $SS_POS_KEY \
                    --ss_gene_symbol_key $SS_GENE_SYMBOL_KEY \
                    --ss_gene_id_key $SS_GENE_ID_KEY \
                    --ss_phenotype_id_key $SS_PHENOTYPE_ID_KEY \
                    --ss_non_effect_allele_key $SS_NONEFFECT_KEY \
                    --ss_effect_allele_key $SS_EFFECT_KEY \
                    --ss_p_key $SS_P_KEY \
                    --ss_beta_key $SS_BETA_KEY \
                    --ss_se_key $SS_SE_KEY \
                    --ss_maf_key $SS_MAF_KEY \
                    --ss_n $QTL_SS_N \
                    --ss_strata_value $TISSUE_NAME \
                    --standardized_chr_key $STANDARDIZED_CHR_KEY \
                    --standardized_pos_key $STANDARDIZED_POS_KEY \
                    --standardized_rsid_key $STANDARDIZED_RSID_KEY \
                    --standardized_variant_id_key $STANDARDIZED_VAR_ID_KEY \
                    --standardized_gene_id_key $STANDARDIZED_GENE_ID_KEY \
                    --standardized_phenotype_id_key $STANDARDIZED_PHENOTYPE_ID_KEY \
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
                    >> "$LOGS_OUTPUT_DIR/standardize_sum_stats.log" 2>&1
            fi
        done

        # Concatenate all standardized QTL summary statistics files into one file per tissue
        awk 'NR==1 || FNR>1' "$TEMP_OUTPUT_DIR"/*standardized.tsv > "$STANDARDIZED_QTL_SUM_STATS_FP"
        
        # Verify concatenation was successful
        if [[ ! -s "$STANDARDIZED_QTL_SUM_STATS_FP" ]]; then
            echo "ERROR: Concatenated file is empty or missing: $STANDARDIZED_QTL_SUM_STATS_FP"
            ((TISSUE_INDEX++))
            continue
        fi
        
        # Clean up temp directory
        rm -r "$TEMP_OUTPUT_DIR"
        echo "Successfully created standardized QTL summary statistics file."
    fi


    CHECK_COVERAGE_RESULTS_FP="${QC_OUTPUT_DIR}/${GWAS_LABEL}_${QTL_LABEL}_coverage_by_locus.tsv"

    if [[ -f "$CHECK_COVERAGE_RESULTS_FP" ]]; then
        echo "Check coverage results file already exists."
        echo "Skipping check_coverage.py step..."
    else
        echo "Check coverage results file not found. Running check_coverage.py..."
         # Run check_coverage.py
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
    fi
   
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
        --susie_results_dir $SUSIE_RESULTS_DIR \
        --gwas_label $GWAS_LABEL \
        --qtl_label $QTL_LABEL \
        --standardized_gene_id_key $STANDARDIZED_PHENOTYPE_ID_KEY \
        > "$LOGS_OUTPUT_DIR/analyze_coloc.log" 2>&1

    # Log Status
    echo "Finished processing tissue: $TISSUE_NAME"
    
    # Increment tissue index
    ((TISSUE_INDEX++))

done < <(tail -n +2 $GTEX_MANIFEST_FILE)

echo ""
echo "Batch processing complete!"