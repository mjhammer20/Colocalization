#!/bin/bash

# GTex Tissue Manifest File Parameters
GTEX_MANIFEST_FILE="/mnt/working/locus_reports/qtl/sqtl/gtex/all_associations/GTEx_Analysis_v11_sQTL_all_associations_prepared/GTEx_tissue_manifest.tsv"

# Initialize tissue counter
TISSUE_COUNTER=0

# Set the parent output directory for GTEx results
PARENT_OUTPUT_DIR="/mnt/output/output/coloc/gtex_sqtl_sub_1kg_v3_hg38_EUR"

# GTex full tissue merged results file path
GTEX_FULL_TISSUE_MERGED_RESULTS_FP="${PARENT_OUTPUT_DIR}/meta6_pd_full_gtex_v11_sqtls_coloc_all_results_merged.tsv"

# Iterate through each line of the GTEx tissue manifest file
while IFS=$'\t' read -r tissue_name sample_size tissue_dir; do
    
    # Extract tissue information
    TISSUE_NAME="${tissue_name}"

    # Set output directory for the current tissue
    TISSUE_OUTPUT_DIR="${PARENT_OUTPUT_DIR}/${TISSUE_NAME}"

    # Set file paths for the current tissue
    SUSIE_RESULTS_FP="${TISSUE_OUTPUT_DIR}/meta6_pd_full_gtex_v11_sqtls_susie_coloc_all_results.tsv"
    COLOC_ABF_RESULTS_FP="${TISSUE_OUTPUT_DIR}/meta6_pd_full_gtex_v11_sqtls_coloc_abf_all_results.tsv"
    TISSUE_MERGED_RESULTS_FP="${TISSUE_OUTPUT_DIR}/meta6_pd_full_gtex_v11_sqtls_coloc_all_results_merged.tsv"

    # Merge the results for the current tissue
    echo "Merging results for tissue: $TISSUE_NAME"
    python3 -u src/merge_results.py \
        --susie_results_fp "$SUSIE_RESULTS_FP" \
        --coloc_abf_results_fp "$COLOC_ABF_RESULTS_FP" \
        --merged_results_fp "$TISSUE_MERGED_RESULTS_FP"

    # Initialize the full tissue merged results file if it doesn't exist
    if [ ! -f "$GTEX_FULL_TISSUE_MERGED_RESULTS_FP" ]; then
        echo "Initializing full tissue merged results file with results from $TISSUE_NAME: $GTEX_FULL_TISSUE_MERGED_RESULTS_FP"
        cp "$TISSUE_MERGED_RESULTS_FP" "$GTEX_FULL_TISSUE_MERGED_RESULTS_FP"
    else
        # Append the current tissue's merged results to the full tissue merged results file
        echo "Appending results from $TISSUE_NAME to full tissue merged results file: $GTEX_FULL_TISSUE_MERGED_RESULTS_FP"
        awk 'NR>1 {print}' "$TISSUE_MERGED_RESULTS_FP" >> "$GTEX_FULL_TISSUE_MERGED_RESULTS_FP"
    fi

done < <(tail -n +2 $GTEX_MANIFEST_FILE)

echo "All tissues processed. Full tissue merged results saved to $GTEX_FULL_TISSUE_MERGED_RESULTS_FP"
