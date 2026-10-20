#!/bin/bash

# GTex Tissue Manifest File Parameters
GTEX_MANIFEST_FILE="/mnt/working/locus_reports/qtl/apaQTL/gtex/manifests/GTEx_tissue_manifest.tsv"

# Set the parent output directory for GTEx results
PARENT_OUTPUT_DIR="/mnt/output/output/coloc/gtex_apaQTL_sub_1kg_v3_hg38_EUR"

# GTex full tissue merged results file path
GTEX_FULL_TISSUE_COLOC_SUMMARY_FP="${PARENT_OUTPUT_DIR}/meta6_pd_full_gtex_v11_apaQTLS_coloc_summary_merged_updated_susie_2.tsv"

# Delete the existing full tissue merged results file if it exists
if [ -f "$GTEX_FULL_TISSUE_COLOC_SUMMARY_FP" ]; then
    echo "Deleting existing full tissue merged results file: $GTEX_FULL_TISSUE_COLOC_SUMMARY_FP"
    rm "$GTEX_FULL_TISSUE_COLOC_SUMMARY_FP"
fi

# Iterate through each line of the GTEx tissue manifest file
while IFS=$'\t' read -r tissue_name sample_size tissue_dir; do
    
    # Extract tissue information
    TISSUE_NAME="${tissue_name}"

    # Set output directory for the current tissue
    TISSUE_OUTPUT_DIR="${PARENT_OUTPUT_DIR}/${TISSUE_NAME}/updated_susie"

    # Set file paths for the current tissue
    COLOC_SUMMARY_FP="${TISSUE_OUTPUT_DIR}/qc/meta6_pd_full_gtex_v11_apaqtls_coloc_summary.tsv"

    # Initialize the full tissue merged results file if it doesn't exist
    if [ ! -f "$GTEX_FULL_TISSUE_COLOC_SUMMARY_FP" ]; then
        echo "Initializing full tissue merged results file with results from $TISSUE_NAME: $GTEX_FULL_TISSUE_COLOC_SUMMARY_FP"
        cp "$COLOC_SUMMARY_FP" "$GTEX_FULL_TISSUE_COLOC_SUMMARY_FP"
    else
        # Append the current tissue's merged results to the full tissue merged results file
        echo "Appending results from $TISSUE_NAME to full tissue merged results file: $GTEX_FULL_TISSUE_COLOC_SUMMARY_FP"
        awk 'NR>1 {print}' "$COLOC_SUMMARY_FP" >> "$GTEX_FULL_TISSUE_COLOC_SUMMARY_FP"
    fi

done < <(tail -n +2 $GTEX_MANIFEST_FILE)

echo "All tissues processed. Full tissue merged results saved to $GTEX_FULL_TISSUE_COLOC_SUMMARY_FP"
