#!/bin/bash

# Set file paths
OUTPUT_DIR="/mnt/output/output/coloc/amp_ad/cortex_meta/1kg_v3_hg38_EUR"
SUSIE_RESULTS_FP="$OUTPUT_DIR/meta6_pd_full_cortex_meta_coloc_abf_all_results.tsv"
COLOC_ABF_RESULTS_FP="$OUTPUT_DIR/meta6_pd_full_cortex_meta_susie_coloc_all_results.tsv"
MERGED_RESULTS_FP="$OUTPUT_DIR/meta6_pd_full_cortex_meta_coloc_all_results_merged_2.tsv"

# Run merge_results.py to merge the results
echo "Merging results from $SUSIE_RESULTS_FP and $COLOC_ABF_RESULTS_FP"
python3 -u src/merge_results.py \
    --susie_results_fp "$SUSIE_RESULTS_FP" \
    --coloc_abf_results_fp "$COLOC_ABF_RESULTS_FP" \
    --merged_results_fp "$MERGED_RESULTS_FP"
echo "Merge completed. Merged results saved to $MERGED_RESULTS_FP"

