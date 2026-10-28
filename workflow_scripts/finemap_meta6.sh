#!/usr/bin/env bash

# Configuration

## LD Manifest File Parameters
LD_OUTPUT_DIR="/mnt/output/output/coloc/ld_1kg_v3_hg38_EUR"
LD_MANIFEST="ld_manifest.tsv"
MANIFEST_LOC_KEY="LOCUS_ID"
MANIFEST_BIM_KEY="BIM"

## Summary Statistics Parameters
SS_PATH="/mnt/output/output/coloc/not_validated/GP2_et_al_2025_PD_case_control_EUR_ALL_hg38_rsID.standardized.tsv"
SS_SAMPLE_SIZE=1827641
SS_CASE_FRACTION=0.04445895
SS_TYPE="cc"
SS_LABEL="META6_PD"

## Output Directory for Finemap Results
SUSIE_RESULTS_DIR="/mnt/output/output/coloc/susie_results_1kg_v3_hg38_EUR_direct"

# Create necessary directories if they don't exist
mkdir -p "$SUSIE_RESULTS_DIR"

# Run the R script for finemapping with susie
echo "Running finemap_susie.R..."
Rscript src/finemap_susie.R \
    --ld_dir "$LD_OUTPUT_DIR" \
    --ld_manifest "$LD_MANIFEST" \
    --manifest_locus_id_key "$MANIFEST_LOC_KEY" \
    --manifest_bim_key "$MANIFEST_BIM_KEY" \
    --ss_fp "$SS_PATH" \
    --ss_sample_size "$SS_SAMPLE_SIZE" \
    --ss_case_fraction "$SS_CASE_FRACTION" \
    --ss_type "$SS_TYPE" \
    --ss_label "$SS_LABEL" \
    --susie_results_dir "$SUSIE_RESULTS_DIR" \
    2>&1 | tee "$SUSIE_RESULTS_DIR/finemap_susie.log"

echo "Finemapping with susie completed. Results saved in $SUSIE_RESULTS_DIR"