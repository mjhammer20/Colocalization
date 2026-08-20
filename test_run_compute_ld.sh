#!/bin/bash

# Configuration 
LOCI_FILE="/mnt/output/output/coloc/meta6_gwas_loci.merged.tsv"
LOCUS_ID_KEY="LOCUS_ID"
LD_REF_BFILE="/mnt/working/ld_reference_data/1000G_EUR/1000G.EUR.QC.hg38"
LD_PANEL="1000G_EUR"
LD_OUTPUT_DIR="/mnt/output/output/coloc/ld"
LD_MANIFEST="ld_manifest.tsv"
CHR_KEY="CHR"
LEFT_BOUND_KEY="LEFT_500KB"
RIGHT_BOUND_KEY="RIGHT_500KB"

# Run compute_ld.py
python3 -u src/compute_ld.py \
    --loci_file "$LOCI_FILE" \
    --locus_id_key "$LOCUS_ID_KEY" \
    --ld_ref_bfile "$LD_REF_BFILE" \
    --ld_panel "$LD_PANEL" \
    --ld_output_dir "$LD_OUTPUT_DIR" \
    --ld_manifest "$LD_MANIFEST" \
    --chr_key "$CHR_KEY" \
    --left_bound_key "$LEFT_BOUND_KEY" \
    --right_bound_key "$RIGHT_BOUND_KEY" \
    2>&1 | tee /mnt/output/output/coloc/logs/compute_ld.log
