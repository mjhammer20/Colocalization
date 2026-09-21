#!/bin/bash

# Configuration 
LOCI_FILE="/mnt/output/output/coloc/meta6_gwas_loci.merged.tsv"
LOCUS_ID_KEY="LOCUS_ID"
LD_REF_BFILE="/mnt/working/ref_panels/1kg_full/1kg_hg38_filtered"
LD_REF_ANCESTRY_KEEP_FILE="/mnt/working/ref_panels/1kg_full/ancestry_keep.txt"
LD_PANEL="1KG_v3_EUR"
LD_OUTPUT_DIR="/mnt/output/output/coloc/ld_1kg_full_EUR"
LD_MANIFEST="ld_manifest.tsv"
CHR_KEY="CHR"
LEFT_BOUND_KEY="LEFT_500KB"
RIGHT_BOUND_KEY="RIGHT_500KB"
MAF_MIN=0.01
GENO_MAX_MISSING=0.05

# Run compute_LD.py
python3 -u src/compute_LD.py \
    --loci_file "$LOCI_FILE" \
    --standardized_locus_id_key "$LOCUS_ID_KEY" \
    --ref_bfile "$LD_REF_BFILE" \
    --ancestry_keep_file "$LD_REF_ANCESTRY_KEEP_FILE" \
    --ld_panel "$LD_PANEL" \
    --ld_output_dir "$LD_OUTPUT_DIR" \
    --ld_manifest "$LD_MANIFEST" \
    --standardized_chr_key "$CHR_KEY" \
    --standardized_left_bound_key "$LEFT_BOUND_KEY" \
    --standardized_right_bound_key "$RIGHT_BOUND_KEY" \
    --maf_min "$MAF_MIN" \
    --geno_max "$GENO_MAX_MISSING" \
    2>&1 | tee /mnt/output/output/coloc/logs/compute_ld_1kg_full_EUR.log
