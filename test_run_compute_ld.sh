#!/bin/bash

# Configuration 
LOCI_FILE="/mnt/output/output/coloc/meta6_gwas_loci.merged.tsv"
LOCUS_ID_KEY="LOCUS_ID"
LD_REF_BFILE="/mnt/working/ref_panels/1kg/ref_panel/ref_panel_gp2_prune_rm_underperform_pos_update"
LD_PANEL="1KG_30x_HGDP_ASHK"
LD_OUTPUT_DIR="/mnt/output/output/coloc/ld"
LD_MANIFEST="ld_manifest.tsv"
CHR_KEY="CHR"
LEFT_BOUND_KEY="LEFT_500KB"
RIGHT_BOUND_KEY="RIGHT_500KB"
MAF_MIN=0.005
GENO_MAX_MISSING=0.1

# Run compute_LD.py
python3 -u src/compute_LD.py \
    --loci_file "$LOCI_FILE" \
    --standardized_locus_id_key "$LOCUS_ID_KEY" \
    --ref_bfile "$LD_REF_BFILE" \
    --ld_panel "$LD_PANEL" \
    --ld_output_dir "$LD_OUTPUT_DIR" \
    --ld_manifest "$LD_MANIFEST" \
    --standardized_chr_key "$CHR_KEY" \
    --standardized_left_bound_key "$LEFT_BOUND_KEY" \
    --standardized_right_bound_key "$RIGHT_BOUND_KEY" \
    --maf_min "$MAF_MIN" \
    --geno_max "$GENO_MAX_MISSING" \
    2>&1 | tee /mnt/output/output/coloc/logs/compute_ld.log
