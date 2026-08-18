#!/bin/bash

# Configuration 

GWAS_SS="/mnt/working/summary_analysis_data/gwas_stats/GP2_et_al_2025_PD_case_control_EUR_ALL_hg38_rsID.txt.gz"
OUTPUT_PREFIX="/mnt/output/output/coloc/meta6_gwas_loci"
SNP_KEY="variantID_hg38"
CHR_KEY="chromosome"
POS_KEY="base_pair_position"

# Run define_loci.py
python3 -u src/define_loci.py --input_file $GWAS_SS --output_file_prefix $OUTPUT_PREFIX --snp_id_col $SNP_KEY --chr_col $CHR_KEY --pos_col $POS_KEY
