# Imports
import os 
import pandas as pd
from pathlib import Path
import coloc_transform
import gc

# Define parameters
regions = ['Astrocytes', 'Bulk', 'Endothelial_cells', 'Excitatory_Neurons', 'Inhibitory_Neurons', 'Microglia', 'Oligodendrocytes', 'OPCs', 'Pericytes']
input_files_dir = Path('/mnt/working/locus_reports/qtl/eqtl/johnson')
relevant_files_suffix = 'full_associations.csv.gz'
output_fn = 'merged_johnson_eqtl_results.tsv.gz'

# Iterate over all relevant files in the input directory and append them to the merged DataFrame
for region in regions:
    region_dir = input_files_dir / region
    # Initialize empty DF to hold merged results
    merged_df = pd.DataFrame()
    for input_file in region_dir.glob(f'*{relevant_files_suffix}'):
        fp = Path(input_file)
        print(f"Merging file: {fp}")
        df = coloc_transform.load_input_data(fp, header_lines=0)
        df['Region'] = region  # Add a new column for the region
        merged_df = pd.concat([merged_df, df], ignore_index=True)

        # Clean up memory
        del df
        gc.collect()
    
    print(f"Completed merging files for region: {region}")

    # Save the merged DataFrame to a new file
    output_fp = input_files_dir / f'{region}_{output_fn}'
    print(f"Saving merged DataFrame to: {output_fp}")
    merged_df.to_csv(output_fp, sep='\t', index=False, compression='gzip')
