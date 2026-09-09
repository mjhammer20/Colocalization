# Imports
import shutil
import dask.dataframe as dd
import pandas as pd
from pathlib import Path

# Configuration
tissue_sample_counts_fp = '/mnt/working/locus_reports/qtl/eqtl/gtex/all_associations/tissue_sample_counts.tsv'
raw_dir = Path('/mnt/working/locus_reports/qtl/eqtl/gtex/all_associations/GTEx_Analysis_v11_eQTL_all_associations')
prepared_dir = Path('/mnt/working/locus_reports/qtl/eqtl/gtex/all_associations/GTEx_Analysis_v11_eQTL_all_associations_prepared')
gtex_tissue_manifest_fn = 'GTEx_tissue_manifest.tsv'
gene_symbol_mapping_fp = '/mnt/working/dbs/ensembl/gene_id_mappings.tsv'

# Function Defs

def process_partition(df):
    """
    Process a single partition of the tissue QTL data.

    Parameters:
        df (pd.DataFrame): A DataFrame containing a partition of the tissue QTL data.

    Returns:
        pd.DataFrame: A DataFrame containing the processed partition with additional columns for gene symbols and
    """
    # Create a copy of the DataFrame to avoid modifying the original
    df = df.copy()

    # Normalize gene_id
    df['gene_id'] = df['gene_id'].astype(str).str.split('.').str[0]
    # Fast dict-based gene symbol lookup

    df['gene_symbol'] = df['gene_id'].map(gene_symbol_map_dict).fillna(df['gene_id'])

    # Split variant_id
    variant_cols = df['variant_id'].astype(str).str.extract(
        r'^([^_]+)_([^_]+)_([^_]+)_([^_]+)(?:_(.+))?$',
        expand=True
    )
    variant_cols.columns = ['chr', 'pos', 'ref', 'alt', 'build']
    df[['chr', 'pos', 'ref', 'alt', 'build']] = variant_cols
    
    return df

def process_tissue(tissue):

    #  Lazy load all chromosome files that match the pattern
    pattern = str(raw_dir / f"{tissue}*.parquet")
    print(f"Reading parquet files for tissue: {tissue} with pattern: {pattern}")
    ddf = dd.read_parquet(pattern, split_row_groups = False) #type: ignore 
    print(f"Successfully read parquet files for tissue: {tissue}. Number of partitions: {ddf.npartitions}")

    # Check if ddf is None (no parquet files found)
    if ddf is None:
        print(f"  Skipping tissue {tissue} (no parquet files)")
        return None

    # Apply transformation to all partitions lazily
    processed = ddf.map_partitions(process_partition)

    # Write the processed data to a parquet file
    temp_output_dir= Path(f"{prepared_dir}/.{tissue}")
    print(f"Writing {tissue} to {temp_output_dir} as partitioned parquet files.")
    processed.to_parquet(str(temp_output_dir), write_index=False, compression='snappy', engine='pyarrow')

    return temp_output_dir


def main():
    # Load tissue sample counts
    tissue_sample_counts = pd.read_table(tissue_sample_counts_fp)

    # Load gene symbol mapping
    gene_symbol_map = pd.read_table(gene_symbol_mapping_fp)
    global gene_symbol_map_dict
    gene_symbol_map_dict = dict(zip(gene_symbol_map['gene_id'], gene_symbol_map['gene_name']))

    # Process each tissue in the tissue sample counts DataFrame
    for tissue in tissue_sample_counts['Tissue']:
        print(f"Processing tissue: {tissue}")

        # Check if any files match the pattern
        tissue_files = list(raw_dir.glob(f"{tissue}*.parquet"))

        if not tissue_files:
            print(f"  No parquet files found for tissue: {tissue}")
            continue
        
        # Process the tissue and get the temporary output directory
        temp_tissue_output_dir = process_tissue(tissue)

        # Check if the temporary output directory was created
        if temp_tissue_output_dir:
            
            # Rename part files based on chromosome from original filenames
            tissue_output_dir = Path(f"{prepared_dir}/{tissue}")
            tissue_output_dir.mkdir(parents=True, exist_ok=True)
            
            for i, orig_file in enumerate(tissue_files):
                # Extract chromosome from filename (e.g., "Adipose_Subcutaneous.chr1.v8.sqtl_eur_all.parquet" -> "1")
                chr_str = orig_file.stem.split('.')[3]
                part_file = sorted(list(temp_tissue_output_dir.glob("part-*.parquet")))[i]
                new_name = tissue_output_dir / f"{chr_str}.parquet"
                shutil.move(str(part_file), str(new_name))
                print(f"  Renamed partition {i} to {new_name.name}")
            
            # Clean up temp dir
            shutil.rmtree(temp_tissue_output_dir)

            print(f"Finished processing tissue: {tissue}. Output written to {tissue_output_dir}")

            # Add file path to DataFrame
            tissue_sample_counts.loc[tissue_sample_counts['Tissue'] == tissue, 'file_path'] = str(tissue_output_dir)

    # Drop rows with missing file paths
    tissue_sample_counts = tissue_sample_counts.dropna(subset=['file_path'])

    # Save the updated tissue sample counts with file paths
    tissue_sample_counts.to_csv(Path(prepared_dir, gtex_tissue_manifest_fn), sep='\t', index=False)

if __name__ == "__main__":
    main()