# Imports
import dask.dataframe as dd
import pandas as pd
from pathlib import Path

# Configuration
tissue_sample_counts_fp = '/mnt/disks/working/locus_reports/qtl/eqtl/gtex/all_associations/tissue_sample_counts.tsv'
raw_dir = Path('/mnt/disks/working/locus_reports/qtl/sqtl/gtex/all_associations/GTEx_Analysis_v11_sQTL_all_associations')
prepared_dir = Path('/mnt/disks/working/locus_reports/qtl/sqtl/gtex/all_associations/GTEx_Analysis_v11_sQTL_all_associations_prepared')
gtex_tissue_manifest_fn = 'GTEx_tissue_manifest.tsv'
gene_symbol_mapping_fp = '/mnt/disks/working/dbs/ensembl/gene_id_mappings.tsv'
variant_id_mapping_fp = '/mnt/disks/working/locus_reports/qtl/eqtl/gtex/all_associations/GTEx_Analysis_2021-02-11_v11_WholeGenomeSeq_953Indiv.lookup_table.txt.gz'

# Function Defs

def process_partition(df, gene_symbol_map_dict, variant_id_map_dict):
    """
    Process a single partition of the tissue QTL data.

    Parameters:
        df (pd.DataFrame): A DataFrame containing a partition of the tissue QTL data.
        gene_symbol_map_dict (dict): Dictionary mapping gene IDs to gene symbols.
        variant_id_map_dict (dict): Dictionary mapping variant IDs to rsIDs.

    Returns:
        pd.DataFrame: A DataFrame containing the processed partition with additional columns for gene symbols and rsIDs.
    """
    # Create a copy of the DataFrame to avoid modifying the original
    df = df.copy()

    # Normalize gene_id
    df['gene_id'] = df['phenotype_id'].astype(str).str.split(':').str[-1].str.split('.').str[0]

    # Fast dict-based gene symbol lookup
    df['gene_symbol'] = df['gene_id'].map(gene_symbol_map_dict).fillna(df['gene_id'])

    # Fast dict-based variant ID lookup
    df['rsid'] = df['variant_id'].map(variant_id_map_dict).fillna(df['variant_id'])

    # Split variant_id
    variant_cols = df['variant_id'].astype(str).str.extract(
        r'^([^_]+)_([^_]+)_([^_]+)_([^_]+)(?:_(.+))?$',
        expand=True
    )
    variant_cols.columns = ['chr', 'pos', 'ref', 'alt', 'build']
    df[['chr', 'pos', 'ref', 'alt', 'build']] = variant_cols
    
    return df

def process_tissue(tissue, gene_symbol_map_dict, variant_id_map_dict):

    #  Lazy load all chromosome files that match the pattern
    pattern = str(raw_dir / f"{tissue}*.parquet")
    print(f"Reading parquet files for tissue: {tissue} with pattern: {pattern}", flush=True)
    ddf = dd.read_parquet(pattern, split_row_groups = False) #type: ignore 
    print(f"Successfully read parquet files for tissue: {tissue}. Number of partitions: {ddf.npartitions}", flush=True)

    # Check if ddf is None (no parquet files found)
    if ddf is None:
        print(f"  Skipping tissue {tissue} (no parquet files)", flush=True)
        return None

    # Process the Dask DataFrame using map_partitions with the provided dictionaries
    print(f"Processing {ddf.npartitions} partitions for tissue: {tissue}", flush=True)
    processed = ddf.map_partitions(process_partition, gene_symbol_map_dict=gene_symbol_map_dict, variant_id_map_dict=variant_id_map_dict)

    # Write the processed data to a parquet file
    tissue_output_dir = Path(f"{prepared_dir}/{tissue}")
    tissue_output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Writing {tissue} to {tissue_output_dir} as partitioned parquet files.", flush=True)
    
    for i, partition in enumerate(processed.to_delayed()):
            print(f"  Writing partition {i+1}/{processed.npartitions} for tissue {tissue}", flush=True)
            
            # Compute this partition (loads into memory) and write it
            df_part = partition.compute()
            
            # Write individual partition file
            output_file = tissue_output_dir / f"part_{i:04d}.parquet"
            df_part.to_parquet(
                str(output_file),
                index=False,
                compression='snappy',
                engine='pyarrow'
            )
            
            print(f"    Wrote {len(df_part)} rows to {output_file.name}", flush=True)
            
            # Explicitly free memory
            del df_part

    return tissue_output_dir


def main():
    print("Loading mappings...", flush=True)

    # Load tissue sample counts
    tissue_sample_counts = pd.read_table(tissue_sample_counts_fp)
    print(f"Loaded {len(tissue_sample_counts)} tissues", flush=True)

    # Load gene symbol mapping
    gene_symbol_map = pd.read_table(gene_symbol_mapping_fp)
    gene_symbol_map_dict = dict(zip(gene_symbol_map['gene_id'], gene_symbol_map['gene_name']))
    print(f"Loaded {len(gene_symbol_map_dict)} gene symbol mappings", flush=True)

    # Load variant ID mapping
    variant_id_map = pd.read_table(variant_id_mapping_fp)
    variant_id_map_dict = dict(zip(variant_id_map['variant_id'], variant_id_map['rs_id_dbSNP155_GRCh38p13']))
    print(f"Loaded {len(variant_id_map_dict)} variant ID mappings", flush=True)

    # Process each tissue WITHOUT a Dask distributed client
    # Dask will still parallelize partitions locally via threading/processes
    for tissue in tissue_sample_counts['Tissue']:
        print(f"Processing tissue: {tissue}", flush=True)

        # Check if any files match the pattern
        tissue_files = list(raw_dir.glob(f"{tissue}*.parquet"))

        if not tissue_files:
            print(f"  No parquet files found for tissue: {tissue}", flush=True)
            continue
        
        # Process the tissue and get the output directory
        tissue_output_dir = process_tissue(tissue, gene_symbol_map_dict, variant_id_map_dict)

        # Check if the output directory was created
        if tissue_output_dir:
            print(f"Finished processing tissue: {tissue}. Output written to {tissue_output_dir}", flush=True)
            tissue_sample_counts.loc[tissue_sample_counts['Tissue'] == tissue, 'file_path'] = str(tissue_output_dir)

    # Drop rows with missing file paths
    tissue_sample_counts = tissue_sample_counts.dropna(subset=['file_path'])

    # Save the updated tissue sample counts with file paths
    tissue_sample_counts.to_csv(Path(prepared_dir, gtex_tissue_manifest_fn), sep='\t', index=False)
    print("Done!", flush=True)


if __name__ == "__main__":
    main()