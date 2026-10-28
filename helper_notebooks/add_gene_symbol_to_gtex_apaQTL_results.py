# Imports
import pandas as pd

# Configuration
gene_symbol_mapping_fp = '/mnt/working/dbs/ensembl/gene_id_mappings.tsv'
gtex_apaQTL_results_fp = '/mnt/working/locus_reports/qtl/apaQTL/gtex/results/meta6_pd_full_gtex_v11_apaQTLs_sub_all_results_merged.tsv'

# Functions
def extract_gene_id_from_phenotype_id(phenotype_id):
    """
    Split the phenotype_id into gene_id and transcript_id.

    Parameters:
        phenotype_id (str): The phenotype_id string to be split.

    Returns:
        gene_id (str): The extracted gene_id from the phenotype_id.
    """
    gene_id = phenotype_id.split('.')[0]
    return gene_id

def map_gene_symbols(df, gene_symbol_map_dict):
    """
    Map gene symbols to the DataFrame based on gene_id.

    Parameters:
        df (pd.DataFrame): The DataFrame containing the QTL results.
        gene_symbol_map_dict (dict): Dictionary mapping gene IDs to gene symbols.

    Returns:
        pd.DataFrame: The DataFrame with an additional column for gene symbols.
    """
    # Initialize phenotype_id column (used as gene_id for colocalization)
    df['phenotype_id'] = df['gene_id']

    # Drop gene_id column prior to id extraction
    df.drop(columns='gene_id')

    # Extract gene_id from phenotype_id
    df['gene_id'] = df['phenotype_id'].apply(extract_gene_id_from_phenotype_id)

    # Map gene symbols using the provided dictionary
    df['gene_symbol'] = df['gene_id'].map(gene_symbol_map_dict).fillna(df['gene_id'])

    return df

def main():
    # Load gene symbol mapping
    gene_symbol_map = pd.read_table(gene_symbol_mapping_fp)
    gene_symbol_map_dict = dict(zip(gene_symbol_map['gene_id'], gene_symbol_map['gene_name']))
    print(f"Loaded {len(gene_symbol_map_dict)} gene symbol mappings", flush=True)

    # Load GTEx APAQTL results
    gtex_apaQTL_results_df = pd.read_csv(gtex_apaQTL_results_fp, sep='\t')

    # Map gene symbols to the results DataFrame
    gtex_apaQTL_results_with_symbols_df = map_gene_symbols(gtex_apaQTL_results_df, gene_symbol_map_dict)

    # Save the updated DataFrame with gene symbols
    output_fp = gtex_apaQTL_results_fp.replace('.tsv', '_w_gene_symbols.tsv')
    gtex_apaQTL_results_with_symbols_df.to_csv(output_fp, sep='\t', index=False)
    print(f"Updated GTEx apaQTL results saved to: {output_fp}")

# Command Line Interface 
if __name__ == "__main__":
    main()