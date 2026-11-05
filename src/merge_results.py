# Imports
import argparse
import numpy as np # type: ignore (silences pylance warning)
import pandas as pd
from pathlib import Path
from helpers import load_input_data

# Define SuSiE Column names 
susie_cols = [
    'GWAS_Strata', 'TISSUE', 'LOCUS_ID', 'ld_panel_note', 'gene_id',
    'method', 'hit1', 'hit2', 'nsnps', 'PP.H0', 'PP.H1', 
    'PP.H2', 'PP.H3', 'PP.H4', 'cred_set_size_95', 'top_snp_h4'
]

# Merge Results
def merge_result(susie_df, coloc_abf_df, susie_cols=susie_cols):
    """
    Merge SuSiE and coloc ABF results on the 'locus_id' column.

    Parameters:
        susie_df (pd.DataFrame): DataFrame containing SuSiE results.
        coloc_abf_df (pd.DataFrame): DataFrame containing coloc ABF results.

    Returns:
        pd.DataFrame: Merged DataFrame containing both SuSiE and coloc ABF results.
    """
    # Ensure that all columns in susie_cols are present in coloc_abf_df
    missing_cols = [col for col in susie_cols if col not in coloc_abf_df.columns]

    # Add missing columns to coloc_abf_df with NaN values
    coloc_abf_df[missing_cols] = np.nan

    # Reorder the columns in coloc_abf_df to match susie_cols
    coloc_abf_df = coloc_abf_df[susie_cols]

    # Concatenate the two DataFrames
    merged_df = pd.concat([susie_df, coloc_abf_df], ignore_index=True, sort=False)
    
    return merged_df

# Def main function
def main(args):

    # Parse command line arguments
    susie_results_fp = Path(args.susie_results_fp)
    coloc_abf_results_fp = Path(args.coloc_abf_results_fp)
    merged_results_fp = Path(args.merged_results_fp)

    # Load results
    try:
        susie_df = load_input_data(susie_results_fp, 0)
    except Exception as e:
        print(f"Error loading SuSiE results from {susie_results_fp}: {e}")
        susie_df = pd.DataFrame(columns=susie_cols)  # Create an empty DataFrame if loading fails

    try:
        coloc_abf_df = load_input_data(coloc_abf_results_fp, 0)
    except Exception as e:
        print(f"Error loading coloc ABF results from {coloc_abf_results_fp}: {e}")
        coloc_abf_df = pd.DataFrame(columns=susie_cols)  # Create an empty DataFrame if loading fails

    # Merge results
    merged_df = merge_result(susie_df, coloc_abf_df)

    # Save merged results
    merged_df.to_csv(merged_results_fp, sep='\t', index=False)


# CLI
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge SuSiE and coloc ABF results.")
    parser.add_argument("--susie_results_fp", type=str, required=True, help="File path to SuSiE results.")
    parser.add_argument("--coloc_abf_results_fp", type=str, required=True, help="File path to coloc ABF results.")
    parser.add_argument("--merged_results_fp", type=str, required=True, help="File path to save merged results.")
    args = parser.parse_args()
    main(args)
