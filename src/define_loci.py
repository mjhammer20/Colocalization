# Imports
import argparse
import pandas as pd
from typing import Optional
from pathlib import Path
from helpers import load_input_data

# -------------------------- Helper Function Definitions -------------------------

def _add_locus_id(df: pd.DataFrame, chr_col: str, left_bound_col: str, right_bound_col: str, locus_id_col: str) -> pd.DataFrame:
    """
    Adds a locus ID to the DataFrame based on chromosome and boundaries.

    Args:
        df (pd.DataFrame): The DataFrame containing the data.
        chr_col (str): The name of the column containing chromosome information.
        left_bound_col (str): The name of the column containing left boundary information.
        right_bound_col (str): The name of the column containing right boundary information.
        locus_id_col (str): The name of the new column to be created for locus IDs.

    Returns:
        pd.DataFrame: The DataFrame with the new locus ID column added.
    """
    df[locus_id_col] = (
        df[chr_col].astype(str) + ":" +
        df[left_bound_col].astype(str) + "-" +
        df[right_bound_col].astype(str)
    )
    return df

def _standardize_df(df: pd.DataFrame, standard_cols: list[str], standard_col_names: list[str]) -> pd.DataFrame:
    """
    Standardizes the DataFrame to include only the specified columns.

    Args:
        df (pd.DataFrame): The DataFrame to be standardized.
        standard_cols (list[str]): List of column names to retain in the standardized DataFrame.

    Returns:
        pd.DataFrame: The standardized DataFrame containing only the specified columns.
    """
    return pd.DataFrame(df[standard_cols].copy(), columns=standard_col_names)

# -------------------------- Class Definition -------------------------

class LocusRanges:
    def __init__(self, input_fp: Path, header_lines: int, chr_col: Optional[str], pos_col: Optional[str], snp_id_col: Optional[str], flank_size: int = 500000):
        self.df = load_input_data(input_fp, header_lines=header_lines)
        self.snp_id_col = snp_id_col
        self.chr_col = chr_col or "CHR"
        self.pos_col = pos_col or "BP"
        self.flank_size = flank_size
        self.left_bound_col = f"LEFT_{flank_size // 1000}KB"
        self.right_bound_col = f"RIGHT_{flank_size // 1000}KB"
        self.locus_id_col = "LOCUS_ID"
        self.merged_df = pd.DataFrame()
        self.n_snps_col = "N_SNPS"
        self.standard_cols = [self.chr_col, self.left_bound_col, self.right_bound_col, self.locus_id_col]
        self.standard_col_names = ["CHR", self.left_bound_col, self.right_bound_col, self.locus_id_col]

    def parse_snp_id(self) -> None:
        """
        Extracts chromosome and position from a SNP identifier string. By default, it assumes the SNP identifier is in the format "chr:pos" or "chr_pos". If the format is not recognized, it returns None for each value.

        Returns:
            None: The DataFrame is modified in place to include the extracted chromosome and position.
        """
        snps = pd.Series(self.df[self.snp_id_col].astype(str))
        self.df[self.chr_col], self.df[self.pos_col] = zip(*[(parts[0], int(parts[1])) if len(parts) > 1 else (None, None) for parts in snps.str.split('[:_]')])


    def define_range(self) -> None:
        """
        Defines a range around a given position in the DataFrame, based on provided flank value.

        Returns:
            None: The DataFrame is modified in place to include the extracted chromosome and position.
        """
        self.df[self.left_bound_col], self.df[self.right_bound_col] = zip(*[(pos - self.flank_size, pos + self.flank_size) if pd.notnull(pos) else (None, None) for pos in self.df[self.pos_col]])


    def merge_overlapping(self) -> None:
        """
        Merges overlapping ranges in the DataFrame based on left and right boundaries.

        Returns:
            None: The DataFrame is modified in place to include the merged ranges and associated SNPs.
        """
        # Sort by left boundary
        df_sorted = self.df.sort_values(by=[self.chr_col, self.left_bound_col]).reset_index(drop=True)


        cummax_end = (
            df_sorted
            .groupby(self.chr_col, sort=False)[self.right_bound_col]
            .transform(lambda s: s.shift(1).cummax())
        )
        new_group = (df_sorted[self.left_bound_col] > cummax_end).fillna(True)

        # Build a global group ID that resets per chromosome
        chr_change = df_sorted[self.chr_col] != df_sorted[self.chr_col].shift(1).fillna(df_sorted[self.chr_col].iloc[0])
        group_id = (new_group | chr_change).cumsum()

        # Aggregate per group
        self.merged_df = df_sorted.groupby(group_id, sort=False).agg(
            chr=(self.chr_col, "first"),
            left_boundary=(self.left_bound_col, "min"),
            right_boundary=(self.right_bound_col, "max"),
            SNPs=(self.snp_id_col, lambda x: ";".join(x.astype(str)))
        ).reset_index(drop=True)

        # Derived columns
        self.merged_df[self.left_bound_col] = self.merged_df["left_boundary"].astype(int)
        self.merged_df[self.right_bound_col] = self.merged_df["right_boundary"].astype(int)
        self.merged_df[self.n_snps_col] = self.merged_df["SNPs"].str.count(";") + 1

        # Drop SNPs column
        self.merged_df.drop(columns=["SNPs"], inplace=True)

    def define_locus_ranges(self) -> tuple[pd.DataFrame, pd.DataFrame]:
        """
        Main function to parse SNP IDs, define ranges, and merge overlapping ranges.

        Returns:
            tuple[pd.DataFrame, pd.DataFrame]: Tuple containing the merged DataFrame with loci information and the original DataFrame with added columns for chromosome, position, and boundaries.
        """
        
        # Parse SNP IDs to extract chromosome and position, if they are not already present in the DataFrame. Raise an error if neither SNP ID column nor both chromosome and position columns are present.
        if self.chr_col not in self.df.columns and self.pos_col not in self.df.columns:
            if self.snp_id_col is not None:
                self.parse_snp_id()
            else:
                raise ValueError("ERROR: Unable to define loci. Either SNP ID column must be provided or both chromosome and position columns must be present in the input file.")

        # Define the range around each SNP based on the flank size
        self.define_range()

        # Merge overlapping ranges to define loci
        self.merge_overlapping()

        # Add locus IDs
        self.df = _add_locus_id(self.df, self.chr_col, self.left_bound_col, self.right_bound_col, self.locus_id_col)
        self.merged_df = _add_locus_id(self.merged_df, self.chr_col, self.left_bound_col, self.right_bound_col, self.locus_id_col)

        # Standardize the DataFrames to include only relevant columns
        full_df = _standardize_df(self.df, self.standard_cols, self.standard_col_names)
        merged_df = _standardize_df(self.merged_df, self.standard_cols + [self.n_snps_col], self.standard_col_names + [self.n_snps_col])

        return full_df, merged_df
    

# -------------------------- Main Function -------------------------

def main(args: argparse.Namespace) -> None:
    """
    Main function to execute the locus definition process.

    Args:
        args (argparse.Namespace): Command-line arguments containing input and output file paths, SNP ID column name, chromosome column name, position column name, and flank size.

    Returns:
        None
    """

    # Create a LocusRanges object
    locus_definer = LocusRanges(input_fp = Path(args.input_file), header_lines=args.header_lines, chr_col = args.chr_col, pos_col = args.pos_col, snp_id_col = args.snp_id_col, flank_size = args.flank_size)

    # Define loci
    full_df, merged_df = locus_definer.define_locus_ranges()

    # Save the results to CSV files
    full_df.to_csv(f'{args.output_file_prefix}.tsv', sep='\t', index=False)
    if len(merged_df) < len(full_df):
        merged_df.to_csv(f'{args.output_file_prefix}.merged.tsv', sep='\t', index=False)


# -------------------------- Command-Line Interface -------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Define locus ranges based on SNP identifiers and flank size.")
    parser.add_argument("--input_file", required=True, help="Path to the input CSV file containing SNP data.")
    parser.add_argument("--header_lines", type=int, default=0, help="Number of header lines to skip in the input file.")
    parser.add_argument("--output_file_prefix", required=True, help="Prefix for the output CSV files, including the path and base name.")
    parser.add_argument("--snp_id_col", default="SNP", help="Column name for SNP identifiers in the input DataFrame.")
    parser.add_argument("--chr_col", default=None, help="Column name for chromosome information. If not provided, it will be extracted from SNP_ID.")
    parser.add_argument("--pos_col", default=None, help="Column name for position information. If not provided, it will be extracted from SNP_ID.")
    parser.add_argument("--flank_size", type=int, default=500000, help="Flank size in base pairs to define the range around each SNP.")

    args = parser.parse_args()
    main(args)