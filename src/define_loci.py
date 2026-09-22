# Imports
import argparse
import pandas as pd
from typing import Optional
from pathlib import Path
from helpers import load_input_data, normalize_chromosome

# -------------------------- Helper Function Definitions -------------------------

def _standardize_df(df: pd.DataFrame, standard_keys: list[str], standard_key_names: list[str]) -> pd.DataFrame:
    """
    Standardizes the DataFrame to include only the specified columns.

    Args:
        df (pd.DataFrame): The DataFrame to be standardized.
        standard_keys (list[str]): List of column names to retain in the standardized DataFrame.
        standard_key_names (list[str]): List of new column names for the standardized DataFrame.

    Returns:
        pd.DataFrame: The standardized DataFrame containing only the specified columns.

    """
    # Create a copy of the DataFrame with only the standard columns
    standard_df = df[standard_keys].copy()
    
    # Rename columns to standard names
    standard_df.columns = standard_key_names

    return standard_df

# -------------------------- Class Definition -------------------------

class LocusRanges:
    """
    A class to define locus ranges based on SNP identifiers and flank size.

    """
    def __init__(self, args: argparse.Namespace) -> None:
        self.df = load_input_data(fp=Path(args.input_file), header_lines=args.header_lines)
        print(f"Loaded GWAS summary statistics from {args.input_file} with {len(self.df)} SNPs.")
        self.snp_id_key = args.snp_id_key
        self.chr_key = args.chr_key
        self.pos_key = args.pos_key
        self.p_key = args.p_key
        self.p_threshold = args.p_threshold
        self.flank_size = args.flank_size
        self.left_bound_key = f"LEFT_{self.flank_size // 1000}KB"
        self.right_bound_key = f"RIGHT_{self.flank_size // 1000}KB"
        self.locus_id_key = args.locus_id_key
        self.merged_df = pd.DataFrame()
        self.n_snps_key = "N_SNPS"
        self.standardized_chr_key = args.standardized_chr_key
        self.standard_keys = [self.chr_key, self.left_bound_key, self.right_bound_key, self.locus_id_key]
        self.standard_key_names = [self.standardized_chr_key, self.left_bound_key, self.right_bound_key, self.locus_id_key]


    def filter_by_p_value(self) -> None:
        """
        Filters the DataFrame to include only rows with p-values below the specified threshold.

        Returns:
            None: The DataFrame is modified in place to include only significant SNPs.

        """
        if self.p_key in self.df.columns:
            self.df = self.df[self.df[self.p_key] < self.p_threshold]
        else:
            raise ValueError(f"ERROR: P-value column '{self.p_key}' not found in the input DataFrame.")


    def parse_snp_id(self) -> None:
        """
        Extracts chromosome and position from a SNP identifier string. By default, it assumes the SNP identifier is in the format "chr:pos" or "chr_pos". If the format is not recognized, it returns None for each value.

        Returns:
            None: The DataFrame is modified in place to include the extracted chromosome and position.

        """
        # Extract chromosome and position from SNP identifier
        snps = pd.Series(self.df[self.snp_id_key].astype(str))

        # Split the SNP identifier by ':' or '_' and extract chromosome and position
        self.df[self.chr_key], self.df[self.pos_key] = zip(*[(parts[0], int(parts[1])) if len(parts) > 1 else (None, None) for parts in snps.str.split('[:_]')])


    def define_range(self) -> None:
        """
        Defines a range around a given position in the DataFrame, based on provided flank value.

        Returns:
            None: The DataFrame is modified in place to include the extracted chromosome and position.

        """
        # Define left and right boundaries based on position and flank size
        self.df[self.left_bound_key], self.df[self.right_bound_key] = zip(*[(pos - self.flank_size, pos + self.flank_size) if pd.notnull(pos) else (None, None) for pos in self.df[self.pos_key]])


    def merge_overlapping(self) -> None:
        """
        Merges overlapping ranges in the DataFrame based on left and right boundaries.

        Returns:
            None: The DataFrame is modified in place to include the merged ranges and associated SNPs.

        """
        # Sort by left boundary
        df_sorted = self.df.sort_values(by=[self.chr_key, self.left_bound_key]).reset_index(drop=True)

        # Identify new groups based on overlapping ranges
        cummax_end = (
            df_sorted
            .groupby(self.chr_key, sort=False)[self.right_bound_key]
            .transform(lambda s: s.shift(1).cummax())
        )
        new_group = (df_sorted[self.left_bound_key] > cummax_end).fillna(True)

        # Build a global group ID that resets per chromosome
        chr_change = df_sorted[self.chr_key] != df_sorted[self.chr_key].shift(1).fillna(df_sorted[self.chr_key].iloc[0])
        group_id = (new_group | chr_change).cumsum()

        # Define aggregation functions for merging
        aggregation_functions = {
            self.chr_key: pd.NamedAgg(column=self.chr_key, aggfunc="first"),
            self.left_bound_key: pd.NamedAgg(column=self.left_bound_key, aggfunc="min"),
            self.right_bound_key: pd.NamedAgg(column=self.right_bound_key, aggfunc="max"),
            self.snp_id_key: pd.NamedAgg(column=self.snp_id_key, aggfunc=lambda x: ";".join(x.astype(str))),
        }

        # Aggregate per group
        self.merged_df = df_sorted.groupby(group_id, sort=False).agg(
            **aggregation_functions
        ).reset_index(drop=True)

        # Derived columns
        self.merged_df[self.left_bound_key] = self.merged_df[self.left_bound_key].astype(int)
        self.merged_df[self.right_bound_key] = self.merged_df[self.right_bound_key].astype(int)
        self.merged_df[self.n_snps_key] = self.merged_df[self.snp_id_key].str.count(";") + 1

        # Drop SNPs column
        self.merged_df.drop(columns=[self.snp_id_key], inplace=True)


    def add_locus_id(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Adds a locus ID to the DataFrame based on chromosome and boundaries.

        Args:
            df (pd.DataFrame): The DataFrame containing the data.

        Returns:
            pd.DataFrame: The DataFrame with the new locus ID column added.

        """
        # Create a locus ID by concatenating chromosome, left boundary, and right boundary
        df[self.locus_id_key] = (
            df[self.chr_key].astype(str) + ":" +
            df[self.left_bound_key].astype(str) + "-" +
            df[self.right_bound_key].astype(str)
        )
        return df


    def define_locus_ranges(self) -> tuple[pd.DataFrame, pd.DataFrame]:
        """
        Main function to parse SNP IDs, define ranges, and merge overlapping ranges.

        Returns:
            tuple[pd.DataFrame, pd.DataFrame]: Tuple containing the merged DataFrame with loci information and the original DataFrame with added columns for chromosome, position, and boundaries.

        """

        # Filter the DataFrame based on the p-value threshold
        self.filter_by_p_value()
        print(f"Significant SNPs after filtering by p-value < {self.p_threshold}: {len(self.df)}")

        # Parse SNP IDs to extract chromosome and position, if they are not already present in the DataFrame. Raise an error if neither SNP ID column nor both chromosome and position columns are present.
        if self.chr_key not in self.df.columns and self.pos_key not in self.df.columns:
            if self.snp_id_key is not None:
                self.parse_snp_id()
            else:
                raise ValueError("ERROR: Unable to define loci. Either SNP ID column must be provided or both chromosome and position columns must be present in the input file.")

        # Normalize chromosome representation
        self.df[self.chr_key] = self.df[self.chr_key].astype(str).map(normalize_chromosome)

        # Define the range around each SNP based on the flank size
        self.define_range()

        # Merge overlapping ranges to define loci
        self.merge_overlapping()

        # Add locus IDs
        self.df = self.add_locus_id(self.df)
        self.merged_df = self.add_locus_id(self.merged_df)

        # Standardize the DataFrames to include only relevant columns
        full_df = _standardize_df(self.df, self.standard_keys, self.standard_key_names)
        merged_df = _standardize_df(self.merged_df, self.standard_keys + [self.n_snps_key], self.standard_key_names + [self.n_snps_key])

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
    locus_definer = LocusRanges(args)

    # Define loci
    full_df, merged_df = locus_definer.define_locus_ranges()

    # Print number of loci defined
    print(f"Number of loci defined: {len(full_df)}")
    print(f"Number of merged loci: {len(merged_df)}")

    # Save the results to CSV files
    full_df.to_csv(f'{args.output_file_prefix}.tsv', sep='\t', index=False)
    if len(merged_df) < len(full_df):
        output_fp = f'{args.output_file_prefix}.merged.tsv'
        print(f"Saving merged loci to {output_fp}")
        merged_df.to_csv(output_fp, sep='\t', index=False)


# -------------------------- Command-Line Interface -------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Define locus ranges based on SNP identifiers and flank size.")
    parser.add_argument("--input_file", type=str, required=True, help="Path to the input CSV file containing SNP data.")
    parser.add_argument("--header_lines", type=int, default=0, help="Number of header lines to skip in the input file.")
    parser.add_argument("--output_file_prefix", type=str, required=True, help="Prefix for the output CSV files, including the path and base name.")
    parser.add_argument("--snp_id_key", type=str, default="SNP", help="Column name for SNP identifiers in the input DataFrame.")
    parser.add_argument("--chr_key", type=str, required=True, help="Column name for chromosome information. If not provided, it will be extracted from SNP_ID.")
    parser.add_argument("--pos_key", type=str, required=True, help="Column name for position information. If not provided, it will be extracted from SNP_ID.")
    parser.add_argument("--p_key", type=str, default="P", help="Column name for p-values in the input DataFrame.")
    parser.add_argument("--p_threshold", type=float, default=5e-6, help="P-value threshold for defining significant loci.")
    parser.add_argument("--flank_size", type=int, default=500000, help="Flank size in base pairs to define the range around each SNP.")
    parser.add_argument("--standardized_chr_key", default="CHR", help="Column name for standardized chromosome information in the output DataFrame.")
    parser.add_argument("--locus_id_key", type=str, default="LOCUS_ID", help="Column name for locus IDs in the output DataFrame.")

    args = parser.parse_args()
    main(args)