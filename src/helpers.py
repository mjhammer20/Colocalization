# Imports
import pandas as pd
from pathlib import Path
from typing import Optional
import numpy as np

# Function to load input data
def load_input_data(
        fp: Path,
        header_lines: int
    ) -> pd.DataFrame:
    """
    Load input data from a file into a pandas DataFrame.

    Args:
        fp (Path): Path to the input file.
        header_lines (int): Number of header lines to skip.

    Returns:
        pd.DataFrame: The loaded DataFrame.
    """

    # Check if the input file exists
    if not fp.exists():
        raise FileNotFoundError(f"Input file {fp} does not exist.")

    try:
        fn = fp.name.lower()
        if fn.endswith('.csv'):
            df = pd.read_csv(fp)
        elif fn.endswith('.csv.gz'):
            df = pd.read_csv(fp, compression='gzip')
        elif fn.endswith('.tsv'):
            df = pd.read_csv(fp, sep='\t')
        elif fn.endswith('.tsv.gz'):
            df = pd.read_csv(fp, sep='\t', compression='gzip')
        elif fn.endswith('.txt'):
            df = pd.read_table(fp, sep='\t', skiprows=header_lines)
            if df.shape[1] == 1:
                df = pd.read_table(fp, sep=',', skiprows=header_lines)
        elif fn.endswith('.txt.gz'):
            df = pd.read_table(fp, sep='\t', compression='gzip', skiprows=header_lines)
            if df.shape[1] == 1:
                df = pd.read_table(fp, sep=',', compression='gzip', skiprows=header_lines)
        elif fn.endswith('.xlsx') or fn.endswith('.xls'):
            df = pd.read_excel(fp, skiprows=header_lines)
        else:
            raise ValueError(f"Unsupported file format for input file {fp}. Supported formats are: .csv, .csv.gz, .tsv, .tsv.gz, .txt, .txt.gz, .xlsx, .xls")

    except Exception as e:
        raise Exception(f"Error loading input file {fp}: {e}")

    return df

# Function to normalize chromosome representation
def normalize_chromosome(chromosome: str) -> str:
    """
    Normalize chromosome representation to a standard format.

    Args:
        chromosome (str): The chromosome identifier (e.g., 'chr1', '1', 'chrX', 'X').

    Returns:
        str: Normalized chromosome representation (e.g., '1', 'X').
    """

    if chromosome.startswith(('chr', )):
        return chromosome[3:]
    
    return chromosome

# Function to create a marker ID in the format chr:bp:ref:alt for each row in the DataFrame
def add_variant_id(
        df: pd.DataFrame,
        chr: pd.Series,
        pos: pd.Series,
        ref: pd.Series,
        alt: pd.Series,
        variant_id_key: str
    ) -> pd.DataFrame:
    """
    Create a marker ID in the format chr:bp:ref:alt for each row in the DataFrame.

    Args:
        df (pd.DataFrame): The input DataFrame.
        chr_col (str): The name of the column containing chromosome information.
        pos_col (str): The name of the column containing position information.
        ref_allele_col (str): The name of the column containing reference allele information.
        alt_allele_col (str): The name of the column containing alternate allele information.
        variant_id_col (str): The name of the column to store the marker ID.

    Returns:
        pd.DataFrame: The DataFrame with an additional column for marker IDs.

    """

    # Build the variant ID as chr:pos:ref:alt; set to "NA" where any component is missing
    all_present = pos.notna() & ref.notna() & alt.notna()
    df[variant_id_key] = np.where(
        all_present,
        chr + ":" + pos.astype(str) + ":" + ref.astype(str) + ":" + alt.astype(str),
        "NA"
    )

    
    return df

def numeric_series(df: pd.DataFrame, colname: Optional[str]) -> pd.Series:
    """
    Convert a specified column in a DataFrame to numeric, coercing errors to NaN.

    Args:
        df (pd.DataFrame): The input DataFrame.
        colname (Optional[str]): The name of the column to convert.

    Returns:
        pd.Series: A pandas Series containing the converted values.

    """
    # Check if the column name is valid and exists in the DataFrame
    if not colname or colname not in df.columns:
        return pd.Series(np.nan, index=df.index)

    # Convert the specified column to numeric, coercing errors to NaN
    ser = pd.to_numeric(df[colname], errors="coerce")

    # Return the Series, ensuring it has the same index as the DataFrame
    if isinstance(ser, pd.Series):
        return ser
    
    # If the result is not a Series (e.g., if it's a DataFrame), return a Series of NaN values
    return pd.Series(ser, index=df.index)