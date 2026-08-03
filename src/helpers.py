# Imports
import pandas as pd
from pathlib import Path

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