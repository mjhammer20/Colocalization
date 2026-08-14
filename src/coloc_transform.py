# Imports
import gc
import time
import xmltodict
import argparse
import json
import pandas as pd
import numpy as np #type: ignore (silences pylance warning)
from functools import lru_cache
from typing import Optional
from scipy.stats import norm #type: ignore (silences pylance warning)
from Bio import Entrez #type: ignore (silences pylance warning)
from pathlib import Path
from http.client import IncompleteRead
from urllib.error import HTTPError, URLError
from xml.parsers.expat import ExpatError
from helpers import load_input_data, normalize_chromosome

# ------------------------- Helper Functions -------------------------

def _numeric_series(df: pd.DataFrame, colname: Optional[str]) -> pd.Series:
    """
    Convert a specified column in a DataFrame to numeric, coercing errors to NaN.

    Args:
        df (pd.DataFrame): The input DataFrame.
        colname (Optional[str]): The name of the column to convert.

    Returns:
        pd.Series: A pandas Series containing the converted values.

    """
    if not colname or colname not in df.columns:
        return pd.Series(np.nan, index=df.index)

    ser = pd.to_numeric(df[colname], errors='coerce')

    # Ensure we always return a pandas Series (pd.to_numeric may return a scalar)
    if isinstance(ser, pd.Series):
        return ser
    return pd.Series(ser, index=df.index)

# ------------------------- Main Functions -------------------------

# Function to load expanded locus ranges from file
def load_expanded_ranges(fp: Path) -> dict:
    """
    Load expanded ranges from a file into a dictionary.

    Args:
        fp (Path): Path to the input file containing expanded ranges in JSON format.
    
    Returns:
        dict: A dictionary with chromosome as keys and a list of [start, end] ranges as values.

    """

    with open(fp, 'r') as f:
        expanded_ranges = json.load(f)

    return expanded_ranges


# Function to filter DataFrame to include only rows that fall within expanded locus ranges
def filter_by_expanded_ranges(
        df: pd.DataFrame,
        chr_col: str,
        pos_col: str,
        expanded_ranges: dict
    ) -> pd.DataFrame:
    """
    Filter a DataFrame to include only rows that fall within specified expanded ranges.

    Args:
        df (pd.DataFrame): The input DataFrame.
        chr_col (str): The name of the column containing chromosome information.
        pos_col (str): The name of the column containing position information.
        expanded_ranges (dict): A dictionary with chromosome as keys and a list of [start, end] ranges as values.

    Returns:
        pd.DataFrame: The filtered DataFrame containing only rows within the specified ranges.

    """

    # Initialize a boolean mask for filtering
    mask = pd.Series(False, index=df.index)

    # Iterate over each chromosome and its corresponding ranges
    for chrom, ranges in expanded_ranges.items():

        # Extract the indices of rows corresponding to the current chromosome
        idx = df.index[df[chr_col] == chrom]

        # Extract the position values for the current chromosome
        pos_vals = df[pos_col].loc[idx]  # type: pd.Series

        # Create a boolean mask for the current chromosome based on the specified ranges
        chrom_mask = pd.Series(False, index=idx)

        # Iterate over each range and update the mask to include positions within the range
        for start, end in ranges:
            chrom_mask |= pos_vals.between(start, end)

        # Update the overall mask with the chromosome-specific mask
        mask.loc[idx] = chrom_mask

    # Return the filtered DataFrame based on the mask
    filtered_df = pd.DataFrame(df[mask])
    return filtered_df


# Function to fetch rsID from NCBI dbSNP based on chromosome and position
@lru_cache(maxsize=65536)
def fetch_rsid(
        chromosome: str,
        position: int
    ) -> str:
    """
    Query NCBI dbSNP to retrieve rsID information based on chromosome and position.

    Args:
        chromosome (str): The chromosome number or identifier.
        position (int): The base pair position on the chromosome.

    Returns:
        Optional[str]: The corresponding rsID if found, otherwise None.

    """

    # Define the search term for the NCBI dbSNP query using the chromosome and position
    search_term = f"{chromosome}[Chromosome] AND {position}[Base Position] AND Homo sapiens[Organism]"

    # Execute the search query
    handle = Entrez.esearch(db="snp", term=search_term)

    # Parse the XML response
    response = xmltodict.parse(handle) # type: ignore (silences pylance warning)

    # Close the handle
    handle.close()

    # Extract rsIDs from the response
    rs_ids = response.get("IdList", None)

    # Return last rsID if available
    if not rs_ids:
        return "NA"
    
    return f'rs{str(rs_ids[-1]).strip()}'


# Function to fetch multiple rsID records from NCBI dbSNP
def fetch_rsid_records(rsids: list) -> dict:
    """
    Query NCBI dbSNP to retrieve allele information based on rsID.

    Args:
        rsids (list): The list of rsIDs to fetch.

    Returns:
        dict: Dictionary containing the fetched rsID records.
    """


    rsids_str = ','.join(rsids)

    # Retry fetching records up to 5 times (in event of transient errors)
    attempts = 0
    while attempts < 5:
        attempts += 1
        handle = None
        try:

            # Fetch the rsID records from NCBI dbSNP
            handle = Entrez.efetch(db="snp", id=rsids_str, retmode="xml")

            # Parse the XML response into a dictionary
            response = xmltodict.parse(handle) # type: ignore (silences pylance warning)

            # Close the handle
            handle.close()

            # Extract the records from the response
            records = response.get('ExchangeSet', {}).get('DocumentSummary', {})

        except (IncompleteRead, HTTPError, URLError, ExpatError) as e:
            
            # Transient error -> retry
            if handle is not None:
                    handle.close()
            if attempts < 5:
                time.sleep(2 * attempts)
                continue
            print(f"fetch_rsid_records: reached max retries, returning empty result for rsids: {rsids[0]} - {rsids[-1]}")
            return {}

        except Exception as e:
            
            # non-transient / unexpected error -> log and return empty dict
            print(f"Error fetching records from dbSNP: {e}")
            if handle is not None:
                handle.close()
            return {}

    return records


# Function to extract effect and non-effect alleles from NCBI dbSNP record
def extract_alleles_from_record(record: dict) -> tuple:
    """
    Extract the effect allele and non-effect allele from the NCBI dbSNP record.

    Args:
        record (dict): The NCBI dbSNP record.

    Returns:
        tuple: A tuple containing the effect allele and non-effect allele, if found, otherwise ("NA", "NA").
    """
    
    # Extract DOCSUM field from the record
    docsum = record.get('DOCSUM', None)
    if not docsum:
        print(f"Error: DOCSUM not found in record for rsID rs{record.get('@uid', 'Unknown')}")
        return "NA", "NA"

    # Convert DOCSUM string to a dictionary for easier access
    docsum_dict = {k.strip(): v.strip() for k, v in (s.split("=", 1) for s in docsum.split('|'))}

    # Extract the SEQ field from the DOCSUM dictionary
    ref_alt = docsum_dict.get('SEQ', None)
    if not ref_alt:
        print(f"Error: SEQ not found in DOCSUM for rsID rs{record.get('@uid', 'Unknown')}")
        return "NA", "NA"
    
    # Convert the SEQ field to a list of alleles and extract the reference and alternate alleles
    ref_alt_list = ref_alt.split('/')
    ref = ref_alt_list[0].replace('[', '')
    alt = ref_alt_list[1].replace(']', '')

    return ref, alt


# Function to extract minor allele frequency (MAF) from NCBI dbSNP record
def extract_maf_from_record(record: dict) -> float:
    """
    Extract the minor allele frequency (MAF) from the NCBI dbSNP record.
    """

    # Extract GLOBAL_MAFS field from the record
    global_mafs = record.get('GLOBAL_MAFS', None)
    if not global_mafs:
        print(f"Error: GLOBAL_MAFS not found in record for rsID rs{record.get('@uid', 'Unknown')}")
        return np.nan

    # Extract MAFs list from the GLOBAL_MAFS field
    mafs_list = global_mafs.get('MAF', None)
    if not mafs_list:
        print(f"Error: MAFs not found in GLOBAL_MAFS for rsID rs{record.get('@uid', 'Unknown')}")
        return np.nan

    # Convert the MAFs list to a dictionary for easier access
    try:
        mafs_dict = {x['STUDY']: x['FREQ'] for x in mafs_list}
    except TypeError:
        
        # Records with single MAF entry may be a dictionary instead of a list, so handle that case
        try:
            mafs_dict = {mafs_list['STUDY']: mafs_list['FREQ']}

        # Handle cases where the MAFs list is not in the expected format
        except Exception as e:
            print(f"Error: MAFs field not in expected format for rsID rs{record.get('@uid', 'Unknown')}: {e}")
            return np.nan

    # Extract MAF for specific studies, prioritizing GnomAD_genomes, then 1000Genomes, and finally any available study
    maf = mafs_dict.get('GnomAD_genomes') or mafs_dict.get('1000Genomes') or list(mafs_dict.values())[0]

    # Parse MAF string to extract the frequency of the alternate allele
    maf_parsed = maf.split("=")
    freqs = maf_parsed[1].split("/")

    return float(freqs[0])

# Function to add dbSNP information (alleles and MAF) to the DataFrame
def add_dbSNP_info(df: pd.DataFrame, ref_col: str, alt_col: str, maf_col: str) -> pd.DataFrame:
    """
    Query NCBI dbSNP to retrieve allele information and MAF based on rsID.

    Args:
        df (pd.DataFrame): The input DataFrame with rsID as the index.
        ref_col (str): The name of the column to store the reference allele.
        alt_col (str): The name of the column to store the alternate allele.
        maf_col (str): The name of the column to store the minor allele frequency (MAF).

    Returns:
        pd.DataFrame: The input DataFrame with additional columns for allele information and MAF.

    """

    rsids_full = df.index.tolist()
    batch_size = 200
    alleles_buf = []
    mafs_buf = []
    for start in range(0, len(rsids_full), batch_size):
        end = min(start + batch_size, len(rsids_full))
        batch_rsids = rsids_full[start:end]

        # Fetch the dbSNP record for the given rsID
        batch_records = fetch_rsid_records(batch_rsids)

        # Extract rsid from records for indexing
        indices = [f"rs{record.get('@uid')}" for record in batch_records]

        # Extract the effect and non-effect alleles from the records and add them to buffer
        batch_alleles = pd.Series([extract_alleles_from_record(record) for record in batch_records], index=indices)
        alleles_buf.append(batch_alleles)

        # Extract the minor allele frequencies (MAFs) from the records
        batch_mafs = pd.Series([extract_maf_from_record(record) for record in batch_records], index=indices)
        mafs_buf.append(batch_mafs)

        # Update the start index for the next batch
        start = end

        # Clean up memory
        del batch_records, batch_alleles, batch_mafs
        gc.collect()

    
    # Concatenate the allele and MAF buffers into single Series
    alleles = pd.concat(alleles_buf)
    mafs = pd.concat(mafs_buf)

    # Split the alleles into separate columns for reference and alternate alleles
    alleles_expanded = pd.DataFrame(
        data=list(alleles),              # list of 2-tuples -> 2 columns
        index=alleles.index,
        columns=[ref_col, alt_col]
    )

    # Add the alleles and MAFs to the DataFrame
    df.loc[mafs.index, maf_col] = mafs
    df = df.join(alleles_expanded, how='left')

    return df


# Function to check if the alleles from the summary statistics match the alleles from dbSNP
def check_allele_match(
        df: pd.DataFrame,
        dbSNP_ref_col: str,
        dbSNP_alt_col: str,
        non_effect_allele_col: str,
        effect_allele_col: str,
        allele_match_col: str
    ) -> pd.DataFrame:
    """
    Check if the alleles from the summary statistics match the alleles from dbSNP.

    Args:
        df (pd.DataFrame): The input DataFrame.
        dbSNP_ref_col (str): The reference allele col from dbSNP.
        dbSNP_alt_col (str): The alternate allele col from dbSNP.
        effect_allele_col (str): The effect allele col from the summary statistics.
        non_effect_allele_col (str): The non-effect allele col from the summary statistics.
        allele_match_col (str): The allele match col from the summary statistics.

    Returns:
        df (pd.DataFrame): The DataFrame with an additional column indicating whether the alleles match.

    """

    # Check if the alleles from the summary statistics match the alleles from dbSNP
    df[allele_match_col] = np.where(
        df[[dbSNP_ref_col, dbSNP_alt_col, non_effect_allele_col, effect_allele_col]].notna().all(axis=1), #type: ignore (silences pylance warning)
        (df[dbSNP_ref_col] == df[non_effect_allele_col]) & (df[dbSNP_alt_col] == df[effect_allele_col]),
        "NA"
    )

    return df


# Function to create a marker ID in the format chr:bp:ref:alt for each row in the DataFrame
def add_variant_id(
        df: pd.DataFrame,
        chr_col: str,
        pos_col: str, 
        ref_allele_col: str,
        alt_allele_col: str,
        variant_id_col: str
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

    # Create the marker ID using vectorized operations for efficiency
    df[variant_id_col] = np.where(
        df[[ref_allele_col, alt_allele_col]].notna().all(axis=1), #type: ignore (silences pylance warning)
        df[chr_col].astype(str) + ':' + df[pos_col].astype(str) + ':' +
        df[ref_allele_col].astype(str) + ':' + df[alt_allele_col].astype(str),
        "NA")
    
    return df


# Function to calculate missing summary statistics (N, Statistic, SE, Beta, Varbeta, P) if not provided
def calculate_missing_summary_statistics(
        df: pd.DataFrame,
        n_col: Optional[str],
        maf_col: Optional[str],
        mac_col: Optional[str],
        p_col: Optional[str],
        statistic_col: Optional[str],
        se_col: Optional[str],
        beta_col: Optional[str],
        var_beta_col: Optional[str],
        sdY_col: Optional[str]
    ) -> pd.DataFrame:
    '''
    Calculate missing summary statistics (N, Statistic, SE, Beta, Varbeta, P) if not provided.

    Args:
        df (pd.DataFrame): The input DataFrame.
        n_col (Optional[str]): The column name for N (sample size).
        maf_col (Optional[str]): The column name for MAF (minor allele frequency).
        mac_col (Optional[str]): The column name for MAC (minor allele count).
        p_col (Optional[str]): The column name for P (p-value).
        statistic_col (Optional[str]): The column name for the test statistic.
        se_col (Optional[str]): The column name for SE (standard error).
        beta_col (Optional[str]): The column name for Beta (effect size).
        var_beta_col (Optional[str]): The column name for VarBeta (variance of effect size).
        sdY_col (Optional[str]): The column name for SDY (standard deviation of Y).

    Returns:
        pd.DataFrame: The DataFrame with missing summary statistics filled in.

    '''

    # Convert relevant columns to numeric, coercing errors to NaN
    s_n = _numeric_series(df, n_col)
    s_maf = _numeric_series(df, maf_col)
    s_mac = _numeric_series(df, mac_col)
    s_p = _numeric_series(df, p_col)
    s_statistic = _numeric_series(df, statistic_col)
    s_se = _numeric_series(df, se_col)
    s_beta = _numeric_series(df, beta_col)
    s_var_beta = _numeric_series(df, var_beta_col)

    # N (sample size)
    if n_col is None:

        # Set default N column name and initialize N series with NaN values
        n_col = 'N'
        s_n = pd.Series(np.nan, index=df.index)

        # Calculate N from MAC and MAF
        m = s_mac.notna() & s_maf.notna() & (s_maf != 0)
        s_n.loc[m] = (s_mac / (2 * s_maf))[m]

        # Add the calculated N to the DataFrame
        df[n_col] = s_n

    # Statistic (used for calculating SE and Beta)
    if statistic_col is None:

        # Set default statistic column name and initialize statistic series with NaN values
        statistic_col = 'STATISTIC'
        se_statistic = pd.Series(np.nan, index=s_p.index)

        # Calculate statistic from p-value using the inverse of the normal cumulative distribution function
        m = s_p.notna() & (s_p > 0) & (s_p < 1)
        se_statistic.loc[m] = (norm.ppf(1 - s_p / 2))[m]
        df[statistic_col] = s_statistic

    # SE
    if not se_col:

        # Set default SE column name and initialize SE series with NaN values
        se_col = 'SE'
        s_se = pd.Series(np.nan, index=df.index)

        # Calculate SE from beta and statistic
        m = s_beta.notna() & s_statistic.notna() & (s_statistic != 0)
        s_se.loc[m] = (s_beta.abs() / s_statistic)[m]

        # Calculate SE from MAF, N, and statistic 
        m = s_se.isna() & s_maf.notna() & s_n.notna() & (s_n > 0) & s_statistic.notna() & (s_statistic != 0)
        s_se.loc[m] = 1 / np.sqrt(2 * s_maf * (1 - s_maf) * (s_n + s_statistic**2))[m]

        # Add the calculated SE to the DataFrame
        df[se_col] = s_se

    # Beta
    if not beta_col:

        # Set default beta column name and initialize beta series with NaN values
        beta_col = 'BETA'
        s_beta = pd.Series(np.nan, index=df.index)

        # statistic * se
        m = s_se.notna() & s_statistic.notna()
        s_beta.loc[m] = (s_statistic * s_se)[m]

        # Add the calculated beta to the DataFrame
        df[beta_col] = s_beta

    # var(beta)
    if not var_beta_col:

        # Set default var(beta) column name and initialize var(beta) series with NaN values
        var_beta_col = 'VARBETA'
        s_var_beta = pd.Series(np.nan, index=df.index)

        # Calculate var(beta) from SE
        m = s_se.notna() & (s_se != 0)
        s_var_beta.loc[m] = (s_se ** 2)[m]

        # Add the calculated var(beta) to the DataFrame
        df[var_beta_col] = s_var_beta

    # sdY
    if not sdY_col:
        sdY_col = 'SDY'
        sdY = pd.Series(np.nan, index=df.index)

        # Calculate sdY from MAF and N
        m = s_maf.notna() & s_n.notna() & (s_n > 0)
        sdY.loc[m] = np.sqrt(s_var_beta * s_n * 2 * s_maf * (1 - s_maf))[m]

        df[sdY_col] = sdY

    # P-value
    if not p_col:
        
        # Set default p-value column name and initialize p-value series with NaN values
        p_col = 'P'
        s_p = pd.Series(np.nan, index=df.index)
        
        # Calculate p-value from statistic using the survival function of the normal distribution
        m = s_p.isna() & s_statistic.notna()
        s_p.loc[m] = 2 * norm.sf(s_statistic.abs())[m]

        # Add the calculated p-value to the DataFrame
        df[p_col] = s_p

    return df


# Function to annotate summary statistics DataFrame with additional information
def annotate_sum_stats(
        df: pd.DataFrame,
        chr_col: str,
        pos_col: str,  
        rsid_col: Optional[str],
        non_effect_allele_col: Optional[str],
        effect_allele_col: Optional[str],
        allele_match_col: str,
        variant_id_col: str,
        n_col: Optional[str],
        maf_col: Optional[str],
        mac_col: Optional[str],
        p_col: Optional[str],
        statistic_col: Optional[str],
        se_col: Optional[str],
        beta_col: Optional[str],
        var_beta_col: Optional[str],
        sdY_col: Optional[str]
    ) -> pd.DataFrame:
    """
    Transform the coloc results DataFrame by selecting relevant columns and renaming them.

    Args:
        df (pd.DataFrame): The input DataFrame containing coloc results.
        chr_col (str): The name of the column containing chromosome information.
        pos_col (str): The name of the column containing position information.
        rsid_col (Optional[str]): The name of the column containing rsID information.
        non_effect_allele_col (Optional[str]): The name of the column containing non effect allele information.
        effect_allele_col (Optional[str]): The name of the column containing effect allele information.
        allele_match_col (str): The name of the column indicating whether alleles match.
        variant_id_col (str): The name of the column containing variant IDs.
        pval_col (Optional[str]): The name of the column containing p-value information.
        beta_col (Optional[str]): The name of the column containing beta coefficient information.
        se_col (Optional[str]): The name of the column containing standard error information.
        n_col (Optional[str]): The name of the column containing sample size information.
        t_col (Optional[str]): The name of the column containing t-statistic information.
        z_col (Optional[str]): The name of the column containing z-score information.
        maf_col (Optional[str]): The name of the column containing minor allele frequency information.
        mac_col (Optional[str]): The name of the column containing minor allele count information.

    Returns:
        pd.DataFrame: The transformed DataFrame with selected and renamed columns.

    """

    # Subset the DataFrame to unique SNPs based on chromosome and position, and fetch rsID if not provided
    if rsid_col is None:
        rsid_col = 'rsID'
        snps_df = df[[chr_col, pos_col]].drop_duplicates()
        snps_df[rsid_col] = snps_df.apply(lambda row: fetch_rsid(row[chr_col], row[pos_col]), axis=1)
        snps_df = snps_df[snps_df[rsid_col] != 'NA']  # Filter out rows where rsID is NA
        snps_df.set_index(rsid_col, inplace=True)
    
    else:
        snps_df = df[[chr_col, pos_col, rsid_col]].drop_duplicates()
        snps_df.set_index(rsid_col, inplace=True)

    # Add dbSNP information (ref allele, alt allele, MAF) based on rsID, skip if rsID is NULL/empty
    dbSNP_ref_col = 'ref_dbSNP'
    dbSNP_alt_col = 'alt_dbSNP'
    dbSNP_maf_col = 'MAF'
    if len(snps_df) > 0:
        snps_df = pd.DataFrame(snps_df)  # Ensure snps_df is a DataFrame
        snps_df = add_dbSNP_info(df = snps_df, ref_col = dbSNP_ref_col, alt_col = dbSNP_alt_col, maf_col = dbSNP_maf_col)

    # Merge the dbSNP information back into the original DataFrame
    if rsid_col not in df.columns:
        annotated_df = pd.merge(df, snps_df, on=[chr_col, pos_col], how='left')
    else:
        annotated_df = pd.merge(df, snps_df, on=[chr_col, pos_col, rsid_col], how='left')

    # Use MAF from dbSNP if maf_col is not provided
    if maf_col is None and dbSNP_maf_col in annotated_df.columns:
        maf_col = dbSNP_maf_col

    # Check provided effect/non-effect alleles from sum stats against dbSNP alleles, if available. If not available, set to NA.
    if non_effect_allele_col is not None and effect_allele_col is not None:
        annotated_df_match_checked = check_allele_match(
            df = annotated_df,
            dbSNP_ref_col = dbSNP_ref_col,
            dbSNP_alt_col = dbSNP_alt_col,
            non_effect_allele_col = non_effect_allele_col,
            effect_allele_col = effect_allele_col,
            allele_match_col = allele_match_col
        )

    # Add variant ID column in the format chr:bp:ref:alt. Use ref/alt from dbSNP if available, otherwise use non-effect/effect alleles if available, else set to NA.
    if allele_match_col in annotated_df.columns:
        annotated_df_var_id = add_variant_id(
            df = annotated_df_match_checked,
            chr_col = chr_col,
            pos_col = pos_col,
            ref_allele_col = dbSNP_ref_col,
            alt_allele_col = dbSNP_alt_col,
            variant_id_col = variant_id_col,
        )
    else:
        annotated_df_var_id = add_variant_id(
            df = annotated_df,
            chr_col = chr_col,
            pos_col = pos_col,
            ref_allele_col = dbSNP_ref_col,
            alt_allele_col = dbSNP_alt_col,
            variant_id_col = variant_id_col,
        )

    # Calculate missing summary statistics (N, SE, Beta, P) if not provided
    annotated_df_full = calculate_missing_summary_statistics(
        df = annotated_df_var_id,
        n_col = n_col,
        maf_col = maf_col,
        mac_col = mac_col,
        p_col = p_col,
        statistic_col = statistic_col,
        se_col = se_col,
        beta_col = beta_col,
        var_beta_col= var_beta_col,
        sdY_col= sdY_col
    )

    return annotated_df_full


# Function to transform the summary statistics DataFrame by selecting relevant columns and renaming them
def transform_sum_stats(
        df: pd.DataFrame,
        chr_col: str,
        pos_col: str,
        variant_id_col: str,
        rsid_col: str,
        gene_id_col: str,
        non_effect_allele_col: str,
        effect_allele_col: str,
        p_col: str,
        beta_col: str,
        var_beta_col: str,
        sdY_col: str,
        se_col: str,
        maf_col: str,
        n_col: str,
        tissue_col: str
    ) -> pd.DataFrame:
    """
    Transform the coloc results DataFrame by selecting relevant columns and renaming them.

    Args:
        df (pd.DataFrame): The input DataFrame containing coloc results.
        chr_col (str): The name of the column containing chromosome information.
        pos_col (str): The name of the column containing position information.
        rsid_col (str): The name of the column containing rsID information.
        gene_id_col (str): The name of the column containing gene ID information.
        non_effect_allele_col (str): The name of the column containing non effect allele information.
        effect_allele_col (str): The name of the column containing effect allele information.
        pval_col (str): The name of the column containing p-value information.
        beta_col (str): The name of the column containing beta coefficient information.
        var_beta_col (str): The name of the column containing variance of beta information.
        sdY_col (str): The name of the column containing standard deviation of Y information.
        se_col (str): The name of the column containing standard error information.
        maf_col (str): The name of the column containing minor allele frequency information.
        n_col (str): The name of the column containing sample size information.
        tissue_col (str): The name of the column containing tissue information.

    Returns:
        pd.DataFrame: The transformed DataFrame with selected and renamed columns.

    """

    # Ensure that the tissue column exists in the DataFrame, if not, create it with default value 'NA'
    if not tissue_col:
        tissue_col = 'TISSUE'
        df[tissue_col] = 'NA'
    
    # Select relevant columns
    transformed_df = pd.DataFrame(df[[chr_col, pos_col, variant_id_col, rsid_col,
                         gene_id_col, non_effect_allele_col, effect_allele_col,
                         p_col, beta_col, var_beta_col, sdY_col, se_col, maf_col,
                         n_col, tissue_col]])

    # Rename columns
    transformed_df.columns = ['CHR', 'BP', 'VARIANT_ID', 'RSID', 
                              'GENE_ID', 'NON_EFFECT_ALLELE', 'EFFECT_ALLELE', 
                              'P', 'BETA', 'VAR_BETA', 'SDY', 'SE', 'MAF', 'N', 'TISSUE']

    return transformed_df


# Function to save the standardized DataFrame and the non-standardized (NA, mismatched alleles) DataFrame to output files
def save_output(annotated_df: pd.DataFrame, annotated_df_mismatched: pd.DataFrame, transformed_df: pd.DataFrame, output_dir: str, input_fp: Path):
    """
    Save the standardized DataFrame and the non-standardized (NA, mismatched alleles) DataFrame to output files.

    Args:
        annotated_df (pd.DataFrame): The annotated DataFrame containing summary statistics.
        annotated_df_mismatched (pd.DataFrame): The DataFrame containing rows with mismatched alleles or NA values.
        transformed_df (pd.DataFrame): The transformed DataFrame containing standardized summary statistics.
        output_dir (str): The directory to save the output files.
        input_fp (Path): The input file path used to derive the output file names.

    Returns:
        None
    """

    # Filter out rows where dbSNP annotation failed (VARIANT_ID is NA) to create the standardized DataFrame
    standardized_df = transformed_df[transformed_df['VARIANT_ID'] != 'NA']

    # Create a DataFrame for rows where dbSNP annotation failed (VARIANT_ID is NA) for retry
    annotated_df_failed_dbsnp = annotated_df[annotated_df.index.isin(transformed_df.index) & ~annotated_df.index.isin(standardized_df.index)]

    # Extract the base name of the input file for naming the output files
    if input_fp.suffix == '.gz':
        input_file_stem = input_fp.stem  # Remove .gz
        input_file_base = Path(input_file_stem).stem  # Remove file extension (e.g., .csv, .tsv, .txt)
    else:
        input_file_base = input_fp.stem  # Remove file extension (e.g., .csv, .tsv, .txt)
    
    # Save the transformed DataFrame to the output file
    if len(standardized_df) > 0:
        standardized_fp = Path(f"{output_dir}/{input_file_base}.standardized.tsv")
        standardized_df.to_csv(standardized_fp, sep='\t', index=False)

    # Save the mismatched alleles DataFrame to a separate file
    if len(annotated_df_failed_dbsnp) > 0:
        failed_dbsnp_fp = Path(f"{output_dir}/{input_file_base}.failed_dbsnp.tsv")
        annotated_df_failed_dbsnp.to_csv(failed_dbsnp_fp, sep='\t', index=False)

    # Save the mismatched alleles DataFrame to a separate file
    if len(annotated_df_mismatched) > 0:
        mismatched_alleles_fp = Path(f"{output_dir}/{input_file_base}.mismatched_alleles.tsv")
        annotated_df_mismatched.to_csv(mismatched_alleles_fp, sep='\t', index=False)


# ------------------------- Full Pipeline -------------------------

# Main function to execute the full pipeline
def main(args: argparse.Namespace):
    """
    Main function to load input data, annotate summary statistics, transform summary statistics, and save the transformed data.

    Args:
        args (argparse.Namespace): Command line arguments.
    """

    # Load input data
    input_fp = Path(args.input_file)
    df = load_input_data(input_fp, args.header_lines)
    print(f"Loaded input data from {input_fp} with shape: {df.shape}")

    # Load expanded ranges from file
    expanded_ranges_fp = Path(args.expanded_ranges_file)
    expanded_ranges = load_expanded_ranges(expanded_ranges_fp)

    # Normalize chromosome column
    df[args.chr_col] = df[args.chr_col].astype(str).map(normalize_chromosome)

    # Filter the DataFrame by expanded ranges
    filtered_df = filter_by_expanded_ranges(
        df=df,
        chr_col=args.chr_col,
        pos_col=args.pos_col,
        expanded_ranges=expanded_ranges
    )

    # Initialize Entrez email (required by NCBI) for API access.
    Entrez.email = args.entrez_email

    # Annotate summary statistics
    allele_match_col = 'Allele_Match'
    variant_id_col = 'Variant_ID'
    annotated_df = annotate_sum_stats(
        df=filtered_df,
        chr_col=args.chr_col,
        pos_col=args.pos_col,
        rsid_col=args.rsid_col,
        non_effect_allele_col=args.non_effect_allele_col,
        effect_allele_col=args.effect_allele_col,
        allele_match_col=allele_match_col,
        variant_id_col=variant_id_col,
        n_col=args.n_col,
        maf_col=args.maf_col,
        mac_col=args.mac_col,
        p_col=args.p_col,
        statistic_col=args.statistic_col,
        se_col=args.se_col,
        beta_col=args.beta_col,
        var_beta_col=args.var_beta_col,
        sdY_col=args.sdY_col
    )

    # Filter out rows with mismatched alleles (Allele_Match == False)
    if allele_match_col in annotated_df.columns:
        annotated_df_mismatched = annotated_df[annotated_df[allele_match_col] == False]
        annotated_df_mismatched = pd.DataFrame(annotated_df_mismatched)  # Ensure annotated_df_mismatched is a DataFrame
        annotated_df = annotated_df[~annotated_df.index.isin(annotated_df_mismatched.index)]
    else:
        annotated_df_mismatched = pd.DataFrame()  # Empty DataFrame if Allele_Match column does not exist

    # Transform summary statistics
    transformed_df = transform_sum_stats(
        df=annotated_df,
        chr_col=args.chr_col,
        pos_col=args.pos_col,
        variant_id_col=variant_id_col,
        rsid_col=args.rsid_col or 'rsID',
        gene_id_col=args.gene_id_col,
        effect_allele_col=args.effect_allele_col or 'dbSNP_alt',
        non_effect_allele_col=args.non_effect_allele_col or 'dbSNP_ref',
        p_col=args.p_col or 'P',
        beta_col=args.beta_col or 'BETA',
        se_col=args.se_col or 'SE',
        n_col=args.n_col or 'N',
        maf_col=args.maf_col or 'MAF',
        var_beta_col=args.var_beta_col or 'VARBETA',
        sdY_col=args.sdY_col or 'SDY',
        tissue_col=args.tissue_col
    )

    # Save the output files
    save_output(
        annotated_df=annotated_df,
        annotated_df_mismatched=annotated_df_mismatched,
        transformed_df=transformed_df,
        output_dir=args.output_dir,
        input_fp=input_fp
    )


# ------------------------- Command Line Interface -------------------------

if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Transform coloc results to standardized format.")
    parser.add_argument("--input_file", type=str, help="Path to the input file containing coloc results.")
    parser.add_argument("--expanded_ranges_file", type=str, help="Path to the file containing expanded ranges for filtering.")
    parser.add_argument("--output_dir", type=str, help="Directory to save the transformed output files.")
    parser.add_argument("--header_lines", type=int, default=0, help="Number of header lines to skip in the input file.")
    parser.add_argument("--chr_col", type=str, required=True, help="Column name for chromosome information.")
    parser.add_argument("--pos_col", type=str, required=True, help="Column name for position information.")
    parser.add_argument("--rsid_col", type=str, default=None, help="Column name for rsID information (optional).")
    parser.add_argument("--entrez_email", type=str, required=True, help="Email address for NCBI Entrez API.")
    parser.add_argument("--gene_id_col", type=str, required=True, help="Column name for gene ID information.")
    parser.add_argument("--non_effect_allele_col", type=str, required=True, help="Column name for non-effect allele information.")
    parser.add_argument("--effect_allele_col", type=str, required=True, help="Column name for effect allele information.")
    parser.add_argument("--p_col", type=str, default=None, help="Column name for p-value information (optional).")
    parser.add_argument("--beta_col", type=str, default=None, help="Column name for beta coefficient information (optional).")
    parser.add_argument("--se_col", type=str, default=None, help="Column name for standard error information (optional).")
    parser.add_argument("--n_col", type=str, default=None, help="Column name for sample size information (optional).")
    parser.add_argument("--statistic_col", type=str, default=None, help="Column name for statistic (z/t) information (optional).")
    parser.add_argument("--maf_col", type=str, default=None, help="Column name for minor allele frequency information (optional).")
    parser.add_argument("--mac_col", type=str, default=None, help="Column name for minor allele count information (optional).")
    parser.add_argument("--var_beta_col", type=str, default=None, help="Column name for variance of beta information (optional).")
    parser.add_argument("--sdY_col", type=str, default=None, help="Column name for standard deviation of Y information (optional).")
    parser.add_argument("--tissue_col", type=str, default=None, help="Column name for tissue information (optional).")

    # Parse arguments
    args = parser.parse_args()

    # Call main function with parsed arguments
    main(args)