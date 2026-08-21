# Imports
import gc
import time
import xmltodict
import argparse
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
from helpers import load_input_data, normalize_chromosome, numeric_series


# ------------------------- Helper Function Definitions -------------------------

@lru_cache(maxsize=65536)
def _fetch_rsid(chromosome: str, position: int) -> str:
    """
    Query NCBI dbSNP to retrieve rsID information based on chromosome and position.

    Args:
        chromosome (str): The chromosome number or identifier.
        position (int): The base pair position on the chromosome.

    Returns:
        str: The corresponding rsID if found, otherwise "NA".

    """
    # Construct the search term for the NCBI dbSNP query
    search_term = f"{chromosome}[Chromosome] AND {position}[Base Position] AND Homo sapiens[Organism]"
    
    # Query dbSNP via Esearch for search term
    handle = Entrez.esearch(db="snp", term=search_term)

    # Parse the response
    response = xmltodict.parse(handle) # type: ignore (silences pylance warning)
    handle.close()

    # Extract the rsIDs from the response
    rs_ids = response.get("IdList", None)

    # Check if any rsIDs were found. If not found return "NA", if found return the last one
    if not rs_ids:
        return "NA"
    return f"rs{str(rs_ids[-1]).strip()}"


def _fetch_rsid_records(rsids: list) -> dict:
    """
    Query NCBI dbSNP to retrieve allele information based on a list of rsIDs.

    Args:
        rsids (list): The list of rsIDs to fetch.

    Returns:
        dict: Dictionary containing the fetched rsID records.

    """
    # Join the rsIDs into a comma-separated string
    rsids_str = ",".join(rsids)
    
    # Attempt dbSNP Query, retry upon exception
    attempts = 0
    while attempts < 5:
        attempts += 1
        handle = None
        try:
            # Fetch the records from dbSNP using Entrez.efetch
            handle = Entrez.efetch(db="snp", id=rsids_str, retmode="xml")

            # Parse the XML response into a dictionary
            response = xmltodict.parse(handle) # type: ignore (silences pylance warning)
            handle.close()

            # Extract the records from the response dictionary
            records = response.get("ExchangeSet", {}).get("DocumentSummary", {})

        # Handle specific exceptions that may result from transient errors during the fetch operation by triggering retry logic
        except (IncompleteRead, HTTPError, URLError, ExpatError) as e:
            if handle is not None:
                handle.close()
            if attempts < 5:
                time.sleep(2 * attempts)
                continue
            print(f"_fetch_rsid_records: reached max retries, returning empty result for rsids: {rsids[0]} - {rsids[-1]}")
            return {}
        
        # Handle any other unexpected exceptions that may occur during the fetch operation by logging the error and returning an empty result
        except Exception as e:
            print(f"Error fetching records from dbSNP: {e}")
            if handle is not None:
                handle.close()
            return {"Failure": {}}

    return records

def _extract_position_from_record(record: dict) -> int:
    """
    Extract the genomic position from an NCBI dbSNP record.
    
    Args:
        record (dict): The NCBI dbSNP record.

    Returns:
        int: The genomic position, or NaN if not found.
    """
    # Extract the CHRPOS field from the record, which contains the genomic position
    position = record.get("CHRPOS", None)

    # Check if CHRPOS is present; if not, log an error and return NaN
    if position is None:
        print(f"Error: CHRPOS not found in record for rsID rs{record.get('@uid', 'Unknown')}")
        return np.nan
    
    # Extract the position. Pre-colon part is the chromosome, post-colon part is the position.
    position = position.split(":")[-1]  

    return int(position)


def _extract_alleles_from_record(record: dict) -> tuple:
    """
    Extract the reference and alternate alleles from an NCBI dbSNP record.

    Args:
        record (dict): The NCBI dbSNP record.

    Returns:
        tuple: A tuple of (ref_allele, alt_allele), or ("NA", "NA") if not found.

    """
    # Extract the DOCSUM field from the record, which contains allele information
    docsum = record.get("DOCSUM", None)
    
    # Check if DOCSUM is present; if not, log an error and return ("NA", "NA")
    if not docsum:
        print(f"Error: DOCSUM not found in record for rsID rs{record.get('@uid', 'Unknown')}")
        return "NA", "NA"

    # Parse the DOCSUM string into a dictionary
    docsum_dict = {k.strip(): v.strip() for k, v in (s.split("=", 1) for s in docsum.split("|"))}
    ref_alt = docsum_dict.get("SEQ", None)
    
    # Check if SEQ is present; if not, log an error and return ("NA", "NA")
    if not ref_alt:
        print(f"Error: SEQ not found in DOCSUM for rsID rs{record.get('@uid', 'Unknown')}")
        return "NA", "NA"
    
    # Extract the reference and alternate alleles from the SEQ string, which is formatted as "[ref]/[alt]"
    ref_alt_list = ref_alt.split("/")
    ref = ref_alt_list[0].replace("[", "")
    alt = ref_alt_list[1].replace("]", "")

    return ref, alt


def _extract_maf_from_record(record: dict) -> float:
    """
    Extract the minor allele frequency (MAF) from an NCBI dbSNP record.

    Args:
        record (dict): The NCBI dbSNP record.

    Returns:
        float: The MAF value, or NaN if not found.

    """
    # Extract the GLOBAL_MAFS field from the record, which contains MAF information
    global_mafs = record.get("GLOBAL_MAFS", None)

    # Check if GLOBAL_MAFS is present; if not, log an error and return NaN
    if not global_mafs:
        print(f"Error: GLOBAL_MAFS not found in record for rsID rs{record.get('@uid', 'Unknown')}")
        return np.nan

    # Extract the MAFs list from the GLOBAL_MAFS field
    mafs_list = global_mafs.get("MAF", None)

    # Check if MAFs list is present; if not, log an error and return NaN
    if not mafs_list:
        print(f"Error: MAFs not found in GLOBAL_MAFS for rsID rs{record.get('@uid', 'Unknown')}")
        return np.nan
    
    # Attempt to convert the MAFs list into a dictionary mapping study names to frequencies
    try:
        mafs_dict = {x["STUDY"]: x["FREQ"] for x in mafs_list}

    # Handle cases where the MAFs list is not in the expected format by attempting to extract a single MAF entry
    except TypeError:
        try:
            mafs_dict = {mafs_list["STUDY"]: mafs_list["FREQ"]}
        
        # Handle any other unexpected exceptions that may occur during the MAF extraction process by logging the error and returning NaN
        except Exception as e:
            print(f"Error: MAFs field not in expected format for rsID rs{record.get('@uid', 'Unknown')}: {e}")
            return np.nan
        
    # Select the MAF from GnomAD_genomes if available, otherwise from 1000Genomes, or fallback to the first available MAF
    maf = mafs_dict.get("GnomAD_genomes") or mafs_dict.get("1000Genomes") or list(mafs_dict.values())[0]

    # Split the MAF string to extract the frequency value and convert it to a float
    maf_parsed = maf.split("=")
    freqs = maf_parsed[1].split("/")

    return float(freqs[0])


# ------------------------- Class Definition -------------------------

class SumStatsTransformer:
    """
    A class to transform and annotate summary statistics data with dbSNP information, allele matching, and derived summary statistics.

    """

    def __init__(self, args: argparse.Namespace):
        self.output_dir = Path(args.output_dir)
        self.sum_stats_fp = Path(args.sum_stats_file)
        print(f"Loading input data from {self.sum_stats_fp}...")
        self.ss_df = load_input_data(self.sum_stats_fp, args.header_lines)
        print(f"Loaded input data from {self.sum_stats_fp} with shape: {self.ss_df.shape}")
        self.loci_fp = Path(self.output_dir, args.loci_file)
        print(f"Loading input data from {self.loci_fp}...")
        self.loci_df = load_input_data(self.loci_fp, 0)
        print(f"Loaded input data from {self.loci_fp} with shape: {self.loci_df.shape}")
        self.standardized_chr_key = args.standardized_chr_key
        self.loci_left_bound_key = args.loci_left_bound_key
        self.loci_right_bound_key = args.loci_right_bound_key
        self.snps_df = pd.DataFrame()
        self.ss_chr_key = args.ss_chr_key
        self.ss_pos_key = args.ss_pos_key
        self.ss_rsid_key = args.ss_rsid_key
        self.ss_non_effect_allele_key = args.ss_non_effect_allele_key
        self.ss_effect_allele_key = args.ss_effect_allele_key
        self.ss_n_key = args.ss_n_key
        self.ss_maf_key = args.ss_maf_key
        self.ss_mac_key = args.ss_mac_key
        self.ss_p_key = args.ss_p_key
        self.ss_statistic_key = args.ss_statistic_key
        self.ss_se_key = args.ss_se_key
        self.ss_beta_key = args.ss_beta_key
        self.ss_var_beta_key = args.ss_var_beta_key
        self.ss_sdy_key = args.ss_sdy_key
        self.ss_gene_id_key = args.ss_gene_id_key
        self.ss_tissue_key = args.ss_tissue_key
        self.annotated_df = pd.DataFrame()
        self.dbSNP_rsid_key = "rsID_dbsnp"
        self.dbSNP_pos_key = "POS_dbsnp"
        self.dbSNP_non_effect_allele_key = "A1_dbsnp"
        self.dbSNP_effect_allele_key = "A2_dbsnp"
        self.dbSNP_maf_key = "MAF_dbsnp"
        self.dbSNP_validation_key = "Allele_Match"
        self.standardized_non_effect_allele_key = args.standardized_non_effect_allele_key
        self.standardized_effect_allele_key = args.standardized_effect_allele_key
        self.standardized_maf_key = args.standardized_maf_key
        self.standardized_pos_key = args.standardized_pos_key
        self.standardized_variant_id_key = args.standardized_variant_id_key
        self.standardized_rsid_key = args.standardized_rsid_key
        self.standardized_gene_id_key = args.standardized_gene_id_key
        self.standardized_non_effect_allele_key = args.standardized_non_effect_allele_key
        self.standardized_effect_allele_key = args.standardized_effect_allele_key
        self.standardized_n_key = args.standardized_n_key
        self.standardized_p_key = args.standardized_p_key
        self.standardized_se_key = args.standardized_se_key
        self.standardized_beta_key = args.standardized_beta_key
        self.standardized_var_beta_key = args.standardized_var_beta_key
        self.standardized_sdy_key = args.standardized_sdy_key
        self.standardized_statistic_key = args.standardized_statistic_key
        self.standardized_tissue_key = args.standardized_tissue_key
        

    def filter_by_expanded_ranges(self):
        """
        Filter summary statistics to include only rows that fall within specified expanded locus ranges.

        Returns:
            None: The DataFrame is modified in place to include only rows within the specified ranges.

        """
        # Initialize a boolean mask for filtering
        mask = pd.Series(False, index=self.ss_df.index)

        # Normalize loci CHR to match ss_df standardized chromosome format
        self.loci_df[self.standardized_chr_key] = self.loci_df[self.standardized_chr_key].astype(str).map(normalize_chromosome)

        # Iterate over each locus and update the mask for rows that fall within the specified range
        for _, row in self.loci_df.iterrows():
            idx = self.ss_df.index[self.ss_df[self.standardized_chr_key] == row[self.standardized_chr_key]]
            pos_vals = numeric_series(pd.DataFrame(self.ss_df.loc[idx]), self.ss_pos_key)
            chrom_mask = pd.Series(False, index=idx)
            chrom_mask |= pos_vals.between(row[self.loci_left_bound_key], row[self.loci_right_bound_key])
            mask.loc[idx] = chrom_mask

        # Filter the summary statistics DataFrame based on the mask
        self.ss_df = pd.DataFrame(self.ss_df[mask])


    def add_dbSNP_info(self):
        """
        Query NCBI dbSNP in batches to retrieve allele information and MAF for each rsID in the index.

        Returns:
            None: The DataFrame is modified in place to include dbSNP allele and MAF information.

        """
        # Get the list of rsIDs from the index of the summary statistics DataFrame
        rsids_full = self.snps_df.index.tolist()

        # Define the batch size for querying dbSNP
        batch_size = 200

        # Initialize buffers to store allele and MAF information for each batch
        alleles_buf = []
        mafs_buf = []
        pos_buf = []

        # Loop through the rsIDs in batches and fetch allele and MAF information from dbSNP
        print(f"Fetching dbSNP information for {len(rsids_full)} rsIDs in batches of {batch_size}...")
        for start in range(0, len(rsids_full), batch_size):

            # Determine the end index for the current batch and extract the corresponding rsIDs
            end = min(start + batch_size, len(rsids_full))
            print(f"Processing batch {start // batch_size + 1} of {((len(rsids_full) - 1) // batch_size) + 1} (rsIDs {start} to {end - 1})...")
            batch_rsids = rsids_full[start:end]

            # Fetch the dbSNP records for the current batch of rsIDs
            batch_records = _fetch_rsid_records(batch_rsids)

            # Create a list of indices corresponding to the fetched records, using the '@uid' field from each record
            indices = [f"rs{record.get('@uid')}" for record in batch_records]

            # Extract position information from the fetched records, store it in a pandas Series, and append it to the positions buffer
            batch_positions = pd.Series(
                [_extract_position_from_record(record) for record in batch_records],
                index=indices
            )
            pos_buf.append(batch_positions)

            # Extract allele information from the fetched records, store it in a pandas Series, and append it to the alleles buffer
            batch_alleles = pd.Series(
                [_extract_alleles_from_record(record) for record in batch_records],
                index=indices
            )
            alleles_buf.append(batch_alleles)

            # Extract MAF information from the fetched records, store it in a pandas Series, and append it to the MAFs buffer
            batch_mafs = pd.Series(
                [_extract_maf_from_record(record) for record in batch_records],
                index=indices
            )
            mafs_buf.append(batch_mafs)

            # Explicitly delete the batch records, positions, alleles, and MAFs to free up memory
            del batch_records, batch_positions, batch_alleles, batch_mafs
            gc.collect()

        # Concatenate the allele and MAF buffers into single Series for alleles and MAFs
        positions = pd.concat(pos_buf)
        alleles = pd.concat(alleles_buf)
        mafs = pd.concat([s for s in mafs_buf if not s.empty])

        # Deduplicate by index (rsID) keeping first occurrence
        positions = positions[~positions.index.duplicated(keep="first")]
        mafs = mafs[~mafs.index.duplicated(keep="first")]
        alleles = alleles[~alleles.index.duplicated(keep="first")]

        # Filter to only rsIDs present in snps_df to avoid index alignment issues
        positions = positions[positions.index.isin(self.snps_df.index)]
        mafs = mafs[mafs.index.isin(self.snps_df.index)]
        alleles = alleles[alleles.index.isin(self.snps_df.index)]

        # Create a DataFrame from the alleles Series, expanding the tuples into separate columns for non-effect and effect alleles
        alleles_expanded = pd.DataFrame(
            data=list(alleles),
            index=alleles.index,
            columns=[self.dbSNP_non_effect_allele_key, self.dbSNP_effect_allele_key]
        )

        # Update the summary statistics DataFrame with the fetched MAFs and join it with the expanded alleles DataFrame
        self.snps_df.loc[positions.index, self.dbSNP_pos_key] = positions
        self.snps_df.loc[mafs.index, self.dbSNP_maf_key] = mafs
        self.snps_df = self.snps_df.join(alleles_expanded, how="left")


    def validate_dbSNP(self):
        """
        Check whether the alleles from the summary statistics match the alleles from dbSNP.

        Returns:
            None: The DataFrame is modified in place to include a new column indicating allele match status.

        """
        # Check if the dbSNP position column and summary statistics position column are present in the annotated DataFrame
        dbsnp_pos_present = self.annotated_df[[self.dbSNP_pos_key]].notna().all(axis=1) #type: ignore
        ss_pos_present = self.annotated_df[[self.ss_pos_key]].notna().all(axis=1) #type: ignore

        # Check if the dbSNP allele columns and summary statistics allele columns are present in the annotated DataFrame
        dbsnp_alleles_present = self.annotated_df[[self.dbSNP_non_effect_allele_key, self.dbSNP_effect_allele_key]].notna().all(axis=1) #type: ignore
        ss_alleles_present = self.annotated_df[[self.ss_non_effect_allele_key, self.ss_effect_allele_key]].notna().all(axis=1) #type: ignore

        # Create a new column in the DataFrame to validate summary statistics against dbSNP
        if self.ss_pos_key is not None and self.ss_non_effect_allele_key is not None and self.ss_effect_allele_key is not None:
            self.annotated_df[self.dbSNP_validation_key] = np.where(
                ~dbsnp_pos_present & ~dbsnp_alleles_present,
                "Failed dbSNP lookup",
                np.where(
                    dbsnp_pos_present & ss_pos_present & dbsnp_alleles_present & ss_alleles_present,
                    np.where(
                        (self.annotated_df[self.dbSNP_pos_key] == self.annotated_df[self.ss_pos_key]) &
                        (self.annotated_df[self.dbSNP_non_effect_allele_key] == self.annotated_df[self.ss_non_effect_allele_key]) &
                        (self.annotated_df[self.dbSNP_effect_allele_key] == self.annotated_df[self.ss_effect_allele_key]),
                        "Validated",
                        np.where(
                            self.annotated_df[self.dbSNP_pos_key] != self.annotated_df[self.ss_pos_key],
                            "Position Mismatch",
                            "Allele Mismatch"
                        )
                    ),
                    "Incomplete summary statistics columns for validation"
                )
            )

    def add_variant_id(self):
        """
        Create a variant ID in the format chr:bp:ref:alt for each row in the DataFrame.

        Returns:
            None: The DataFrame is modified in place to include a new column with the variant ID.

        """
        # Check if dbSNP alleles and summary statistics alleles are present for each row in the DataFrame
        dbsnp_alleles_present = self.annotated_df[[self.dbSNP_non_effect_allele_key, self.dbSNP_effect_allele_key]].notna().all(axis=1) #type: ignore
        ss_alleles_present = self.annotated_df[[self.ss_non_effect_allele_key, self.ss_effect_allele_key]].notna().all(axis=1) #type: ignore

        # Check if dbSNP position and summary statistics position are present for each row in the DataFrame
        dbsnp_pos_present = self.annotated_df[self.dbSNP_pos_key].notna().all(axis=1) #type: ignore
        ss_pos_present = self.annotated_df[self.ss_pos_key].notna().all(axis=1) #type: ignore

        # Create a new column in the DataFrame to store the variant ID, which is constructed from the chromosome, position, non-effect allele, and effect allele. Use dbSNP alleles if available, otherwise use the alleles from the summary statistics if available or "NA".
        self.annotated_df[self.standardized_variant_id_key] = np.where(
            ~dbsnp_alleles_present & ~ss_alleles_present & ~dbsnp_pos_present & ~ss_pos_present,
            "NA",
            np.where(
                dbsnp_alleles_present & dbsnp_pos_present,
                self.annotated_df[self.ss_chr_key].astype(str) + ":" + self.annotated_df[self.dbSNP_pos_key].astype(str) + ":" +
                self.annotated_df[self.dbSNP_non_effect_allele_key].astype(str) + ":" + self.annotated_df[self.dbSNP_effect_allele_key].astype(str),
                np.where(
                    dbsnp_pos_present,
                    self.annotated_df[self.ss_chr_key].astype(str) + ":" + self.annotated_df[self.dbSNP_pos_key].astype(str) + ":" +
                    self.annotated_df[self.ss_non_effect_allele_key].astype(str) + ":" + self.annotated_df[self.ss_effect_allele_key].astype(str),
                    np.where(
                        dbsnp_alleles_present,
                        self.annotated_df[self.ss_chr_key].astype(str) + ":" + self.annotated_df[self.ss_pos_key].astype(str) + ":" +
                        self.annotated_df[self.dbSNP_non_effect_allele_key].astype(str) + ":" + self.annotated_df[self.dbSNP_effect_allele_key].astype(str),
                        self.annotated_df[self.ss_chr_key].astype(str) + ":" + self.annotated_df[self.ss_pos_key].astype(str) + ":" +
                        self.annotated_df[self.ss_non_effect_allele_key].astype(str) + ":" + self.annotated_df[self.ss_effect_allele_key].astype(str)
                    )
                )
            )
        )


    def calculate_missing_summary_statistics(self):
        """
        Calculate missing summary statistics (N, statistic, SE, beta, var(beta), SDY, P)
        from available columns, filling in defaults where possible.

        Returns:
            pd.DataFrame: The DataFrame with missing summary statistics filled in.

        """
        # Extract numeric series for relevant summary statistics columns
        s_n = numeric_series(self.annotated_df, self.ss_n_key)
        s_maf = numeric_series(self.annotated_df, self.ss_maf_key)
        s_mac = numeric_series(self.annotated_df, self.ss_mac_key)
        s_p = numeric_series(self.annotated_df, self.ss_p_key)
        s_statistic = numeric_series(self.annotated_df, self.ss_statistic_key)
        s_se = numeric_series(self.annotated_df, self.ss_se_key)
        s_beta = numeric_series(self.annotated_df, self.ss_beta_key)
        s_var_beta = numeric_series(self.annotated_df, self.ss_var_beta_key)

        # Calculate N from MAC and MAF if both are available and MAF is not zero
        if self.ss_n_key is None:
            s_n = pd.Series(np.nan, index=self.annotated_df.index)
            m = s_mac.notna() & s_maf.notna() & (s_maf != 0)
            s_n.loc[m] = (s_mac / (2 * s_maf))[m]
            self.annotated_df[self.standardized_n_key] = s_n

        # Calculate statistic from P-value if P-value is available and within (0, 1)
        if self.ss_statistic_key is None:
            s_statistic = pd.Series(np.nan, index=s_p.index)
            m = s_p.notna() & (s_p > 0) & (s_p < 1)
            s_statistic.loc[m] = (norm.ppf(1 - s_p / 2))[m]
            self.annotated_df[self.standardized_statistic_key] = s_statistic

        # Calculate standardized error (SE) from beta and statistic if both are available and statistic is not zero, or from MAF and N if both are available and N is greater than zero
        if not self.ss_se_key:
            s_se = pd.Series(np.nan, index=self.annotated_df.index)
            m = s_beta.notna() & s_statistic.notna() & (s_statistic != 0)
            s_se.loc[m] = (s_beta.abs() / s_statistic)[m]
            m = s_se.isna() & s_maf.notna() & s_n.notna() & (s_n > 0) & s_statistic.notna() & (s_statistic != 0)
            s_se.loc[m] = 1 / np.sqrt(2 * s_maf * (1 - s_maf) * (s_n + s_statistic**2))[m]
            self.annotated_df[self.standardized_se_key] = s_se

        # Calculate beta from statistic and SE if both are available
        if not self.ss_beta_key:
            s_beta = pd.Series(np.nan, index=self.annotated_df.index)
            m = s_se.notna() & s_statistic.notna()
            s_beta.loc[m] = (s_statistic * s_se)[m]
            self.annotated_df[self.standardized_beta_key] = s_beta

        # Calculate variance of beta from SE if SE is available and not zero
        if not self.ss_var_beta_key:
            s_var_beta = pd.Series(np.nan, index=self.annotated_df.index)
            m = s_se.notna() & (s_se != 0)
            s_var_beta.loc[m] = (s_se ** 2)[m]
            self.annotated_df[self.standardized_var_beta_key] = s_var_beta

        # Calculate SDY from variance of beta, N, and MAF if all are available and N is greater than zero
        if not self.ss_sdy_key:
            sdy = pd.Series(np.nan, index=self.annotated_df.index)
            m = s_maf.notna() & s_n.notna() & (s_n > 0)
            sdy.loc[m] = np.sqrt(s_var_beta * s_n * 2 * s_maf * (1 - s_maf))[m]
            self.annotated_df[self.standardized_sdy_key] = sdy

        # Calculate P-value from statistic if statistic is available and not zero, using a two-tailed test
        if not self.ss_p_key:
            s_p = pd.Series(np.nan, index=self.annotated_df.index)
            m = s_p.isna() & s_statistic.notna()
            s_p.loc[m] = 2 * norm.sf(s_statistic.abs())[m]
            self.annotated_df[self.standardized_p_key] = s_p


    def annotate(self):
        """
        Annotate a summary statistics DataFrame with dbSNP allele information, allele
        match flags, variant IDs, and any missing derived summary statistics.

        Returns:
            None: The DataFrame is modified in place to include the new annotations.

        """
        # Subset to unique SNPs and fetch rsID if not provided
        if self.ss_rsid_key is None:
            self.snps_df = self.ss_df[[self.standardized_chr_key, self.ss_pos_key]].drop_duplicates()
            self.snps_df[self.dbSNP_rsid_key] = self.snps_df.apply(
                lambda row: _fetch_rsid(row[self.standardized_chr_key], row[self.ss_pos_key]), axis=1
            )
            self.snps_df = self.snps_df[self.snps_df[self.dbSNP_rsid_key] != "NA"]
            self.snps_df.set_index(self.dbSNP_rsid_key, inplace=True)

        # If rsID is provided, subset to unique SNPs based on chromosome, position, and rsID
        else:
            self.snps_df = self.ss_df[[self.standardized_chr_key, self.ss_pos_key, self.ss_rsid_key]].drop_duplicates()
            self.snps_df.set_index(self.ss_rsid_key, inplace=True)

        # Add dbSNP allele and MAF information
        if len(self.snps_df) > 0:
            self.add_dbSNP_info()

        # Merge dbSNP info back into the original DataFrame
        merge_keys = [self.standardized_chr_key, self.ss_pos_key] if self.ss_rsid_key not in self.ss_df.columns \
            else [self.standardized_chr_key, self.ss_pos_key, self.ss_rsid_key]
        self.annotated_df = pd.merge(self.ss_df, self.snps_df, on=merge_keys, how="left")

        # Use dbSNP MAF if maf_key was not provided
        if self.ss_maf_key is None and self.dbSNP_maf_key in self.annotated_df.columns:
            self.ss_maf_key = self.dbSNP_maf_key

        # Check allele match if effect/non-effect allele columns are available
        self.validate_dbSNP()

        # Add variant ID 
        self.add_variant_id()

        # Calculate missing summary statistics
        self.calculate_missing_summary_statistics()


    def transform(self):
        """
        Transform the annotated DataFrame into the standardized column layout expected
        by downstream colocalization tools.

        Returns:
            None: The DataFrame is modified in place to include only the standardized columns.

        """
        # Ensure tissue column exists
        if not self.ss_tissue_key:
            self.annotated_df[self.standardized_tissue_key] = "NA"

        # Ensure gene ID column exists
        if not self.ss_gene_id_key:
            self.annotated_df[self.ss_gene_id_key] = "NA"

        # Select columns to match the standardized layout
        self.transformed_df = pd.DataFrame(self.annotated_df[[
            self.standardized_chr_key,
            self.dbSNP_pos_key or self.ss_pos_key,
            self.standardized_variant_id_key,
            self.ss_rsid_key or self.dbSNP_rsid_key,
            self.ss_gene_id_key,
            self.dbSNP_non_effect_allele_key or self.ss_non_effect_allele_key,
            self.dbSNP_effect_allele_key or self.ss_effect_allele_key,
            self.ss_p_key or self.standardized_p_key,
            self.ss_beta_key or self.standardized_beta_key,
            self.ss_var_beta_key or self.standardized_var_beta_key,
            self.ss_sdy_key or self.standardized_sdy_key,
            self.ss_se_key or self.standardized_se_key,
            self.ss_maf_key or self.dbSNP_maf_key,
            self.ss_n_key or self.standardized_n_key,
            self.ss_tissue_key or self.standardized_tissue_key,
            self.dbSNP_validation_key
        ]])

        # Rename columns to standardized names
        self.transformed_df.columns = [
            self.standardized_chr_key, self.standardized_pos_key, self.standardized_variant_id_key, self.standardized_rsid_key,
            self.standardized_gene_id_key, self.standardized_non_effect_allele_key, self.standardized_effect_allele_key,
            self.standardized_p_key, self.standardized_beta_key, self.standardized_var_beta_key, self.standardized_sdy_key,
            self.standardized_se_key, self.standardized_maf_key, self.standardized_n_key, self.standardized_tissue_key, self.dbSNP_validation_key
        ]


    def save_output(self):
        """
        Save the standardized, failed-dbSNP, and mismatched-allele DataFrames to output files.

        Returns:
            None: The DataFrames are saved to TSV files in the specified output directory.

        """
        # Filter out rows with "NA" in the standardized variant ID column for the standardized output
        print(f"Number of SNPs that were unable to be standardized: {len(self.transformed_df[self.transformed_df[self.standardized_variant_id_key] == 'NA'])}")
        standardized_df = self.transformed_df[self.transformed_df[self.standardized_variant_id_key] != "NA"]
        print(f"Number of SNPs validated against dbSNP: {len(standardized_df[standardized_df[self.dbSNP_validation_key] == 'Validated'])}")

        # Filter out rows that failed dbSNP lookup for the failed dbSNP output
        annotated_df_failed_dbsnp = self.annotated_df[
            self.annotated_df.index.isin(self.transformed_df.index) &
            self.annotated_df.index.isin(standardized_df[standardized_df[self.dbSNP_validation_key] == "Failed dbSNP lookup"].index)
        ]
        print(f"Number of SNPs that failed dbSNP lookup: {len(annotated_df_failed_dbsnp)}")

        # Filter out rows with mismatched alleles for the mismatched alleles output
        annotated_df_mismatched = self.annotated_df[self.annotated_df[self.dbSNP_validation_key].isin(["Position Mismatch", "Allele Mismatch"])]
        print(f"Number of SNPs with mismatched alleles: {len(annotated_df_mismatched)}")

        # Determine the base name for the output files based on the input file name, handling .gz suffix if present
        if self.sum_stats_fp.suffix == ".gz":
            input_file_base = Path(self.sum_stats_fp.stem).stem
        else:
            input_file_base = self.sum_stats_fp.stem

        # Save the standardized DataFrame to a TSV file if it contains any rows
        if len(standardized_df) > 0:
            standardized_fp = Path(f"{self.output_dir}/{input_file_base}.standardized.tsv")
            standardized_df.to_csv(standardized_fp, sep="\t", index=False)
            print(f"Saved standardized summary statistics to {standardized_fp}")

        # Save the failed dbSNP DataFrame to a TSV file if it contains any rows
        if len(annotated_df_failed_dbsnp) > 0:
            failed_dbsnp_fp = Path(f"{self.output_dir}/{input_file_base}.failed_dbsnp.tsv")
            annotated_df_failed_dbsnp.to_csv(failed_dbsnp_fp, sep="\t", index=False)
            print(f"Saved failed dbSNP annotations to {failed_dbsnp_fp}")

        # Save the mismatched alleles DataFrame to a TSV file if it contains any rows
        if len(annotated_df_mismatched) > 0:
            mismatched_alleles_fp = Path(f"{self.output_dir}/{input_file_base}.mismatched_alleles.tsv")
            annotated_df_mismatched.to_csv(mismatched_alleles_fp, sep="\t", index=False)
            print(f"Saved mismatched alleles to {mismatched_alleles_fp}")


    def run(self):
        """
        Execute the full pipeline: filter -> annotate -> split mismatches -> transform -> save.

        Returns:
            None: The pipeline is executed in sequence, and output files are saved to the specified directory.

        """
        # Normalize chromosome column
        print("Normalizing chromosome column in summary statistics...")
        self.ss_df[self.standardized_chr_key] = self.ss_df[self.ss_chr_key].astype(str).map(normalize_chromosome)

        # Filter to locus regions
        print("Filtering summary statistics to specified locus regions...")
        self.filter_by_expanded_ranges()
        print(f"Filtered summary stats shape: {self.ss_df.shape}")

        # Annotate
        print("Annotating summary statistics with dbSNP information...")
        self.annotate()

        # Transform to standardized layout
        print("Transforming annotated summary statistics to standardized layout...")
        self.transform()

        # Save output files
        print("Saving output files...")
        self.save_output()


# ------------------------- Full Pipeline -------------------------

def main(args: argparse.Namespace) -> None:
    """
    Main function to load input data and execute the full summary statistics
    transformation pipeline via SumStatsTransformer.

    Args:
        args (argparse.Namespace): Parsed command-line arguments.

    """
    # Initialize Entrez email (required by NCBI) for API access
    Entrez.email = args.entrez_email

    # Instantiate transformer with user-supplied column names
    transformer = SumStatsTransformer(args=args)

    # Run the full pipeline
    transformer.run()


# ------------------------- Command Line Interface -------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Transform summary statistics to standardized format.")
    parser.add_argument("--sum_stats_file", type=str, required=True, help="Path to the input file containing summary statistics.")
    parser.add_argument("--loci_file", type=str, required=True, help="Name of the file containing loci information. Output from define_loci.")
    parser.add_argument("--output_dir", type=str, required=True, help="Directory to save the output files.")
    parser.add_argument("--header_lines", type=int, default=0, help="Number of header lines to skip in the input file.")
    parser.add_argument("--entrez_email", type=str, required=True, help="Email address for NCBI Entrez API.")
    parser.add_argument("--ss_chr_key", type=str, required=True, help="Column name for chromosome.")
    parser.add_argument("--ss_pos_key", type=str, required=True, help="Column name for genomic position.")
    parser.add_argument("--ss_rsid_key", type=str, default=None, help="Column name for rsID (optional).")
    parser.add_argument("--ss_gene_id_key", type=str, default=None, help="Column name for gene ID.")
    parser.add_argument("--ss_non_effect_allele_key", type=str, default=None, help="Column name for non-effect allele (optional).")
    parser.add_argument("--ss_effect_allele_key", type=str, default=None, help="Column name for effect allele (optional).")
    parser.add_argument("--ss_p_key", type=str, default=None, help="Column name for p-value (optional).")
    parser.add_argument("--ss_beta_key", type=str, default=None, help="Column name for beta coefficient (optional).")
    parser.add_argument("--ss_se_key", type=str, default=None, help="Column name for standard error (optional).")
    parser.add_argument("--ss_n_key", type=str, default=None, help="Column name for sample size (optional).")
    parser.add_argument("--ss_statistic_key", type=str, default=None, help="Column name for test statistic (optional).")
    parser.add_argument("--ss_maf_key", type=str, default=None, help="Column name for minor allele frequency (optional).")
    parser.add_argument("--ss_mac_key", type=str, default=None, help="Column name for minor allele count (optional).")
    parser.add_argument("--ss_var_beta_key", type=str, default=None, help="Column name for variance of beta (optional).")
    parser.add_argument("--ss_sdy_key", type=str, default=None, help="Column name for standard deviation of Y (optional).")
    parser.add_argument("--ss_tissue_key", type=str, default=None, help="Column name for tissue label (optional).")
    parser.add_argument("--standardized_chr_key", type=str, required=True, help="Standardized column name for chromosome.")
    parser.add_argument("--standardized_pos_key", type=str, required=True, help="Standardized column name for position.")
    parser.add_argument("--standardized_rsid_key", type=str, required=True, help="Standardized column name for rsID.")
    parser.add_argument("--standardized_variant_id_key", type=str, required=True, help="Standardized column name for variant ID.")
    parser.add_argument("--standardized_gene_id_key", type=str, required=True, help="Standardized column name for gene ID.")
    parser.add_argument("--standardized_non_effect_allele_key", type=str, required=True, help="Standardized column name for non-effect allele.")
    parser.add_argument("--standardized_effect_allele_key", type=str, required=True, help="Standardized column name for effect allele.")
    parser.add_argument("--standardized_p_key", type=str, required=True, help="Standardized column name for p-value.")
    parser.add_argument("--standardized_beta_key", type=str, required=True, help="Standardized column name for beta coefficient.")
    parser.add_argument("--standardized_se_key", type=str, required=True, help="Standardized column name for standard error.")
    parser.add_argument("--standardized_n_key", type=str, required=True, help="Standardized column name for sample size.")
    parser.add_argument("--standardized_statistic_key", type=str, required=True, help="Standardized column name for test statistic.")
    parser.add_argument("--standardized_maf_key", type=str, required=True, help="Standardized column name for minor allele frequency.")
    parser.add_argument("--standardized_var_beta_key", type=str, required=True, help="Standardized column name for variance of beta.")
    parser.add_argument("--standardized_sdy_key", type=str, required=True, help="Standardized column name for standard deviation of Y.")
    parser.add_argument("--standardized_tissue_key", type=str, required=True, help="Standardized column name for tissue label.")
    parser.add_argument("--loci_left_bound_key", type=str, required=True, help="Column name for left bound of loci.")
    parser.add_argument("--loci_right_bound_key", type=str, required=True, help="Column name for right bound of loci.")

    args = parser.parse_args()
    main(args)
