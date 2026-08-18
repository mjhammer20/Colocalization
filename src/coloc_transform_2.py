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


# ------------------------- Module-Level Helper Functions -------------------------
# These functions are stateless utilities used internally by SumStatsTransformer.
# They are defined at module level so they can be cached/shared across instances.

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

    if isinstance(ser, pd.Series):
        return ser
    return pd.Series(ser, index=df.index)


@lru_cache(maxsize=65536)
def fetch_rsid(chromosome: str, position: int) -> str:
    """
    Query NCBI dbSNP to retrieve rsID information based on chromosome and position.

    Args:
        chromosome (str): The chromosome number or identifier.
        position (int): The base pair position on the chromosome.

    Returns:
        str: The corresponding rsID if found, otherwise "NA".

    """
    search_term = f"{chromosome}[Chromosome] AND {position}[Base Position] AND Homo sapiens[Organism]"
    handle = Entrez.esearch(db="snp", term=search_term)
    response = xmltodict.parse(handle) # type: ignore (silences pylance warning)
    handle.close()
    rs_ids = response.get("IdList", None)
    if not rs_ids:
        return "NA"
    return f'rs{str(rs_ids[-1]).strip()}'


def fetch_rsid_records(rsids: list) -> dict:
    """
    Query NCBI dbSNP to retrieve allele information based on a list of rsIDs.

    Args:
        rsids (list): The list of rsIDs to fetch.

    Returns:
        dict: Dictionary containing the fetched rsID records.

    """
    rsids_str = ','.join(rsids)
    attempts = 0
    while attempts < 5:
        attempts += 1
        handle = None
        try:
            handle = Entrez.efetch(db="snp", id=rsids_str, retmode="xml")
            response = xmltodict.parse(handle) # type: ignore (silences pylance warning)
            handle.close()
            records = response.get('ExchangeSet', {}).get('DocumentSummary', {})
        except (IncompleteRead, HTTPError, URLError, ExpatError) as e:
            if handle is not None:
                handle.close()
            if attempts < 5:
                time.sleep(2 * attempts)
                continue
            print(f"fetch_rsid_records: reached max retries, returning empty result for rsids: {rsids[0]} - {rsids[-1]}")
            return {}
        except Exception as e:
            print(f"Error fetching records from dbSNP: {e}")
            if handle is not None:
                handle.close()
            return {}
    return records


def extract_alleles_from_record(record: dict) -> tuple:
    """
    Extract the reference and alternate alleles from an NCBI dbSNP record.

    Args:
        record (dict): The NCBI dbSNP record.

    Returns:
        tuple: A tuple of (ref_allele, alt_allele), or ("NA", "NA") if not found.

    """
    docsum = record.get('DOCSUM', None)
    if not docsum:
        print(f"Error: DOCSUM not found in record for rsID rs{record.get('@uid', 'Unknown')}")
        return "NA", "NA"
    docsum_dict = {k.strip(): v.strip() for k, v in (s.split("=", 1) for s in docsum.split('|'))}
    ref_alt = docsum_dict.get('SEQ', None)
    if not ref_alt:
        print(f"Error: SEQ not found in DOCSUM for rsID rs{record.get('@uid', 'Unknown')}")
        return "NA", "NA"
    ref_alt_list = ref_alt.split('/')
    ref = ref_alt_list[0].replace('[', '')
    alt = ref_alt_list[1].replace(']', '')
    return ref, alt


def extract_maf_from_record(record: dict) -> float:
    """
    Extract the minor allele frequency (MAF) from an NCBI dbSNP record.

    Args:
        record (dict): The NCBI dbSNP record.

    Returns:
        float: The MAF value, or NaN if not found.

    """
    global_mafs = record.get('GLOBAL_MAFS', None)
    if not global_mafs:
        print(f"Error: GLOBAL_MAFS not found in record for rsID rs{record.get('@uid', 'Unknown')}")
        return np.nan
    mafs_list = global_mafs.get('MAF', None)
    if not mafs_list:
        print(f"Error: MAFs not found in GLOBAL_MAFS for rsID rs{record.get('@uid', 'Unknown')}")
        return np.nan
    try:
        mafs_dict = {x['STUDY']: x['FREQ'] for x in mafs_list}
    except TypeError:
        try:
            mafs_dict = {mafs_list['STUDY']: mafs_list['FREQ']}
        except Exception as e:
            print(f"Error: MAFs field not in expected format for rsID rs{record.get('@uid', 'Unknown')}: {e}")
            return np.nan
    maf = mafs_dict.get('GnomAD_genomes') or mafs_dict.get('1000Genomes') or list(mafs_dict.values())[0]
    maf_parsed = maf.split("=")
    freqs = maf_parsed[1].split("/")
    return float(freqs[0])


# ------------------------- Class Definition -------------------------

class SumStatsTransformer:
    """
    Encapsulates the full summary statistics annotation and transformation pipeline.

    Column name arguments provided at construction are stored as instance variables
    and reused across all pipeline methods, eliminating the need to pass them
    repeatedly to individual functions.
    """

    def __init__(self, args: argparse.Namespace):
        self.output_dir = Path(args.output_dir)
        self.sum_stats_fp = Path(args.sum_stats_file)
        self.ss_df = load_input_data(self.sum_stats_fp, args.header_lines)
        print(f"Loaded input data from {self.sum_stats_fp} with shape: {self.ss_df.shape}")
        self.loci_fp = Path(self.output_dir, args.loci_file)
        self.loci_df = load_input_data(self.loci_fp, 0)
        print(f"Loaded input data from {self.loci_fp} with shape: {self.loci_df.shape}")
        self.standardized_chr_col = args.standardized_chr_col
        self.loci_left_bound_col = args.loci_left_bound_col
        self.loci_right_bound_col = args.loci_right_bound_col
        self.snps_df = pd.DataFrame()
        self.ss_chr_col = args.ss_chr_col
        self.ss_pos_col = args.ss_pos_col
        self.ss_rsid_col = args.ss_rsid_col
        self.ss_non_effect_allele_col = args.ss_non_effect_allele_col
        self.ss_effect_allele_col = args.ss_effect_allele_col
        self.ss_n_col = args.ss_n_col
        self.ss_maf_col = args.ss_maf_col
        self.ss_mac_col = args.ss_mac_col
        self.ss_p_col = args.ss_p_col
        self.ss_statistic_col = args.ss_statistic_col
        self.ss_se_col = args.ss_se_col
        self.ss_beta_col = args.ss_beta_col
        self.ss_var_beta_col = args.ss_var_beta_col
        self.ss_sdY_col = args.ss_sdY_col
        self.ss_gene_id_col = args.ss_gene_id_col
        self.ss_tissue_col = args.ss_tissue_col
        self.annotated_df = pd.DataFrame()
        self.allele_match_col = "Allele_Match"
        self.standardized_ref_col = args.standardized_ref_col
        self.standardized_alt_col = args.standardized_alt_col
        self.standardized_maf_col = args.standardized_maf_col
        self.standardized_pos_col = args.standardized_pos_col
        self.standardized_variant_id_col = args.standardized_variant_id_col
        self.standardized_rsid_col = args.standardized_rsid_col
        self.standardized_gene_id_col = args.standardized_gene_id_col
        self.standardized_non_effect_allele_col = args.standardized_non_effect_allele_col
        self.standardized_effect_allele_col = args.standardized_effect_allele_col
        self.standardized_n_col = args.standardized_n_col
        self.standardized_p_col = args.standardized_p_col
        self.standardized_se_col = args.standardized_se_col
        self.standardized_beta_col = args.standardized_beta_col
        self.standardized_var_beta_col = args.standardized_var_beta_col
        self.standardized_sdY_col = args.standardized_sdY_col
        self.standardized_statistic_col = args.standardized_statistic_col
        self.standardized_tissue_col = args.standardized_tissue_col
        self.annotated_df_mismatched = pd.DataFrame()
        

    def filter_by_expanded_ranges(self):
        """
        Filter summary statistics to include only rows that fall within specified expanded locus ranges.

        Returns:
            None: The DataFrame is modified in place to include only rows within the specified ranges.

        """
        mask = pd.Series(False, index=self.ss_df.index)

        for _, row in self.loci_df.iterrows():
            idx = self.ss_df.index[self.ss_df[self.standardized_chr_col] == row[self.standardized_chr_col]]
            pos_vals = self.ss_df[self.ss_pos_col].loc[idx]  # type: pd.Series
            chrom_mask = pd.Series(False, index=idx)
            chrom_mask |= pos_vals.between(row[self.loci_left_bound_col], row[self.loci_right_bound_col])
            mask.loc[idx] = chrom_mask

        self.ss_df = pd.DataFrame(self.ss_df[mask])

    def add_dbSNP_info(self):
        """
        Query NCBI dbSNP in batches to retrieve allele information and MAF for each rsID in the index.

        Returns:
            None: The DataFrame is modified in place to include dbSNP allele and MAF information.

        """
        rsids_full = self.snps_df.index.tolist()
        batch_size = 200
        alleles_buf = []
        mafs_buf = []

        for start in range(0, len(rsids_full), batch_size):
            end = min(start + batch_size, len(rsids_full))
            batch_rsids = rsids_full[start:end]

            batch_records = fetch_rsid_records(batch_rsids)
            indices = [f"rs{record.get('@uid')}" for record in batch_records]

            batch_alleles = pd.Series(
                [extract_alleles_from_record(record) for record in batch_records],
                index=indices
            )
            alleles_buf.append(batch_alleles)

            batch_mafs = pd.Series(
                [extract_maf_from_record(record) for record in batch_records],
                index=indices
            )
            mafs_buf.append(batch_mafs)

            del batch_records, batch_alleles, batch_mafs
            gc.collect()

        alleles = pd.concat(alleles_buf)
        mafs = pd.concat(mafs_buf)

        alleles_expanded = pd.DataFrame(
            data=list(alleles),
            index=alleles.index,
            columns=[self.standardized_ref_col, self.standardized_alt_col]
        )

        self.snps_df.loc[mafs.index, self.standardized_maf_col] = mafs
        self.snps_df = self.ss_df.join(alleles_expanded, how='left')


    def check_allele_match(self):
        """
        Check whether the alleles from the summary statistics match the alleles from dbSNP.

        Returns:
            None: The DataFrame is modified in place to include a new column indicating allele match status.

        """
        self.annotated_df[self.allele_match_col] = np.where(
            self.annotated_df[[self.standardized_ref_col, self.standardized_alt_col, self.ss_non_effect_allele_col, self.ss_effect_allele_col]].notna().all(axis=1), #type: ignore (silences pylance warning)
            (self.annotated_df[self.standardized_ref_col] == self.annotated_df[self.ss_non_effect_allele_col]) & (self.annotated_df[self.standardized_alt_col] == self.annotated_df[self.ss_effect_allele_col]),
            "NA"
        )

    def add_variant_id(self):
        """
        Create a variant ID in the format chr:bp:ref:alt for each row in the DataFrame.

        Returns:
            None: The DataFrame is modified in place to include a new column with the variant ID.

        """
        self.annotated_df[self.standardized_variant_id_col] = np.where(
            self.annotated_df[[self.standardized_ref_col, self.standardized_alt_col]].notna().all(axis=1), #type: ignore (silences pylance warning)
            self.annotated_df[self.standardized_chr_col].astype(str) + ':' + self.annotated_df[self.ss_pos_col].astype(str) + ':' +
            self.annotated_df[self.standardized_ref_col].astype(str) + ':' + self.annotated_df[self.standardized_alt_col].astype(str),
            "NA"
        )

    def calculate_missing_summary_statistics(self):
        """
        Calculate missing summary statistics (N, statistic, SE, beta, var(beta), SDY, P)
        from available columns, filling in defaults where possible.

        Args:
            df (pd.DataFrame): The input DataFrame.

        Returns:
            pd.DataFrame: The DataFrame with missing summary statistics filled in.

        """
        s_n = _numeric_series(self.annotated_df, self.ss_n_col)
        s_maf = _numeric_series(self.annotated_df, self.ss_maf_col)
        s_mac = _numeric_series(self.annotated_df, self.ss_mac_col)
        s_p = _numeric_series(self.annotated_df, self.ss_p_col)
        s_statistic = _numeric_series(self.annotated_df, self.ss_statistic_col)
        s_se = _numeric_series(self.annotated_df, self.ss_se_col)
        s_beta = _numeric_series(self.annotated_df, self.ss_beta_col)
        s_var_beta = _numeric_series(self.annotated_df, self.ss_var_beta_col)

        # N
        if self.ss_n_col is None:
            s_n = pd.Series(np.nan, index=self.annotated_df.index)
            m = s_mac.notna() & s_maf.notna() & (s_maf != 0)
            s_n.loc[m] = (s_mac / (2 * s_maf))[m]
            self.annotated_df[self.standardized_n_col] = s_n

        # Statistic
        if self.ss_statistic_col is None:
            s_statistic = pd.Series(np.nan, index=s_p.index)
            m = s_p.notna() & (s_p > 0) & (s_p < 1)
            s_statistic.loc[m] = (norm.ppf(1 - s_p / 2))[m]
            self.annotated_df[self.standardized_statistic_col] = s_statistic

        # SE
        if not self.ss_se_col:
            s_se = pd.Series(np.nan, index=self.annotated_df.index)
            m = s_beta.notna() & s_statistic.notna() & (s_statistic != 0)
            s_se.loc[m] = (s_beta.abs() / s_statistic)[m]
            m = s_se.isna() & s_maf.notna() & s_n.notna() & (s_n > 0) & s_statistic.notna() & (s_statistic != 0)
            s_se.loc[m] = 1 / np.sqrt(2 * s_maf * (1 - s_maf) * (s_n + s_statistic**2))[m]
            self.annotated_df[self.standardized_se_col] = s_se

        # Beta
        if not self.ss_beta_col:
            s_beta = pd.Series(np.nan, index=self.annotated_df.index)
            m = s_se.notna() & s_statistic.notna()
            s_beta.loc[m] = (s_statistic * s_se)[m]
            self.annotated_df[self.standardized_beta_col] = s_beta

        # var(beta)
        if not self.ss_var_beta_col:
            s_var_beta = pd.Series(np.nan, index=self.annotated_df.index)
            m = s_se.notna() & (s_se != 0)
            s_var_beta.loc[m] = (s_se ** 2)[m]
            self.annotated_df[self.standardized_var_beta_col] = s_var_beta

        # SDY
        if not self.ss_sdY_col:
            sdY = pd.Series(np.nan, index=self.annotated_df.index)
            m = s_maf.notna() & s_n.notna() & (s_n > 0)
            sdY.loc[m] = np.sqrt(s_var_beta * s_n * 2 * s_maf * (1 - s_maf))[m]
            self.annotated_df[self.standardized_sdY_col] = sdY

        # P-value
        if not self.ss_p_col:
            s_p = pd.Series(np.nan, index=self.annotated_df.index)
            m = s_p.isna() & s_statistic.notna()
            s_p.loc[m] = 2 * norm.sf(s_statistic.abs())[m]
            self.annotated_df[self.standardized_p_col] = s_p


    def annotate(self):
        """
        Annotate a summary statistics DataFrame with dbSNP allele information, allele
        match flags, variant IDs, and any missing derived summary statistics.

        Returns:
            None: The DataFrame is modified in place to include the new annotations.

        """
        # Subset to unique SNPs and fetch rsID if not provided
        if self.ss_rsid_col is None:
            self.snps_df = self.ss_df[[self.standardized_chr_col, self.ss_pos_col]].drop_duplicates()
            self.snps_df[self.standardized_rsid_col] = self.snps_df.apply(
                lambda row: fetch_rsid(row[self.standardized_chr_col], row[self.ss_pos_col]), axis=1
            )
            self.snps_df = self.snps_df[self.snps_df[self.standardized_rsid_col] != 'NA']
            self.snps_df.set_index(self.standardized_rsid_col, inplace=True)
        else:
            self.snps_df = self.ss_df[[self.standardized_chr_col, self.ss_pos_col, self.ss_rsid_col]].drop_duplicates()
            self.snps_df.set_index(self.ss_rsid_col, inplace=True)

        # Add dbSNP allele and MAF information
        if len(self.snps_df) > 0:
            self.add_dbSNP_info()

        # Merge dbSNP info back into the original DataFrame
        merge_cols = [self.standardized_chr_col, self.ss_pos_col] if self.ss_rsid_col not in self.ss_df.columns \
            else [self.standardized_chr_col, self.ss_pos_col, self.ss_rsid_col]
        self.annotated_df = pd.merge(self.ss_df, self.snps_df, on=merge_cols, how='left')

        # Use dbSNP MAF if maf_col was not provided
        if self.ss_maf_col is None and self.standardized_maf_col in self.annotated_df.columns:
            self.ss_maf_col = self.standardized_maf_col

        # Check allele match if effect/non-effect allele columns are available
        if self.ss_non_effect_allele_col is not None and self.ss_effect_allele_col is not None:
            self.check_allele_match()

        # Add variant ID
        if self.allele_match_col in self.annotated_df.columns:
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
        if not self.ss_tissue_col:
            self.annotated_df[self.standardized_tissue_col] = 'NA'

        self.transformed_df = pd.DataFrame(self.annotated_df[[
            self.standardized_chr_col,
            self.ss_pos_col,
            self.ss_rsid_col or self.standardized_variant_id_col,
            self.ss_rsid_col or self.standardized_rsid_col,
            self.ss_gene_id_col,
            self.standardized_ref_col or self.ss_non_effect_allele_col,
            self.standardized_alt_col or self.ss_effect_allele_col,
            self.ss_p_col or self.standardized_p_col,
            self.ss_beta_col or self.standardized_beta_col,
            self.ss_var_beta_col or self.standardized_var_beta_col,
            self.ss_sdY_col or self.standardized_sdY_col,
            self.ss_se_col or self.standardized_se_col,
            self.ss_maf_col or self.standardized_maf_col,
            self.ss_n_col or self.standardized_n_col,
            self.ss_tissue_col or self.standardized_tissue_col
        ]])

        self.transformed_df.columns = [
            self.standardized_chr_col, self.standardized_pos_col, self.standardized_variant_id_col, self.standardized_rsid_col,
            self.standardized_gene_id_col, self.standardized_non_effect_allele_col, self.standardized_effect_allele_col,
            self.standardized_p_col, self.standardized_beta_col, self.standardized_var_beta_col, self.standardized_sdY_col,
            self.standardized_se_col, self.standardized_maf_col, self.standardized_n_col, self.standardized_tissue_col
        ]


    def save_output(self):
        """
        Save the standardized, failed-dbSNP, and mismatched-allele DataFrames to output files.

        Returns:
            None: The DataFrames are saved to TSV files in the specified output directory.

        """
        standardized_df = self.transformed_df[self.transformed_df[self.standardized_variant_id_col] != 'NA']
        annotated_df_failed_dbsnp = self.annotated_df[
            self.annotated_df.index.isin(self.transformed_df.index) &
            ~self.annotated_df.index.isin(standardized_df.index)
        ]

        if self.sum_stats_fp.suffix == '.gz':
            input_file_base = Path(self.sum_stats_fp.stem).stem
        else:
            input_file_base = self.sum_stats_fp.stem

        if len(standardized_df) > 0:
            standardized_fp = Path(f"{self.output_dir}/{input_file_base}.standardized.tsv")
            standardized_df.to_csv(standardized_fp, sep='\t', index=False)

        if len(annotated_df_failed_dbsnp) > 0:
            failed_dbsnp_fp = Path(f"{self.output_dir}/{input_file_base}.failed_dbsnp.tsv")
            annotated_df_failed_dbsnp.to_csv(failed_dbsnp_fp, sep='\t', index=False)

        if len(self.annotated_df_mismatched) > 0:
            mismatched_alleles_fp = Path(f"{self.output_dir}/{input_file_base}.mismatched_alleles.tsv")
            self.annotated_df_mismatched.to_csv(mismatched_alleles_fp, sep='\t', index=False)

    def run(self):
        """
        Execute the full pipeline: filter -> annotate -> split mismatches -> transform -> save.

        Returns:
            None: The pipeline is executed in sequence, and output files are saved to the specified directory.

        """
        # Normalize chromosome column
        self.ss_df[self.standardized_chr_col] = self.ss_df[self.ss_chr_col].astype(str).map(normalize_chromosome)

        # Filter to locus regions
        self.filter_by_expanded_ranges()
        print(f"Filtered summary stats shape: {self.ss_df.shape}")

        # Annotate
        self.annotate()

        # Split out mismatched alleles
        if self.allele_match_col in self.annotated_df.columns:
            self.annotated_df_mismatched = pd.DataFrame(
                self.annotated_df[self.annotated_df[self.allele_match_col] == False]
            )
            self.annotated_df = self.annotated_df[~self.annotated_df.index.isin(self.annotated_df_mismatched.index)]

        # Transform to standardized layout
        self.transform()

        # Save output files
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
    parser.add_argument("--ss_chr_col", type=str, required=True, help="Column name for chromosome.")
    parser.add_argument("--ss_pos_col", type=str, required=True, help="Column name for genomic position.")
    parser.add_argument("--ss_rsid_col", type=str, default=None, help="Column name for rsID (optional).")
    parser.add_argument("--ss_gene_id_col", type=str, required=True, help="Column name for gene ID.")
    parser.add_argument("--ss_non_effect_allele_col", type=str, default=None, help="Column name for non-effect allele (optional).")
    parser.add_argument("--ss_effect_allele_col", type=str, default=None, help="Column name for effect allele (optional).")
    parser.add_argument("--ss_p_col", type=str, default=None, help="Column name for p-value (optional).")
    parser.add_argument("--ss_beta_col", type=str, default=None, help="Column name for beta coefficient (optional).")
    parser.add_argument("--ss_se_col", type=str, default=None, help="Column name for standard error (optional).")
    parser.add_argument("--ss_n_col", type=str, default=None, help="Column name for sample size (optional).")
    parser.add_argument("--ss_statistic_col", type=str, default=None, help="Column name for test statistic (optional).")
    parser.add_argument("--ss_maf_col", type=str, default=None, help="Column name for minor allele frequency (optional).")
    parser.add_argument("--ss_mac_col", type=str, default=None, help="Column name for minor allele count (optional).")
    parser.add_argument("--ss_var_beta_col", type=str, default=None, help="Column name for variance of beta (optional).")
    parser.add_argument("--ss_sdY_col", type=str, default=None, help="Column name for standard deviation of Y (optional).")
    parser.add_argument("--ss_tissue_col", type=str, default=None, help="Column name for tissue label (optional).")
    parser.add_argument("--standardized_chr_col", type=str, required=True, help="Standardized column name for chromosome.")
    parser.add_argument("--standardized_pos_col", type=str, required=True, help="Standardized column name for position.")
    parser.add_argument("--standardized_rsid_col", type=str, required=True, help="Standardized column name for rsID.")
    parser.add_argument("--standardized_variant_id_col", type=str, required=True, help="Standardized column name for variant ID.")
    parser.add_argument("--standardized_gene_id_col", type=str, required=True, help="Standardized column name for gene ID.")
    parser.add_argument("--standardized_non_effect_allele_col", type=str, required=True, help="Standardized column name for non-effect allele.")
    parser.add_argument("--standardized_effect_allele_col", type=str, required=True, help="Standardized column name for effect allele.")
    parser.add_argument("--standardized_p_col", type=str, required=True, help="Standardized column name for p-value.")
    parser.add_argument("--standardized_beta_col", type=str, required=True, help="Standardized column name for beta coefficient.")
    parser.add_argument("--standardized_se_col", type=str, required=True, help="Standardized column name for standard error.")
    parser.add_argument("--standardized_n_col", type=str, required=True, help="Standardized column name for sample size.")
    parser.add_argument("--standardized_statistic_col", type=str, required=True, help="Standardized column name for test statistic.")
    parser.add_argument("--standardized_maf_col", type=str, required=True, help="Standardized column name for minor allele frequency.")
    parser.add_argument("--standardized_var_beta_col", type=str, required=True, help="Standardized column name for variance of beta.")
    parser.add_argument("--standardized_sdY_col", type=str, required=True, help="Standardized column name for standard deviation of Y.")
    parser.add_argument("--standardized_tissue_col", type=str, required=True, help="Standardized column name for tissue label.")
    parser.add_argument("--loci_left_bound_col", type=str, required=True, help="Column name for left bound of loci.")
    parser.add_argument("--loci_right_bound_col", type=str, required=True, help="Column name for right bound of loci.")

    args = parser.parse_args()
    main(args)
