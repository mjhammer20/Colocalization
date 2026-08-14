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


def load_expanded_ranges(fp: Path) -> dict:
    """
    Load expanded locus ranges from a JSON file into a dictionary.

    Args:
        fp (Path): Path to the JSON file containing expanded ranges.

    Returns:
        dict: A dictionary with chromosome as keys and a list of [start, end] ranges as values.

    """
    with open(fp, 'r') as f:
        expanded_ranges = json.load(f)
    return expanded_ranges


# ------------------------- Class Definition -------------------------

class SumStatsTransformer:
    """
    Encapsulates the full summary statistics annotation and transformation pipeline.

    Column name arguments provided at construction are stored as instance variables
    and reused across all pipeline methods, eliminating the need to pass them
    repeatedly to individual functions.

    Attributes:
        chr_col (str): Column name for chromosome.
        pos_col (str): Column name for genomic position.
        rsid_col (Optional[str]): Column name for rsID.
        non_effect_allele_col (Optional[str]): Column name for the non-effect allele.
        effect_allele_col (Optional[str]): Column name for the effect allele.
        n_col (Optional[str]): Column name for sample size.
        maf_col (Optional[str]): Column name for minor allele frequency.
        mac_col (Optional[str]): Column name for minor allele count.
        p_col (Optional[str]): Column name for p-value.
        statistic_col (Optional[str]): Column name for test statistic.
        se_col (Optional[str]): Column name for standard error.
        beta_col (Optional[str]): Column name for effect size (beta).
        var_beta_col (Optional[str]): Column name for variance of beta.
        sdY_col (Optional[str]): Column name for standard deviation of Y.
        gene_id_col (Optional[str]): Column name for gene ID.
        tissue_col (Optional[str]): Column name for tissue label.

    Default derived column names (set in __init__):
        allele_match_col (str): Column name for allele match flag. Default: "Allele_Match".
        variant_id_col (str): Column name for variant ID. Default: "Variant_ID".
        dbSNP_ref_col (str): Column name for dbSNP reference allele. Default: "ref_dbSNP".
        dbSNP_alt_col (str): Column name for dbSNP alternate allele. Default: "alt_dbSNP".
        dbSNP_maf_col (str): Column name for dbSNP MAF. Default: "MAF".
        default_rsid_col (str): Fallback rsID column name. Default: "rsID".
        default_n_col (str): Fallback N column name. Default: "N".
        default_p_col (str): Fallback p-value column name. Default: "P".
        default_se_col (str): Fallback SE column name. Default: "SE".
        default_beta_col (str): Fallback beta column name. Default: "BETA".
        default_var_beta_col (str): Fallback var(beta) column name. Default: "VARBETA".
        default_sdY_col (str): Fallback SDY column name. Default: "SDY".
        default_statistic_col (str): Fallback statistic column name. Default: "STATISTIC".
        default_tissue_col (str): Fallback tissue column name. Default: "TISSUE".
        default_maf_col (str): Fallback MAF column name. Default: "MAF".
        default_ref_col (str): Fallback reference allele column name. Default: "ref_dbSNP".
        default_alt_col (str): Fallback alternate allele column name. Default: "alt_dbSNP".

    """

    def __init__(
        self,
        chr_col: str,
        pos_col: str,
        rsid_col: Optional[str] = None,
        non_effect_allele_col: Optional[str] = None,
        effect_allele_col: Optional[str] = None,
        n_col: Optional[str] = None,
        maf_col: Optional[str] = None,
        mac_col: Optional[str] = None,
        p_col: Optional[str] = None,
        statistic_col: Optional[str] = None,
        se_col: Optional[str] = None,
        beta_col: Optional[str] = None,
        var_beta_col: Optional[str] = None,
        sdY_col: Optional[str] = None,
        gene_id_col: Optional[str] = None,
        tissue_col: Optional[str] = None
    ):
        # User-supplied column names
        self.chr_col               = chr_col
        self.pos_col               = pos_col
        self.rsid_col              = rsid_col
        self.non_effect_allele_col = non_effect_allele_col
        self.effect_allele_col     = effect_allele_col
        self.n_col                 = n_col
        self.maf_col               = maf_col
        self.mac_col               = mac_col
        self.p_col                 = p_col
        self.statistic_col         = statistic_col
        self.se_col                = se_col
        self.beta_col              = beta_col
        self.var_beta_col          = var_beta_col
        self.sdY_col               = sdY_col
        self.gene_id_col           = gene_id_col
        self.tissue_col            = tissue_col

        # Derived / default column names
        self.allele_match_col    = "Allele_Match"
        self.variant_id_col      = "Variant_ID"
        self.dbSNP_ref_col       = "ref_dbSNP"
        self.dbSNP_alt_col       = "alt_dbSNP"
        self.dbSNP_maf_col       = "MAF"
        self.default_rsid_col    = "rsID"
        self.default_n_col       = "N"
        self.default_p_col       = "P"
        self.default_se_col      = "SE"
        self.default_beta_col    = "BETA"
        self.default_var_beta_col = "VARBETA"
        self.default_sdY_col     = "SDY"
        self.default_statistic_col = "STATISTIC"
        self.default_tissue_col  = "TISSUE"
        self.default_maf_col     = "MAF"
        self.default_ref_col     = "ref_dbSNP"
        self.default_alt_col     = "alt_dbSNP"

    # ----------------------- Private Helper Methods -----------------------

    def _filter_by_expanded_ranges(self, df: pd.DataFrame, expanded_ranges: dict) -> pd.DataFrame:
        """
        Filter a DataFrame to include only rows that fall within specified expanded locus ranges.

        Args:
            df (pd.DataFrame): The input DataFrame.
            expanded_ranges (dict): A dictionary with chromosome as keys and a list of
                [start, end] ranges as values.

        Returns:
            pd.DataFrame: The filtered DataFrame containing only rows within the specified ranges.

        """
        mask = pd.Series(False, index=df.index)

        for chrom, ranges in expanded_ranges.items():
            idx = df.index[df[self.chr_col] == chrom]
            pos_vals = df[self.pos_col].loc[idx]  # type: pd.Series
            chrom_mask = pd.Series(False, index=idx)
            for start, end in ranges:
                chrom_mask |= pos_vals.between(start, end)
            mask.loc[idx] = chrom_mask

        return pd.DataFrame(df[mask])

    def _add_dbSNP_info(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Query NCBI dbSNP in batches to retrieve allele information and MAF for each rsID in the index.

        Args:
            df (pd.DataFrame): The input DataFrame with rsID as the index.

        Returns:
            pd.DataFrame: The input DataFrame with additional columns for reference allele,
                alternate allele, and MAF sourced from dbSNP.

        """
        rsids_full = df.index.tolist()
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
            columns=[self.dbSNP_ref_col, self.dbSNP_alt_col]
        )

        df.loc[mafs.index, self.dbSNP_maf_col] = mafs
        df = df.join(alleles_expanded, how='left')

        return df

    def _check_allele_match(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Check whether the alleles from the summary statistics match the alleles from dbSNP.

        Args:
            df (pd.DataFrame): The input DataFrame containing allele columns from both
                the summary statistics and dbSNP.

        Returns:
            pd.DataFrame: The DataFrame with an additional allele match column indicating
                whether the non-effect/effect alleles match the dbSNP ref/alt alleles.

        """
        df[self.allele_match_col] = np.where(
            df[[self.dbSNP_ref_col, self.dbSNP_alt_col,
                self.non_effect_allele_col, self.effect_allele_col]].notna().all(axis=1), #type: ignore (silences pylance warning)
            (df[self.dbSNP_ref_col] == df[self.non_effect_allele_col]) &
            (df[self.dbSNP_alt_col] == df[self.effect_allele_col]),
            "NA"
        )
        return df

    def _add_variant_id(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Create a variant ID in the format chr:bp:ref:alt for each row in the DataFrame.

        Args:
            df (pd.DataFrame): The input DataFrame containing chromosome, position,
                and allele columns.

        Returns:
            pd.DataFrame: The DataFrame with an additional variant ID column.

        """
        df[self.variant_id_col] = np.where(
            df[[self.dbSNP_ref_col, self.dbSNP_alt_col]].notna().all(axis=1), #type: ignore (silences pylance warning)
            df[self.chr_col].astype(str) + ':' + df[self.pos_col].astype(str) + ':' +
            df[self.dbSNP_ref_col].astype(str) + ':' + df[self.dbSNP_alt_col].astype(str),
            "NA"
        )
        return df

    def _calculate_missing_summary_statistics(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Calculate missing summary statistics (N, statistic, SE, beta, var(beta), SDY, P)
        from available columns, filling in defaults where possible.

        Args:
            df (pd.DataFrame): The input DataFrame.

        Returns:
            pd.DataFrame: The DataFrame with missing summary statistics filled in.

        """
        s_n         = _numeric_series(df, self.n_col)
        s_maf       = _numeric_series(df, self.maf_col)
        s_mac       = _numeric_series(df, self.mac_col)
        s_p         = _numeric_series(df, self.p_col)
        s_statistic = _numeric_series(df, self.statistic_col)
        s_se        = _numeric_series(df, self.se_col)
        s_beta      = _numeric_series(df, self.beta_col)
        s_var_beta  = _numeric_series(df, self.var_beta_col)

        # N
        if self.n_col is None:
            self.n_col = self.default_n_col
            s_n = pd.Series(np.nan, index=df.index)
            m = s_mac.notna() & s_maf.notna() & (s_maf != 0)
            s_n.loc[m] = (s_mac / (2 * s_maf))[m]
            df[self.n_col] = s_n

        # Statistic
        if self.statistic_col is None:
            self.statistic_col = self.default_statistic_col
            s_statistic = pd.Series(np.nan, index=s_p.index)
            m = s_p.notna() & (s_p > 0) & (s_p < 1)
            s_statistic.loc[m] = (norm.ppf(1 - s_p / 2))[m]
            df[self.statistic_col] = s_statistic

        # SE
        if not self.se_col:
            self.se_col = self.default_se_col
            s_se = pd.Series(np.nan, index=df.index)
            m = s_beta.notna() & s_statistic.notna() & (s_statistic != 0)
            s_se.loc[m] = (s_beta.abs() / s_statistic)[m]
            m = s_se.isna() & s_maf.notna() & s_n.notna() & (s_n > 0) & s_statistic.notna() & (s_statistic != 0)
            s_se.loc[m] = 1 / np.sqrt(2 * s_maf * (1 - s_maf) * (s_n + s_statistic**2))[m]
            df[self.se_col] = s_se

        # Beta
        if not self.beta_col:
            self.beta_col = self.default_beta_col
            s_beta = pd.Series(np.nan, index=df.index)
            m = s_se.notna() & s_statistic.notna()
            s_beta.loc[m] = (s_statistic * s_se)[m]
            df[self.beta_col] = s_beta

        # var(beta)
        if not self.var_beta_col:
            self.var_beta_col = self.default_var_beta_col
            s_var_beta = pd.Series(np.nan, index=df.index)
            m = s_se.notna() & (s_se != 0)
            s_var_beta.loc[m] = (s_se ** 2)[m]
            df[self.var_beta_col] = s_var_beta

        # SDY
        if not self.sdY_col:
            self.sdY_col = self.default_sdY_col
            sdY = pd.Series(np.nan, index=df.index)
            m = s_maf.notna() & s_n.notna() & (s_n > 0)
            sdY.loc[m] = np.sqrt(s_var_beta * s_n * 2 * s_maf * (1 - s_maf))[m]
            df[self.sdY_col] = sdY

        # P-value
        if not self.p_col:
            self.p_col = self.default_p_col
            s_p = pd.Series(np.nan, index=df.index)
            m = s_p.isna() & s_statistic.notna()
            s_p.loc[m] = 2 * norm.sf(s_statistic.abs())[m]
            df[self.p_col] = s_p

        return df

    # ----------------------- Public Pipeline Methods -----------------------

    def annotate(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Annotate a summary statistics DataFrame with dbSNP allele information, allele
        match flags, variant IDs, and any missing derived summary statistics.

        Args:
            df (pd.DataFrame): The raw input summary statistics DataFrame.

        Returns:
            pd.DataFrame: The annotated DataFrame with additional columns for rsID
                (if not provided), dbSNP alleles, MAF, allele match status, variant ID,
                and any imputed summary statistics.

        """
        # Subset to unique SNPs and fetch rsID if not provided
        if self.rsid_col is None:
            self.rsid_col = self.default_rsid_col
            snps_df = df[[self.chr_col, self.pos_col]].drop_duplicates()
            snps_df[self.rsid_col] = snps_df.apply(
                lambda row: fetch_rsid(row[self.chr_col], row[self.pos_col]), axis=1
            )
            snps_df = snps_df[snps_df[self.rsid_col] != 'NA']
            snps_df.set_index(self.rsid_col, inplace=True)
        else:
            snps_df = df[[self.chr_col, self.pos_col, self.rsid_col]].drop_duplicates()
            snps_df.set_index(self.rsid_col, inplace=True)

        # Add dbSNP allele and MAF information
        if len(snps_df) > 0:
            snps_df = pd.DataFrame(snps_df)
            snps_df = self._add_dbSNP_info(snps_df)

        # Merge dbSNP info back into the original DataFrame
        merge_cols = [self.chr_col, self.pos_col] if self.rsid_col not in df.columns \
            else [self.chr_col, self.pos_col, self.rsid_col]
        annotated_df = pd.merge(df, snps_df, on=merge_cols, how='left')

        # Use dbSNP MAF if maf_col was not provided
        if self.maf_col is None and self.dbSNP_maf_col in annotated_df.columns:
            self.maf_col = self.dbSNP_maf_col

        # Check allele match if effect/non-effect allele columns are available
        if self.non_effect_allele_col is not None and self.effect_allele_col is not None:
            annotated_df = self._check_allele_match(annotated_df)

        # Add variant ID
        annotated_df = self._add_variant_id(annotated_df)

        # Calculate missing summary statistics
        annotated_df = self._calculate_missing_summary_statistics(annotated_df)

        return annotated_df

    def transform(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Transform an annotated DataFrame into the standardized column layout expected
        by downstream colocalization tools.

        Args:
            df (pd.DataFrame): The annotated summary statistics DataFrame.

        Returns:
            pd.DataFrame: The transformed DataFrame with standardized column names:
                CHR, BP, VARIANT_ID, RSID, GENE_ID, NON_EFFECT_ALLELE, EFFECT_ALLELE,
                P, BETA, VAR_BETA, SDY, SE, MAF, N, TISSUE.

        """
        # Ensure tissue column exists
        if not self.tissue_col:
            self.tissue_col = self.default_tissue_col
            df[self.tissue_col] = 'NA'

        transformed_df = pd.DataFrame(df[[
            self.chr_col,
            self.pos_col,
            self.variant_id_col,
            self.rsid_col or self.default_rsid_col,
            self.gene_id_col or '',
            self.non_effect_allele_col or self.default_ref_col,
            self.effect_allele_col or self.default_alt_col,
            self.p_col or self.default_p_col,
            self.beta_col or self.default_beta_col,
            self.var_beta_col or self.default_var_beta_col,
            self.sdY_col or self.default_sdY_col,
            self.se_col or self.default_se_col,
            self.maf_col or self.default_maf_col,
            self.n_col or self.default_n_col,
            self.tissue_col
        ]])

        transformed_df.columns = [
            'CHR', 'BP', 'VARIANT_ID', 'RSID',
            'GENE_ID', 'NON_EFFECT_ALLELE', 'EFFECT_ALLELE',
            'P', 'BETA', 'VAR_BETA', 'SDY', 'SE', 'MAF', 'N', 'TISSUE'
        ]

        return transformed_df

    def save_output(
        self,
        annotated_df: pd.DataFrame,
        annotated_df_mismatched: pd.DataFrame,
        transformed_df: pd.DataFrame,
        output_dir: str,
        input_fp: Path
    ) -> None:
        """
        Save the standardized, failed-dbSNP, and mismatched-allele DataFrames to output files.

        Args:
            annotated_df (pd.DataFrame): The annotated DataFrame containing summary statistics.
            annotated_df_mismatched (pd.DataFrame): The DataFrame containing rows with
                mismatched alleles.
            transformed_df (pd.DataFrame): The transformed DataFrame with standardized columns.
            output_dir (str): The directory to save the output files.
            input_fp (Path): The input file path used to derive output file names.

        Returns:
            None

        """
        standardized_df = transformed_df[transformed_df['VARIANT_ID'] != 'NA']
        annotated_df_failed_dbsnp = annotated_df[
            annotated_df.index.isin(transformed_df.index) &
            ~annotated_df.index.isin(standardized_df.index)
        ]

        if input_fp.suffix == '.gz':
            input_file_base = Path(input_fp.stem).stem
        else:
            input_file_base = input_fp.stem

        if len(standardized_df) > 0:
            standardized_fp = Path(f"{output_dir}/{input_file_base}.standardized.tsv")
            standardized_df.to_csv(standardized_fp, sep='\t', index=False)

        if len(annotated_df_failed_dbsnp) > 0:
            failed_dbsnp_fp = Path(f"{output_dir}/{input_file_base}.failed_dbsnp.tsv")
            annotated_df_failed_dbsnp.to_csv(failed_dbsnp_fp, sep='\t', index=False)

        if len(annotated_df_mismatched) > 0:
            mismatched_alleles_fp = Path(f"{output_dir}/{input_file_base}.mismatched_alleles.tsv")
            annotated_df_mismatched.to_csv(mismatched_alleles_fp, sep='\t', index=False)

    def run(
        self,
        df: pd.DataFrame,
        expanded_ranges: dict,
        output_dir: str,
        input_fp: Path
    ) -> None:
        """
        Execute the full pipeline: filter -> annotate -> split mismatches -> transform -> save.

        Args:
            df (pd.DataFrame): The raw input summary statistics DataFrame.
            expanded_ranges (dict): Expanded locus ranges for filtering.
            output_dir (str): Directory to save output files.
            input_fp (Path): Input file path used to derive output file names.

        Returns:
            None

        """
        # Filter to locus regions
        filtered_df = self._filter_by_expanded_ranges(df, expanded_ranges)
        print(f"Filtered DataFrame shape: {filtered_df.shape}")

        # Annotate
        annotated_df = self.annotate(filtered_df)

        # Split out mismatched alleles
        if self.allele_match_col in annotated_df.columns:
            annotated_df_mismatched = pd.DataFrame(
                annotated_df[annotated_df[self.allele_match_col] == False]
            )
            annotated_df = annotated_df[~annotated_df.index.isin(annotated_df_mismatched.index)]
        else:
            annotated_df_mismatched = pd.DataFrame()

        # Transform to standardized layout
        transformed_df = self.transform(annotated_df)

        # Save output files
        self.save_output(
            annotated_df=annotated_df,
            annotated_df_mismatched=annotated_df_mismatched,
            transformed_df=transformed_df,
            output_dir=output_dir,
            input_fp=input_fp
        )


# ------------------------- Full Pipeline -------------------------

def main(args: argparse.Namespace) -> None:
    """
    Main function to load input data and execute the full summary statistics
    transformation pipeline via SumStatsTransformer.

    Args:
        args (argparse.Namespace): Parsed command-line arguments.

    """
    # Load input data
    input_fp = Path(args.input_file)
    df = load_input_data(input_fp, args.header_lines)
    print(f"Loaded input data from {input_fp} with shape: {df.shape}")

    # Load expanded locus ranges
    expanded_ranges = load_expanded_ranges(Path(args.expanded_ranges_file))

    # Normalize chromosome column
    df[args.chr_col] = df[args.chr_col].astype(str).map(normalize_chromosome)

    # Initialize Entrez email (required by NCBI) for API access
    Entrez.email = args.entrez_email

    # Instantiate transformer with user-supplied column names
    transformer = SumStatsTransformer(
        chr_col=args.chr_col,
        pos_col=args.pos_col,
        rsid_col=args.rsid_col,
        non_effect_allele_col=args.non_effect_allele_col,
        effect_allele_col=args.effect_allele_col,
        n_col=args.n_col,
        maf_col=args.maf_col,
        mac_col=args.mac_col,
        p_col=args.p_col,
        statistic_col=args.statistic_col,
        se_col=args.se_col,
        beta_col=args.beta_col,
        var_beta_col=args.var_beta_col,
        sdY_col=args.sdY_col,
        gene_id_col=args.gene_id_col,
        tissue_col=args.tissue_col
    )

    # Run the full pipeline
    transformer.run(
        df=df,
        expanded_ranges=expanded_ranges,
        output_dir=args.output_dir,
        input_fp=input_fp
    )


# ------------------------- Command Line Interface -------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Transform summary statistics to standardized format.")
    parser.add_argument("--input_file", type=str, required=True, help="Path to the input file containing summary statistics.")
    parser.add_argument("--expanded_ranges_file", type=str, required=True, help="Path to the JSON file containing expanded locus ranges.")
    parser.add_argument("--output_dir", type=str, required=True, help="Directory to save the output files.")
    parser.add_argument("--header_lines", type=int, default=0, help="Number of header lines to skip in the input file.")
    parser.add_argument("--entrez_email", type=str, required=True, help="Email address for NCBI Entrez API.")
    parser.add_argument("--chr_col", type=str, required=True, help="Column name for chromosome.")
    parser.add_argument("--pos_col", type=str, required=True, help="Column name for genomic position.")
    parser.add_argument("--rsid_col", type=str, default=None, help="Column name for rsID (optional).")
    parser.add_argument("--gene_id_col", type=str, required=True, help="Column name for gene ID.")
    parser.add_argument("--non_effect_allele_col", type=str, default=None, help="Column name for non-effect allele (optional).")
    parser.add_argument("--effect_allele_col", type=str, default=None, help="Column name for effect allele (optional).")
    parser.add_argument("--p_col", type=str, default=None, help="Column name for p-value (optional).")
    parser.add_argument("--beta_col", type=str, default=None, help="Column name for beta coefficient (optional).")
    parser.add_argument("--se_col", type=str, default=None, help="Column name for standard error (optional).")
    parser.add_argument("--n_col", type=str, default=None, help="Column name for sample size (optional).")
    parser.add_argument("--statistic_col", type=str, default=None, help="Column name for test statistic (optional).")
    parser.add_argument("--maf_col", type=str, default=None, help="Column name for minor allele frequency (optional).")
    parser.add_argument("--mac_col", type=str, default=None, help="Column name for minor allele count (optional).")
    parser.add_argument("--var_beta_col", type=str, default=None, help="Column name for variance of beta (optional).")
    parser.add_argument("--sdY_col", type=str, default=None, help="Column name for standard deviation of Y (optional).")
    parser.add_argument("--tissue_col", type=str, default=None, help="Column name for tissue label (optional).")

    args = parser.parse_args()
    main(args)
