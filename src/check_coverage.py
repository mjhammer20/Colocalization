# Imports 
import argparse
import re
import numpy as np # type: ignore (silences pylance warning)
import pandas as pd
from pathlib import Path
from typing import Optional
from helpers import load_input_data, normalize_chromosome, numeric_series, add_variant_id

# ------------------------- Helper Function Definitions -------------------------

def _set_without_na(x):
    """
    Convert a pandas Series or list to a set, dropping any NA values.
    
    Args:
        x (pd.Series or list): The input data to convert.
 
    Returns:
        set: A set containing the non-NA values from the input.

    """
    # Handle the case where the input is None
    if x is None:
        return set()
    
    # Convert the input to a pandas Series, drop NA values, and convert to a set of strings
    s = pd.Series(x).dropna()

    # Normalize numeric types to int string to avoid "int" vs "float" mismatches
    try:
        s = s.astype(float).astype(int).astype(str)

    # Handle any remaining non-numeric types
    except (ValueError, TypeError):
        s = s.astype(str)
    return set(s)

def _pct(n: float, d: float) -> float:
    """
    Compute the proportion n/d, returning NaN for invalid denominators.

    Args:
        n (float): The numerator.
        d (float): The denominator.

    Returns:
        float: The proportion n/d, or NaN if d is None, 0, or NaN.

    """
    # Calculate the proportion, handling cases where the denominator is invalid
    if d is None or d == 0 or pd.isna(d):
        return np.nan
    return n / d

def _summarize_overlap(s1: set, s2: set) -> tuple[int, int, int, float, float]:
    """
    Summarize the overlap between two sets, returning counts and proportions.

    Args:
        s1 (set): The first set.
        s2 (set): The second set.

    Returns:
        tuple[int, int, int, float, float]: A tuple containing:
            - n_s1 (int): Number of elements in s1.
            - n_s2 (int): Number of elements in s2.
            - n_overlap (int): Number of elements in both s1 and s2.
            - pct_s1_in_s2 (float): Proportion of s1 that overlaps with s2.
            - pct_s2_in_s1 (float): Proportion of s2 that overlaps with s1.

    """
    # Compute counts and proportions
    n_s1 = len(s1)
    n_s2 = len(s2)
    n_overlap = len(s1 & s2)
    pct_s1_in_s2 = _pct(n_overlap, n_s1)
    pct_s2_in_s1 = _pct(n_overlap, n_s2)

    return n_s1, n_s2, n_overlap, pct_s1_in_s2, pct_s2_in_s1

def _compute_coverage_metrics(snp_1: set, pos_1: set, var_1: set, source1: str, snp_2: set, pos_2: set, var_2: set, source2: str) -> dict:
    """
    Compute pairwise coverage metrics between two sources based on SNP IDs, positions, and alleles.

    Args:
        snp_1 (set): Set of SNP IDs from source1.
        pos_1 (set): Set of SNP positions from source1.
        var_1 (set): Set of SNP alleles from source1.
        source1 (str): Name label for the first source (used in metric keys).
        snp_2 (set): Set of SNP IDs from source2.
        pos_2 (set): Set of SNP positions from source2.
        var_2 (set): Set of SNP alleles from source2.
        source2 (str): Name label for the second source (used in metric keys).

    Returns:
        dict: A dictionary containing counts and proportions of overlapping SNPs
              between source1 and source2, keyed by source name.

    """
    # Compute coverage metrics for GWAS and LD based on SNP IDs
    n_snp_df1, n_snp_df2, n_snp_overlap, pct_df1_snp_overlap_df2, pct_df2_snp_overlap_df1 = _summarize_overlap(snp_1, snp_2)
    
    # Compute coverage metrics for GWAS and LD based on SNP positions
    n_pos_df1, n_pos_df2, n_pos_overlap, pct_df1_pos_overlap_df2, pct_df2_pos_overlap_df1 = _summarize_overlap(pos_1, pos_2)

    # Compute coverage metrics for GWAS and LD based on SNP alleles
    n_var_df1, n_var_df2, n_var_overlap, pct_df1_var_overlap_df2, pct_df2_var_overlap_df1 = _summarize_overlap(var_1, var_2)

    # Merge the coverage metrics into a single dictionary
    coverage_metrics = {
        f"n_snp_{source1}": n_snp_df1,
        f"n_snp_{source2}": n_snp_df2,
        f"n_snp_overlap_{source1}_{source2}": n_snp_overlap,
        f"pct_{source1}_snp_overlap_{source2}": pct_df1_snp_overlap_df2,
        f"pct_{source2}_snp_overlap_{source1}": pct_df2_snp_overlap_df1,
        f"n_pos_{source1}": n_pos_df1,
        f"n_pos_{source2}": n_pos_df2,
        f"n_pos_overlap_{source1}_{source2}": n_pos_overlap,
        f"pct_{source1}_pos_overlap_{source2}": pct_df1_pos_overlap_df2,
        f"pct_{source2}_pos_overlap_{source1}": pct_df2_pos_overlap_df1,
        f"n_var_{source1}": n_var_df1,
        f"n_var_{source2}": n_var_df2,
        f"n_var_overlap_{source1}_{source2}": n_var_overlap,
        f"pct_{source1}_var_overlap_{source2}": pct_df1_var_overlap_df2,
        f"pct_{source2}_var_overlap_{source1}": pct_df2_var_overlap_df1
    }

    return coverage_metrics

def _rank_class(overlap_class: str) -> float:
    """
    Return a numeric rank for an overlap class label (larger is better).

    Args:
        overlap_class (str): One of "Low", "Medium", "High", "No variants", or "Missing".

    Returns:
        float: Numeric rank (0 = Low, 1 = Medium, 2 = High, NaN for unrankable classes).

    """
    # Define a mapping of overlap class labels to numeric ranks
    ranks = {
        "Low": 0,
        "Medium": 1,
        "High": 2,
        "No variants": np.nan,
        "Missing": np.nan,
    }

    return ranks.get(overlap_class, np.nan)

def _determine_overall_overlap_class(ld_missing: bool, gwas_qtl_class: str, gwas_ld_class: str, qtl_ld_class: str) -> str:
    """
    Determine the overall overlap classification across GWAS, QTL, and LD.

    Args:
        ld_missing (bool): Whether the LD BIM file is missing or empty.
        gwas_qtl_class (str): Overlap class for GWAS vs. QTL.
        gwas_ld_class (str): Overlap class for GWAS vs. LD.
        qtl_ld_class (str): Overlap class for QTL vs. LD.

    Returns:
        str: Overall overlap class ("Low", "Medium", "High", or "No variants").

    """
    # If LD is missing, the overall overlap class is automatically "Low"
    if ld_missing:
                overall_overlap_class = "Low"
    
    # If LD is present, determine the overall overlap class based on the minimum rank of the individual overlap classes
    else:
        classes = [gwas_ld_class, qtl_ld_class, gwas_qtl_class]
        numeric_ranks = [_rank_class(c) for c in classes if not pd.isna(_rank_class(c))]
        
        if not numeric_ranks:
            overall_overlap_class = "No variants"
        
        else:
            overall_overlap_class = {0: "Low", 1: "Medium", 2: "High"}[int(min(numeric_ranks))]

    return overall_overlap_class


def _determine_review_priority(ld_missing: bool, overall_overlap_class: str, missing_critical_lead: bool) -> str:
    """
    Determine the manual review priority for a locus based on overlap quality.

    Args:
        ld_missing (bool): Whether the LD BIM file is missing or empty.
        overall_overlap_class (str): The overall overlap class ("Low", "Medium", "High", or "No variants").
        missing_critical_lead (bool): Whether a critical lead variant is absent from LD.

    Returns:
        str: Review priority label ("High", "Medium", or "Low").

    """
    # Assign review priority based on the overlap class and presence of critical lead variants
    if ld_missing or overall_overlap_class == "Low" or missing_critical_lead:
        review_priority = "High"
    elif overall_overlap_class == "Medium":
        review_priority = "Medium"
    elif overall_overlap_class == "High":
        review_priority = "Low"
    else:
        review_priority = "Medium"

    return review_priority

def _provide_review_priority_reasons(ld_missing: bool, gwas_ld_class: str, qtl_ld_class: str, gwas_qtl_class: str, n_qtl: int, n_gwas: int, gwas_lead_in_ld: dict, qtl_lead_in_ld: dict) -> str:
    """
    Build a human-readable string explaining the review priority for a locus.

    Args:
        ld_missing (bool): Whether the LD BIM file is missing or empty.
        gwas_ld_class (str): Overlap class for GWAS vs. LD.
        qtl_ld_class (str): Overlap class for QTL vs. LD.
        gwas_qtl_class (str): Overlap class for GWAS vs. QTL.
        n_qtl (int): Number of QTL variants in the locus.
        n_gwas (int): Number of GWAS variants in the locus.
        gwas_lead_in_ld (dict): Dictionary of booleans indicating whether the lead GWAS variant is present in LD.
        qtl_lead_in_ld (dict): Dictionary of booleans indicating whether the lead QTL variant is present in LD.

    Returns:
        str: A semicolon-separated string of reasons for the assigned review priority.

    """
    # Build a list of reasons for the review priority
    reasons = []
    if ld_missing:
        reasons.append("missing/empty LD BIM; PLINK likely excluded all variants")
    if gwas_ld_class in ["Low", "Missing"]:
        reasons.append("low/missing GWAS->LD coverage")
    if qtl_ld_class in ["Low", "Missing"]:
        reasons.append("low/missing QTL->LD coverage")
    if gwas_qtl_class == "Low":
        reasons.append("low GWAS<->QTL overlap")
    if n_gwas > 0 and not (gwas_lead_in_ld["gwas_lead_allele_in_ld"] or gwas_lead_in_ld["gwas_lead_pos_in_ld"]):
        reasons.append("lead GWAS variant missing from LD")
    if n_qtl > 0 and not (qtl_lead_in_ld["qtl_lead_allele_in_ld"] or qtl_lead_in_ld["qtl_lead_pos_in_ld"]):
        reasons.append("lead QTL variant missing from LD")
    if n_qtl == 0:
        reasons.append("no QTL variants in locus/tissue")
    if not reasons:
        reasons.append("coverage looks acceptable")

    return "; ".join(reasons)


# ------------------------------ Class Definition -------------------------------

class CoverageCheck:
    """
    Class to check coverage of GWAS and QTL variants in LD regions.
    """
    def __init__(self, args: argparse.Namespace):
        """
        Initialize the CoverageChecker with paths, column keys, and overlap thresholds.

        Args:
            args (argparse.Namespace): Parsed command-line arguments containing paths,
                column name keys, and overlap classification thresholds.

        """
        self.out_ld_dir = args.out_ld_dir
        self.ld_manifest = load_input_data(Path(self.out_ld_dir, args.ld_manifest), 0)
        self.standardized_locus_id_key = args.standardized_locus_id_key
        self.standardized_left_bound_key = args.standardized_left_bound_key
        self.standardized_right_bound_key = args.standardized_right_bound_key
        self.standardized_chr_key = args.standardized_chr_key
        self.standardized_snp_key = args.standardized_snp_key
        self.standardized_var_key = args.standardized_var_key
        self.standardized_pos_key = args.standardized_pos_key
        self.standardized_non_effect_allele_key = args.standardized_non_effect_allele_key
        self.standardized_effect_allele_key = args.standardized_effect_allele_key
        self.standardized_p_key = args.standardized_p_key
        self.manifest_bim_key = args.manifest_bim_key
        self.bim_cm_key = "cm"
        self.bim_missing_key = "missing_bim"
        self.bim_source_path_key = "source_path"
        self.gwas = load_input_data(Path(args.gwas_fp), 0)
        self.gwas_strata_key = args.gwas_strata_key
        self.qtl = load_input_data(Path(args.qtl_fp), 0)
        self.qtl_strata_key = args.qtl_strata_key
        self.high_overlap_min = args.high_overlap_min
        self.medium_overlap_min = args.medium_overlap_min
        self.results_ld_missing_key = "ld_missing"
        self.results_n_var_gwas_key = "n_var_gwas"
        self.results_n_var_qtl_key = "n_var_qtl"
        self.results_pct_gwas_var_overlap_ld_key = "pct_gwas_var_overlap_ld"
        self.results_pct_qtl_var_overlap_ld_key = "pct_qtl_var_overlap_ld"
        self.results_pct_gwas_var_overlap_qtl_key = "pct_gwas_var_overlap_qtl"
        self.results_pct_qtl_var_overlap_gwas_key = "pct_qtl_var_overlap_gwas"
        self.results_review_priority_key = "review_priority"
        self.results_overall_overlap_class_key = "overall_overlap_class"
        self.results_review_priority_reason_key = "review_priority_reason"
        self.output_priority_rank_key = "priority_rank"
        self.output_overlap_rank_key = "overlap_rank"
        self.out_qc_dir = args.out_qc_dir


    def load_bim(self, bim_path: Path) -> pd.DataFrame:
        """
        Load a PLINK .bim file, returning an empty DataFrame with metadata if the file is missing or empty.

        Args:
            bim_path (Path): Path to the .bim file to load.

        Returns:
            pd.DataFrame: A DataFrame containing the .bim file contents with additional
                derived columns for locus (chr, pos) and allele (chr, pos, a1, a2) tuples.
                The DataFrame"s attrs dict contains "missing_bim" (bool) and "source_path" (str).

        """
        # Check if the .bim file exists; if not, return an empty DataFrame with metadata
        if not bim_path.exists():
            bim = pd.DataFrame(
                columns=[
                    self.standardized_chr_key,
                    self.standardized_snp_key,
                    self.bim_cm_key,
                    self.standardized_pos_key,
                    self.standardized_non_effect_allele_key,
                    self.standardized_effect_allele_key,
                    self.standardized_locus_id_key,
                    self.standardized_var_key
                ]
            )
            bim.attrs[self.bim_missing_key] = True
            bim.attrs[self.bim_source_path_key] = str(bim_path)
            return bim

        # If the .bim file exists, attempt to load it
        try:
            bim = pd.read_csv(
                bim_path,
                sep=r"\s+",
                header=None,
                names = [
                    self.standardized_chr_key,
                    self.standardized_snp_key,
                    self.bim_cm_key,
                    self.standardized_pos_key,
                    self.standardized_non_effect_allele_key,
                    self.standardized_effect_allele_key,
                ]
            )

        # Handle the case where the .bim file is empty by returning an empty DataFrame with metadata
        except pd.errors.EmptyDataError:
            bim = pd.DataFrame(
                columns=[
                    self.standardized_chr_key,
                    self.standardized_snp_key,
                    self.bim_cm_key,
                    self.standardized_pos_key,
                    self.standardized_non_effect_allele_key,
                    self.standardized_effect_allele_key,
                    self.standardized_locus_id_key,
                    self.standardized_var_key
                ]
            )
            bim.attrs[self.bim_missing_key] = True
            bim.attrs[self.bim_source_path_key] = str(bim_path)
            return bim
        
        # Normalize chromosome representation
        bim[self.standardized_chr_key] = bim[self.standardized_chr_key].astype(str).map(normalize_chromosome)

        # Ensure SNP ID column is a string
        bim[self.standardized_snp_key] = bim[self.standardized_snp_key].astype(str)

        # Ensure postion colum in numeric
        bim[self.standardized_pos_key] = numeric_series(bim, self.standardized_pos_key)
        
        # Create derived column for variant ID (chr:pos:a1:a2) using the add_variant_id helper function
        bim = add_variant_id(bim, self.standardized_chr_key, self.standardized_pos_key, self.standardized_non_effect_allele_key, self.standardized_effect_allele_key, self.standardized_var_key)

        # Set metadata attributes indicating that the .bim file was successfully loaded
        bim.attrs[self.bim_missing_key] = False
        bim.attrs[self.bim_source_path_key] = str(bim_path)

        return bim


    def get_locus_subsets(self, locus_id: str) -> tuple[dict, dict]:
        """
        Extract per-stratum GWAS and QTL subsets for a given locus.

        Args:
            locus_id (str): The locus identifier used to look up boundaries in the LD manifest.

        Returns:
            tuple[dict, dict]: A tuple of (gwas_subsets, qtl_subsets), where each is a dict
                mapping stratum label -> subset DataFrame restricted to the locus coordinates.

        """
        # Extract chromosome and position boundaries for the locus
        chrom = str(self.ld_manifest.loc[self.ld_manifest[self.standardized_locus_id_key] == locus_id, self.standardized_chr_key].values[0])
        left_bound = int(self.ld_manifest.loc[self.ld_manifest[self.standardized_locus_id_key] == locus_id, self.standardized_left_bound_key].values[0])
        right_bound = int(self.ld_manifest.loc[self.ld_manifest[self.standardized_locus_id_key] == locus_id, self.standardized_right_bound_key].values[0])

        # Ensure position columns in GWAS and QTL DataFrames are numeric for filtering
        self.gwas[self.standardized_pos_key] = numeric_series(self.gwas, self.standardized_pos_key)
        self.qtl[self.standardized_pos_key] = numeric_series(self.qtl, self.standardized_pos_key)

        # Filter the GWAS DataFrame to include only variants within the locus boundaries
        gwas_filtered = self.gwas[
            (self.gwas[self.standardized_chr_key].astype(str) == chrom) &
            (self.gwas[self.standardized_pos_key] >= left_bound) &
            (self.gwas[self.standardized_pos_key] <= right_bound)
        ]

        # Create subsets of the GWAS DataFrame based on the locus boundaries and stratification keys
        if self.gwas_strata_key:
            gwas_subsets = {key: sub for key, sub in gwas_filtered.groupby(self.gwas_strata_key)}
        else:
            gwas_subsets = {"Full": gwas_filtered}

        # Filter the QTL DataFrame to include only variants within the locus boundaries
        qtl_filtered = self.qtl[
            (self.qtl[self.standardized_chr_key].astype(str) == chrom) &
            (self.qtl[self.standardized_pos_key] >= left_bound) &
            (self.qtl[self.standardized_pos_key] <= right_bound)
        ]

        # Create subsets of the QTL DataFrame based on the locus boundaries and stratification keys
        if self.qtl_strata_key:
            qtl_subsets = {key: sub for key, sub in qtl_filtered.groupby(self.qtl_strata_key)}

        else:
            qtl_subsets = {"Bulk": qtl_filtered}

        return gwas_subsets, qtl_subsets


    def get_lead_variant_info(self, df: pd.DataFrame) -> dict:
        """
        Extract key fields for the lead (most significant) variant in a DataFrame.

        Args:
            df (pd.DataFrame): A DataFrame of variants for a single locus/stratum,
                containing standardized ID, position, allele, and p-value columns.

        Returns:
            dict: A dictionary with keys for standardized ID, position, allele, and p-value
                of the lead variant, or NaN values if the DataFrame is empty.

        """
        # Sort the DataFrame by p-value and select the lead variant (lowest p-value)
        if len(df) == 0:
            return {self.standardized_snp_key: np.nan, self.standardized_pos_key: np.nan, self.standardized_var_key: np.nan, self.standardized_p_key: np.nan}
        lead = df.sort_values(self.standardized_p_key).iloc[0]
        return {
            self.standardized_snp_key: lead[self.standardized_snp_key],
            self.standardized_pos_key: lead.get(self.standardized_pos_key, np.nan),
            self.standardized_var_key: lead.get(self.standardized_var_key, np.nan),
            self.standardized_p_key: lead[self.standardized_p_key],
        }

    def evaluate_lead_overlap(self, lead: dict, lead_source: str, target_ids: set, target_pos: set, target_allele: set, target_source: str) -> dict:
        """
        Evaluates the overlap of a lead variant with a target set of variants based on SNP IDs, positions, and alleles.

        Args:
            lead (dict): A dictionary containing lead variant information with keys "id", "pos", and "allele".
            ids (set): A set of SNP IDs to compare against.
            pos (set): A set of SNP positions to compare against.
            allele (set): A set of SNP alleles to compare against.

        Returns:
            dict: A dictionary containing boolean values indicating whether the lead variant overlaps with the provided sets.
        """
        # Evaluate whether the lead variant's SNP ID, position, and allele are present in the target sets
        return {
            f"{lead_source}_lead_id_in_{target_source}": lead.get(self.standardized_snp_key) in target_ids,
            f"{lead_source}_lead_pos_in_{target_source}": lead.get(self.standardized_pos_key) in target_pos,
            f"{lead_source}_lead_allele_in_{target_source}": lead.get(self.standardized_var_key) in target_allele
        }


    def classify_overlap(self, value: float, n_query: Optional[int] = None) -> str:
        """
        Classify a coverage proportion into a discrete overlap class label.

        Args:
            value (float): The coverage proportion to classify (e.g. pct overlap).
            n_query (int, optional): The number of query variants. If 0, returns
                "No variants" regardless of value. Defaults to None.

        Returns:
            str: Overlap class label — one of "No variants", "Missing", "Low",
                "Medium", or "High".

        """
        # Classify the overlap based on the provided value and number of query variants
        if n_query is not None and n_query == 0:
            return "No variants"
        if pd.isna(value):
            return "Missing"
        if value >= self.high_overlap_min:
            return "High"
        if value >= self.medium_overlap_min:
            return "Medium"
        return "Low"


    def evaluate_coverage(self) -> pd.DataFrame:
        """
        Compute variant coverage metrics for every locus x GWAS stratum x QTL stratum combination.

        For each locus in the LD manifest, loads the corresponding PLINK .bim file and
        computes pairwise overlap metrics between GWAS, QTL, and LD variant sets. Classifies
        each comparison as "Low", "Medium", "High", "No variants", or "Missing", determines
        an overall overlap class and review priority, and records human-readable reasons for
        the assigned priority.

        Returns:
            pd.DataFrame: A DataFrame with one row per locus x GWAS stratum x QTL stratum
                combination, containing coverage counts, overlap proportions, overlap class
                labels, lead variant overlap flags, overall overlap class, review priority,
                and priority reasons.

        """
        # Initialize list to store coverage rows
        coverage_rows = []

        # Iterate over each locus in the LD manifest
        for _, locus_row in self.ld_manifest.iterrows():
            
            # Extract the locus ID from the current row using the specified key
            locus_id = locus_row[self.standardized_locus_id_key]

            # Load the corresponding .bim file for the locus
            bim_path = Path(self.out_ld_dir, locus_row[self.manifest_bim_key]) #type: ignore silences pylance warning
            bim = self.load_bim(bim_path)

            # Check if LD is missing based on the existence of the .bim file and its content
            ld_missing = bool(bim.attrs.get(self.bim_missing_key, False)) or len(bim) == 0

            # Extract SNP IDs, positions, and alleles from the .bim file if not ld_missing, otherwise set them to empty sets
            ld_snp = _set_without_na(bim[self.standardized_snp_key]) if not ld_missing else set()
            ld_pos = _set_without_na(bim[self.standardized_pos_key]) if not ld_missing else set()
            ld_var = _set_without_na(bim[self.standardized_var_key]) if not ld_missing else set()

            # Extract locus specific subset of GWAS and QTL data
            gwas_subsets, qtl_subsets = self.get_locus_subsets(locus_id) #type: ignore silences pylance warning

            # Initialize dict to store ld coverage metrics for each stratum x locus combination
            qtl_ld_coverage_metrics = {}

            # Intialize dict to store lead qtl overlap with ld evaluation for each stratum x locus combination
            qtl_lead_in_ld = {}

            # Initalize dict to store qtl overlap classification for each stratum x locus combination
            qtl_ld_class = {}

            # Iterate over each GWAS stratum and compute coverage separately
            for key, gwas_sub in gwas_subsets.items():
                
                # Extract the stratum name from the key tuple
                gwas_stratum = key

                # Extract lead variant information for GWAS
                gwas_lead = self.get_lead_variant_info(gwas_sub)

                # Extract SNP IDs, positions, and alleles for the GWAS subset if it exists, otherwise set them to empty sets
                gwas_snp= _set_without_na(gwas_sub[self.standardized_snp_key]) if len(gwas_sub) else set()
                gwas_pos = _set_without_na(gwas_sub[self.standardized_pos_key]) if len(gwas_sub) else set()
                gwas_var = _set_without_na(gwas_sub[self.standardized_var_key]) if len(gwas_sub) else set()

                # Compute coverage metrics for GWAS stratum and LD (Will be unique for each stratum x locus combination, so no need to cache and check for existence). 
                gwas_ld_coverage_metrics = _compute_coverage_metrics(
                    snp_1 = gwas_snp,
                    pos_1 = gwas_pos,
                    var_1 = gwas_var,
                    source1 = "gwas",
                    snp_2 = ld_snp,
                    pos_2 = ld_pos,
                    var_2 = ld_var,
                    source2 = "ld"
                )

                # Evaluate lead variant overlap for GWAS and LD
                gwas_lead_in_ld = self.evaluate_lead_overlap(
                    lead = gwas_lead,
                    lead_source = "gwas",
                    target_ids = ld_snp,
                    target_pos = ld_pos,
                    target_allele = ld_var,
                    target_source = "ld"
                )

                # Classify the overlap between GWAS and LD based on the computed coverage metrics. If LD is missing, classify as "Missing", otherwise use the classify_overlap method to determine the overlap class based on the percentage of GWAS variant overlap with LD and the number of GWAS variants.
                gwas_ld_class = "Missing" if ld_missing else self.classify_overlap(
                    value = gwas_ld_coverage_metrics[self.results_pct_gwas_var_overlap_ld_key],
                    n_query = gwas_ld_coverage_metrics[self.results_n_var_gwas_key]
                )

                # Iterate over each QTL stratum and compute coverage separately
                for key, qtl_sub in qtl_subsets.items():

                    # Extract the stratum name from the key tuple
                    qtl_stratum = key

                    # Extract lead variant information for QTL
                    qtl_lead = self.get_lead_variant_info(qtl_sub)

                    # Extract SNP IDs, positions, and alleles for the QTL subset if it exists, otherwise set them to empty sets
                    qtl_snp = _set_without_na(qtl_sub[self.standardized_snp_key]) if len(qtl_sub) else set()
                    qtl_pos = _set_without_na(qtl_sub[self.standardized_pos_key]) if len(qtl_sub) else set()
                    qtl_var = _set_without_na(qtl_sub[self.standardized_var_key]) if len(qtl_sub) else set()

                    # Compute coverage metrics for QTL stratum and LD. Cache the results to avoid redundant computations for the same stratum x locus combination.
                    if not qtl_stratum in qtl_ld_coverage_metrics.keys():
                        qtl_ld_coverage_metrics[qtl_stratum] = _compute_coverage_metrics(
                            snp_1 = qtl_snp,
                            pos_1 = qtl_pos,
                            var_1 = qtl_var,
                            source1 = "qtl",
                            snp_2 = ld_snp,
                            pos_2 = ld_pos,
                            var_2 = ld_var,
                            source2 = "ld"
                        )

                    # Evaluate lead variant overlap for QTL and LD. Cache the results to avoid redundant computations for the same stratum x locus combination.
                    if not qtl_stratum in qtl_lead_in_ld.keys():
                        qtl_lead_in_ld[qtl_stratum] = self.evaluate_lead_overlap(
                            lead = qtl_lead,
                            lead_source = "qtl",
                            target_ids = ld_snp,
                            target_pos = ld_pos,
                            target_allele = ld_var,
                            target_source = "ld"
                        )

                    # Classify the overlap between QTL and LD based on the computed coverage metrics. If LD is missing, classify as "Missing", otherwise use the classify_overlap method to determine the overlap class based on the percentage of QTL variant overlap with LD and the number of QTL variants.
                    if not qtl_stratum in qtl_ld_class.keys():
                        qtl_ld_class[qtl_stratum] = "Missing" if ld_missing else self.classify_overlap(
                            value = qtl_ld_coverage_metrics[qtl_stratum][self.results_pct_qtl_var_overlap_ld_key],
                            n_query = qtl_ld_coverage_metrics[qtl_stratum][self.results_n_var_qtl_key]
                        )

                    # Compute coverage metrics for GWAS stratum and QTL stratum (Will be unique for each stratum x stratum combination, so no need to cache and check for existence)
                    gwas_qtl_coverage_metrics = _compute_coverage_metrics(
                        snp_1 = gwas_snp,
                        pos_1 = gwas_pos,
                        var_1 = gwas_var,
                        source1 = "gwas",
                        snp_2 = qtl_snp,
                        pos_2 = qtl_pos,
                        var_2 = qtl_var,
                        source2 = "qtl"
                    )

                    # Evaluate lead variant overlap for GWAS and QTL
                    gwas_lead_in_qtl = self.evaluate_lead_overlap(
                        lead = gwas_lead,
                        lead_source = "gwas",
                        target_ids = qtl_snp,
                        target_pos = qtl_pos,
                        target_allele = qtl_var,
                        target_source = "qtl"
                    )
                    qtl_lead_in_gwas = self.evaluate_lead_overlap(
                        lead = qtl_lead,
                        lead_source = "qtl",
                        target_ids = gwas_snp,
                        target_pos = gwas_pos,
                        target_allele = gwas_var,
                        target_source = "gwas"
                    )

                    # Classify the overlap between GWAS and QTL based on the computed coverage metrics. 
                    # If either GWAS or QTL is missing, classify as "Missing", otherwise use the classify_overlap method to determine the overlap class based on the percentage of GWAS variant overlap with QTL and the number of GWAS variants.
                    gwas_qtl_pct_min = np.nanmin([gwas_qtl_coverage_metrics[self.results_pct_gwas_var_overlap_qtl_key], gwas_qtl_coverage_metrics[self.results_pct_qtl_var_overlap_gwas_key]]) if gwas_qtl_coverage_metrics[self.results_n_var_gwas_key] and gwas_qtl_coverage_metrics[self.results_n_var_qtl_key] else np.nan
                    gwas_qtl_n_min = min(gwas_qtl_coverage_metrics[self.results_n_var_gwas_key], gwas_qtl_coverage_metrics[self.results_n_var_qtl_key]) if gwas_qtl_coverage_metrics[self.results_n_var_gwas_key] and gwas_qtl_coverage_metrics[self.results_n_var_qtl_key] else 0
                    gwas_qtl_class = self.classify_overlap(
                        value=gwas_qtl_pct_min,
                        n_query=gwas_qtl_n_min
                    )

                    # Determine the overall overlap class based on the individual overlap classes for GWAS-LD, QTL-LD, and GWAS-QTL. 
                    # If LD is missing, classify as "Low", otherwise use the _determine_overall_overlap_class method to determine the overall overlap class based on the individual overlap classes.
                    overall_overlap_class = _determine_overall_overlap_class(
                        ld_missing = ld_missing,
                        gwas_qtl_class = gwas_qtl_class,
                        gwas_ld_class = gwas_ld_class,
                        qtl_ld_class = qtl_ld_class[qtl_stratum]
                    )

                    # Determine the review priority based on the overall overlap class and whether any critical lead variants are missing. 
                    # If LD is missing or if there are GWAS or QTL variants but their lead variants are not in LD, set the review priority to "High", otherwise set it to "Low".
                    missing_critical_lead = (
                        ld_missing or
                        (gwas_ld_coverage_metrics[self.results_n_var_gwas_key] > 0 and not (gwas_lead_in_ld["gwas_lead_allele_in_ld"] or gwas_lead_in_ld["gwas_lead_pos_in_ld"])) or
                        (qtl_ld_coverage_metrics[qtl_stratum][self.results_n_var_qtl_key] > 0 and not (qtl_lead_in_ld[qtl_stratum]["qtl_lead_allele_in_ld"] or qtl_lead_in_ld[qtl_stratum]["qtl_lead_pos_in_ld"]))
                    )
                    review_priority = _determine_review_priority(
                        ld_missing = ld_missing,
                        overall_overlap_class = overall_overlap_class,
                        missing_critical_lead = missing_critical_lead
                    )

                    # Provide reasons for the review priority based on the overlap classes and whether any critical lead variants are missing. 
                    # If LD is missing, add a reason indicating that the LD BIM file is missing or empty. If the GWAS-LD or QTL-LD overlap classes are "Low" or "Missing", add reasons indicating low or missing coverage. 
                    # If the GWAS-QTL overlap class is "Low", add a reason indicating low overlap. If there are GWAS or QTL variants but their lead variants are not in LD, add reasons indicating that the lead variants are missing from LD. 
                    # If there are no QTL variants in the locus/tissue, add a reason indicating that. If none of these conditions apply, add a reason indicating that coverage looks acceptable.
                    review_priority_reasons = _provide_review_priority_reasons(
                        ld_missing = ld_missing,
                        gwas_ld_class = gwas_ld_class,
                        qtl_ld_class = qtl_ld_class[qtl_stratum],
                        gwas_qtl_class = gwas_qtl_class,
                        n_qtl = qtl_ld_coverage_metrics[qtl_stratum][self.results_n_var_qtl_key],
                        n_gwas = gwas_ld_coverage_metrics[self.results_n_var_gwas_key],
                        gwas_lead_in_ld = gwas_lead_in_ld,
                        qtl_lead_in_ld = qtl_lead_in_ld[qtl_stratum]
                    )

                    # Append the results to the coverage rows
                    coverage_rows.append({
                        self.standardized_locus_id_key: locus_id,
                        self.standardized_chr_key: locus_row[self.standardized_chr_key],
                        self.standardized_left_bound_key: locus_row[self.standardized_left_bound_key],
                        self.standardized_right_bound_key: locus_row[self.standardized_right_bound_key],
                        self.results_ld_missing_key: ld_missing,
                        self.manifest_bim_key: locus_row[self.manifest_bim_key],
                        f"{self.manifest_bim_key}_path": str(bim_path),
                        "gwas_stratum": gwas_stratum,
                        "qtl_stratum": qtl_stratum,
                        **gwas_ld_coverage_metrics,
                        "gwas_ld_overlap_class": gwas_ld_class,
                        "lead_gwas_variant": gwas_lead[self.standardized_snp_key],
                        "lead_gwas_p": gwas_lead[self.standardized_p_key],
                        **gwas_lead_in_ld,
                        **qtl_ld_coverage_metrics[qtl_stratum],
                        "qtl_ld_overlap_class": qtl_ld_class[qtl_stratum],
                        "lead_qtl_variant": qtl_lead[self.standardized_snp_key],
                        "lead_qtl_p": qtl_lead[self.standardized_p_key],
                        **qtl_lead_in_ld[qtl_stratum],
                        **gwas_qtl_coverage_metrics,
                        "gwas_qtl_overlap_class": gwas_qtl_class,
                        **gwas_lead_in_qtl,
                        **qtl_lead_in_gwas,
                        self.results_overall_overlap_class_key: overall_overlap_class,
                        self.results_review_priority_key: review_priority,
                        self.results_review_priority_reason_key: review_priority_reasons
                    })

        # Return the final coverage DataFrame
        return pd.DataFrame(coverage_rows)


    def save_results(self, results_df: pd.DataFrame) -> None:
        # Sort so worst / most actionable rows appear first.
        priority_order = {"High": 0, "Medium": 1, "Low": 2}
        overlap_order = {"Low": 0, "Medium": 1, "High": 2, "No variants": 3, "Missing": 4}
        results_df["_priority_sort"] = results_df[self.results_review_priority_key].map(priority_order).fillna(9)
        results_df["_overlap_sort"] = results_df[self.results_overall_overlap_class_key].map(overlap_order).fillna(9)
        sort_cols = [
            "_priority_sort",
            "_overlap_sort",    
            self.results_ld_missing_key,
            self.results_pct_gwas_var_overlap_ld_key,
            self.results_pct_qtl_var_overlap_ld_key,
            self.results_pct_gwas_var_overlap_qtl_key,
        ]
        coverage_sorted = results_df.sort_values(sort_cols, ascending=[True, True, False, True, True, True], na_position="first")
        coverage_sorted = coverage_sorted.drop(columns=["_priority_sort", "_overlap_sort"])

        # Save the full coverage results to a TSV file in the output QC directory.
        all_out_fp = Path(self.out_qc_dir, "coverage_by_locus.tsv")
        coverage_sorted.to_csv(all_out_fp, sep="\t", index=False)
        print(f"Saved coverage results for {len(coverage_sorted)} locus x stratum combinations to {all_out_fp}")

        # Save a filtered file containing only loci with "High" or "Medium" review priority to a TSV file in the output QC directory.
        problem = coverage_sorted[coverage_sorted[self.results_review_priority_key].isin(["High", "Medium"])].copy()
        problem_out_fp = Path(self.out_qc_dir, "problem_loci_ranked.tsv")
        problem.to_csv(problem_out_fp, sep="\t", index=False)
        print(f"Saved problem loci results for {len(problem)} locus x stratum combinations to {problem_out_fp}")

        # Collapsed locus-level summary: one row per locus with worst priority/overlap observed across strata.
        tmp = results_df.copy()
        tmp[self.output_priority_rank_key] = tmp[self.results_review_priority_key].map(priority_order)
        tmp[self.output_overlap_rank_key] = tmp[self.results_overall_overlap_class_key].map(overlap_order)
        locus_summary = (
            tmp
            .groupby([self.standardized_locus_id_key, self.standardized_chr_key, self.standardized_left_bound_key, self.standardized_right_bound_key], as_index=False)
            .agg(
                any_ld_missing=(self.results_ld_missing_key, "max"),
                worst_review_priority_rank=(self.output_priority_rank_key, "min"),
                worst_overlap_rank=(self.output_overlap_rank_key, "min"),
                n_rows=(self.standardized_locus_id_key, "size"),
                n_high_priority_rows=(self.results_review_priority_key, lambda x: int((x == "High").sum())),
                n_medium_priority_rows=(self.results_review_priority_key, lambda x: int((x == "Medium").sum())),
                min_pct_gwas_in_ld_by_allele=(self.results_pct_gwas_var_overlap_ld_key, "min"),
                min_pct_qtl_in_ld_by_allele=(self.results_pct_qtl_var_overlap_ld_key, "min"),
                min_pct_gwas_in_qtl_by_allele=(self.results_pct_gwas_var_overlap_qtl_key, "min"),
            )
        )

        # Map the worst priority and overlap ranks back to their corresponding string labels for easier interpretation in the summary.
        reverse_priority = {v: k for k, v in priority_order.items()}
        reverse_overlap = {v: k for k, v in overlap_order.items()}
        locus_summary["worst_review_priority"] = locus_summary["worst_review_priority_rank"].map(reverse_priority)
        locus_summary["worst_overlap_class"] = locus_summary["worst_overlap_rank"].map(reverse_overlap)
        locus_summary = locus_summary.sort_values(
            [
                "worst_review_priority_rank",
                "any_ld_missing",
                "worst_overlap_rank",
                "min_pct_gwas_in_ld_by_allele",
                "min_pct_qtl_in_ld_by_allele",
            ],
            ascending=[True, False, True, True, True],
            na_position="first"
        )

        # Save the locus-level summary to a TSV file in the output QC directory.
        locus_out_fp = Path(self.out_qc_dir, "coverage_summary_by_locus.tsv")
        locus_summary.to_csv(locus_out_fp, sep="\t", index=False)
        print(f"Saved locus summary results for {len(locus_summary)} locus x stratum combinations to {locus_out_fp}")

        # Save a filtered file containing only loci with missing or empty LD to a TSV file in the output QC directory.
        missing_ld = (
            results_df.loc[results_df[self.results_ld_missing_key], [
                self.standardized_locus_id_key, self.standardized_chr_key, self.standardized_left_bound_key, self.standardized_right_bound_key,
                self.manifest_bim_key, f"{self.manifest_bim_key}_path", self.results_review_priority_reason_key
            ]]
            .drop_duplicates()
        )
        missing_ld_out_fp = Path(self.out_qc_dir, "missing_or_empty_ld_loci.tsv")
        missing_ld.to_csv(missing_ld_out_fp, sep="\t", index=False)
        print(f"Saved missing/empty LD loci results for {len(missing_ld)} loci to {missing_ld_out_fp}")


    def run(self) -> None:
        """
        Run the coverage evaluation process and save the results.

        Returns:
            None: The function saves the coverage evaluation results to a CSV file specified in self.out_qc_dir.

        """
        # Evaluate coverage
        results = self.evaluate_coverage()

        # Save results to CSV
        self.save_results(results)
    

# ------------------------------ Main Function -------------------------------

def main(args: argparse.Namespace) -> None:
    """
    Main function to execute coverage evaluation process.

    Args:
        args (argparse.Namespace): Parsed command-line arguments containing paths to input files,
            output directory, and other parameters.

    Returns:
        None: The function saves the coverage evaluation results to a CSV file specified in args.output_file

    """
    # Create output QC directory if it doesn't exist
    qc_dir = Path(args.out_qc_dir)
    qc_dir.mkdir(parents=True, exist_ok=True)

    # Create CoverageCheck instance
    coverage_checker = CoverageCheck(args)

    # Run Coverage Evaluation
    coverage_checker.run()


# ------------------------------ Command-Line Interface -------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Evaluate variant coverage across GWAS, QTL, and LD datasets.")
    parser.add_argument("--out_ld_dir", required=True, help="Directory containing LD .bim files and the LD manifest.")
    parser.add_argument("--ld_manifest", required=True, help="Path to the LD manifest file.")
    parser.add_argument("--standardized_locus_id_key", default="LOCUS_ID", help="Column name in the LD manifest that contains locus identifiers.")
    parser.add_argument("--standardized_chr_key", default="CHR", help="Column name in the LD manifest that contains chromosome information.")
    parser.add_argument("--standardized_left_bound_key", default="LEFT_500KB", help="Column name in the LD manifest that contains left boundary positions.")
    parser.add_argument("--standardized_right_bound_key", default="RIGHT_500KB", help="Column name in the LD manifest that contains right boundary positions.")
    parser.add_argument("--standardized_snp_key", required=True, help="Column name in the GWAS/QTL files that contains standardized rsIDs.")
    parser.add_argument("--standardized_pos_key", required=True, help="Column name in the GWAS/QTL files that contains standardized position information.")
    parser.add_argument("--standardized_var_key", required=True, help="Column name in the GWAS/QTL files that contains standardized variant IDs")
    parser.add_argument("--standardized_p_key", required=True, help="Column name in the GWAS/QTL files that contains standardized p-value information.")
    parser.add_argument("--standardized_non_effect_allele_key", required=True, help="Column name in the GWAS/QTL files that contains standardized non-effect allele information.")
    parser.add_argument("--standardized_effect_allele_key", required=True, help="Column name in the GWAS/QTL files that contains standardized effect allele information.")
    parser.add_argument("--manifest_bim_key", default="BIM", help="Column name in the LD manifest that contains .bim file names.")
    parser.add_argument("--gwas_fp", required=True, help="Path to the GWAS summary statistics file.")
    parser.add_argument("--gwas_strata_key", default=None, help="Column name in the GWAS file that contains stratum identifiers (optional).")
    parser.add_argument("--qtl_fp", required=True, help="Path to the QTL summary statistics file.")
    parser.add_argument("--qtl_strata_key", default=None, help="Column name in the QTL file that contains stratum identifiers (optional).")
    parser.add_argument("--high_overlap_min", type=float, default=0.90, help="Minimum coverage proportion to classify as 'High' overlap.")
    parser.add_argument("--medium_overlap_min", type=float, default=0.70, help="Minimum coverage proportion to classify as 'Medium' overlap.")
    parser.add_argument("--out_qc_dir", required=True, help="Path to the output directory where coverage evaluation results will be saved.")

    args = parser.parse_args()
    main(args)