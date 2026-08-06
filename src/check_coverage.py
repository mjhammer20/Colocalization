# Imports 
import argparse
import re
import numpy as np # type: ignore (silences pylance warning)
import pandas as pd
from pathlib import Path
from helpers import load_input_data

# ------------------------- Helper Function Definitions -------------------------

def _first_existing(df, candidates, required=True):
    match = next((c for c in candidates if c in df.columns), None)
    if match is None and required:
        raise ValueError(f"None of these columns found: {candidates}\nAvailable columns: {list(df.columns)}")
    return match


def _set_without_na(x):
    if x is None:
        return set()
    return set(pd.Series(x).dropna().astype(str))

def _pct(n, d):
    if d is None or d == 0 or pd.isna(d):
        return np.nan
    return n / d

def _class_rank(label):
    """
    Larger is better. Used to define the worst class.
    """
    ranks = {
        "Low": 0,
        "Medium": 1,
        "High": 2,
        "No variants": np.nan,
        "Missing": np.nan,
    }
    return ranks.get(label, np.nan)


# ------------------------- Class Definition -------------------------

class CoverageChecker:
    def __init__(self, args: argparse.Namespace):
        self.out_ld_dir = args.out_ld_dir
        self.ld_manifest = args.ld_manifest
        self.ld_manifest_path = Path(self.out_ld_dir) / self.ld_manifest
        self.loc_key = args.loc_key
        self.allele_key = args.allele_key
        self.snp_key = args.snp_key
        self.chr_key = args.chr_key
        self.pos_key = args.pos_key
        self.high_overlap_min = args.high_overlap_min
        self.medium_overlap_min = args.medium_overlap_min
        self.left_bound_col = args.left_bound_col
        self.right_bound_col = args.right_bound_col
        self.locus_id_col = args.locus_id_col

    def resolve_manifest_path(self, x):
        """
        Resolve paths in ld_manifest.tsv.
        Handles absolute paths, relative paths, and basenames relative to OUT_LD_DIR or LD_MANIFEST.parent.
        Returns None if the file does not exist.
        """
        p = Path(str(x))
        candidates = [
            p,
            Path(self.out_ld_dir) / p,
            self.ld_manifest_path.parent / p,
            Path(self.out_ld_dir) / p.name,
            self.ld_manifest_path.parent / p.name,
        ]
        for cand in candidates:
            if cand.exists():
                return cand
        return None


    def read_bim(self, bim_path, missing_ok=True):
        """
        Read PLINK .bim and add matching keys.
        If the .bim is missing/empty, return an empty BIM-like dataframe with attrs['missing_bim']=True.
        """
        resolved = self.resolve_manifest_path(bim_path)

        if resolved is None:
            if not missing_ok:
                raise FileNotFoundError(f"Could not find BIM from manifest entry: {bim_path}")
            bim = pd.DataFrame(columns=["chr", "snp", "cm", "bp", "a1", "a2", "pos_key", "allele_key"])
            bim.attrs["missing_bim"] = True
            bim.attrs["source_path"] = str(bim_path)
            return bim

        try:
            bim = pd.read_csv(
                resolved,
                sep=r"\s+",
                header=None,
                names=["chr", "snp", "cm", "bp", "a1", "a2"],
            )
        except pd.errors.EmptyDataError:
            bim = pd.DataFrame(columns=["chr", "snp", "cm", "bp", "a1", "a2"])

        if len(bim) == 0:
            bim[self.loc_key] = []
            bim[self.allele_key] = []
            bim.attrs["missing_bim"] = True
            bim.attrs["source_path"] = str(resolved)
            return bim

        bim[self.snp_key] = bim["snp"].astype(str)
        bim[self.pos_key] = pd.to_numeric(bim["bp"], errors="coerce")
        bim[self.loc_key] = [self.pos_key(c, p) for c, p in zip(bim["chr"], bim["bp"])]
        bim[self.allele_key] = [
            self.allele_key(c, p, a1, a2)
            for c, p, a1, a2 in zip(bim["chr"], bim["bp"], bim["a1"], bim["a2"])
        ]
        bim.attrs["missing_bim"] = False
        bim.attrs["source_path"] = str(resolved)
        return bim

    def classify_overlap(self, value, n_query=None):
        """
        High/Medium/Low label for coverage proportion.
        If there are no query variants, return No variants.
        """
        if n_query is not None and n_query == 0:
            return "No variants"
        if pd.isna(value):
            return "Missing"
        if value >= self.high_overlap_min:
            return "High"
        if value >= self.medium_overlap_min:
            return "Medium"
        return "Low"
    

    def load_bims(self):
        """
        Load BIMs for all loci in the manifest.
        Returns a dictionary mapping locusID to BIM DataFrame.
        """

        manifest = load_input_data(self.ld_manifest_path, header_lines=0)

        bim_col = _first_existing(manifest, ["gwas_bim", "bim", "ld_bim", "bim_path"])

        # Standardize locus coordinates from manifest.
        manifest[self.left_bound_col] = pd.Series(pd.to_numeric(manifest[self.left_bound_col], errors="coerce")).astype("Int64")
        manifest[self.right_bound_col] = pd.Series(pd.to_numeric(manifest[self.right_bound_col], errors="coerce")).astype("Int64")

        bim_by_locus = {}
        missing_bim_loci = []

        for _, row in manifest.iterrows():
            locus_id = row[self.locus_id_col]
            bim = self.read_bim(row[bim_col], missing_ok=True)
            bim_by_locus[locus_id] = bim
            if bim.attrs.get("missing_bim", False) or len(bim) == 0:
                missing_bim_loci.append(locus_id)

        print(f"Loaded {len(bim_by_locus)} locus BIM entries")
        print(f"Missing/empty BIM loci: {len(missing_bim_loci)}")
        if missing_bim_loci:
            print("First missing/empty BIM loci:", missing_bim_loci[:10])


