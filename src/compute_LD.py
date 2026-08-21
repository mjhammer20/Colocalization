# Imports
import shlex
import subprocess
import shutil
import argparse
from pathlib import Path
import pandas as pd
import numpy as np #type: ignore (silences pylance warning)
from helpers import load_input_data, normalize_chromosome, numeric_series

# ------------------------- Helper Function Definitions -------------------------

def _safe_tag(locus_id: str) -> str:
    """
    Create a filesystem-safe tag for a locus ID, e.g. 1:7937247-8937247 becomes 1_7937247_8937247

    Args:
        locus_id (str): The locus ID to be converted.

    Returns:
        str: A filesystem-safe version of the locus ID.

    """
    # Filesystem-safe locus tag, e.g. chr1:7937247-8937247 becomes chr1_7937247_8937247
    return str(locus_id).replace(":", "_").replace("-", "_")


def _run(cmd: list) -> int:
    """
    Run a shell command and return the exit code.

    Args:
        cmd (list): The command to run as a list of strings.

    Returns:
        int: The exit code of the command.

    """
    # Convert all elements to strings and join them into a single command string
    cmd = [str(x) for x in cmd]
    cmd_str = " ".join(shlex.quote(x) for x in cmd)

    # Prepend the command with module load for plink/1.9
    full_cmd = f"ml plink/1.9; {cmd_str}"

    # Print the command being run for debugging purposes
    print("  $", full_cmd)

    # Run the command in a bash shell and capture output
    res = subprocess.run(["bash", "-lc", full_cmd], capture_output=True, text=True)
    if res.returncode != 0:
        print("    [stderr]", res.stderr.strip()[-1000:])
        print("    [stdout]", res.stdout.strip()[-1000:])
    
    return res.returncode


#------------------------- Main Function Definitions -------------------------

def compute_ld_for_locus(bfile: str, ld_output_dir: str, tag: str, row: pd.Series, standardized_chr_key: str, standardized_left_bound_key: str, standardized_right_bound_key: str, maf_min: float, geno_max: float) -> int:
    """
    Extracts one locus from `bfile` and writes a signed-r square LD matrix, returning the number of variants in the locus. The .gwas.* naming is kept for downstream compatibility

    Args:
        bfile (str): Path to the reference PLINK binary fileset (prefix).
        ld_output_dir (str): Directory to write output LD matrices.
        tag (str): Filesystem-safe tag for the locus.
        row (pd.Series): A row from the loci DataFrame containing locus information.
        standardized_chr_key (str): Column name for chromosome in the loci DataFrame.
        standardized_left_bound_key (str): Column name for left boundary in the loci DataFrame.
        standardized_right_bound_key (str): Column name for right boundary in the loci DataFrame.
        maf_min (float): Minimum minor allele frequency for filtering.
        geno_max (float): Maximum missing genotype rate for filtering.

    Produces:
      {tag}.gwas.ld: signed r LD matrix from PLINK --r square
      {tag}.gwas.bim: variant order + alleles

    Returns:
        int: Number of variants in the locus, or -1 if an error occurred.
    """
    # Create paths for the region and LD output files
    region = str(Path(ld_output_dir) / f"{tag}.gwas.region")
    ldpref = str(Path(ld_output_dir) / f"{tag}.gwas")

    # Extract locus
    rc = _run([
        "plink",
        "--bfile", bfile,
        "--chr", row[standardized_chr_key],
        "--from-bp", row[standardized_left_bound_key],
        "--to-bp", row[standardized_right_bound_key],
        "--maf", str(maf_min),
        "--geno", str(geno_max),
        "--keep-allele-order",
        "--make-bed",
        "--out", region
    ])
    if rc != 0:
        return -1

    # Compute signed LD matrix
    ld_cmd = [
        "plink",
        "--bfile", region,
        "--keep-allele-order",
        "--r", "square",
        "--out", ldpref
    ]
    bim = Path(region + ".bim")
    if not bim.exists():
        print(f"    region empty for {tag}")
        return 0
    nvar = sum(1 for _ in open(bim))
    rc = _run(ld_cmd)
    if rc != 0:
        return -1

    # Copy BIM next to LD so R knows SNP order + alleles
    shutil.copyfile(region + ".bim", ldpref + ".bim")

    # Remove temporary region files
    for ext in (".bed", ".bim", ".fam", ".log", ".nosex"):
        p = Path(region + ext)
        if p.exists():
            p.unlink()

    return nvar


# ------------------------- Full Pipeline -------------------------
def main(args: argparse.Namespace):
    """
    Main function to compute LD matrices for loci.

    """
    # Load loci from the input file
    loci = load_input_data(Path(args.loci_file), header_lines=args.header_lines)
    
    # Create the LD output directory if it doesn't exist
    ld_output_dir = Path(args.ld_output_dir)
    ld_output_dir.mkdir(parents=True, exist_ok=True)

    # Initialize an empty list to hold manifest entries
    manifest = []

    # Loop through each locus and compute LD
    for _, row in loci.iterrows():
        # Create a filesystem-safe tag for the locus ID
        tag = _safe_tag(str(row[args.standardized_locus_id_key]))
        print(f"\n[{row[args.standardized_locus_id_key]}]")

        # Compute LD for the locus and get the number of variants
        n = compute_ld_for_locus(args.ref_bfile, args.ld_output_dir, tag, row, args.standardized_chr_key, args.standardized_left_bound_key, args.standardized_right_bound_key, args.maf_min, args.geno_max)

        # Append the locus information and LD file paths to the manifest
        manifest.append({
            args.standardized_locus_id_key: row[args.standardized_locus_id_key],
            "tag": tag,
            args.standardized_chr_key: row[args.standardized_chr_key],
            args.standardized_left_bound_key: row[args.standardized_left_bound_key],
            args.standardized_right_bound_key: row[args.standardized_right_bound_key],
            "ld_panel": args.ld_panel,
            "ld_bfile": args.ref_bfile,
            args.ld_file_key: f"{tag}.gwas.ld",
            args.bim_file_key: f"{tag}.gwas.bim",
            "n": n,
            "same_panel": True,
            "note": "Same EUR LD reference used for GWAS and QTL."
        })

    # 

    # Convert the manifest list to a DataFrame and write it to a TSV file
    man = pd.DataFrame(manifest)

    # Normalize chromosome representation in the manifest
    man[args.standardized_chr_key] = man[args.standardized_chr_key].astype(str).map(normalize_chromosome)

    # Ensure left and right boundaries are numeric
    man[args.standardized_left_bound_key] = numeric_series(man, args.standardized_left_bound_key)
    man[args.standardized_right_bound_key] = numeric_series(man, args.standardized_right_bound_key)

    man_path = Path(args.ld_output_dir) / args.ld_manifest
    man.to_csv(man_path, sep="\t", index=False)
    print(f"\nWrote LD manifest: {man_path}")


# ------------------------- Command-Line Interface -------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute LD matrices for loci.")
    parser.add_argument("--loci_file", required=True, help="Path to the input loci file (TSV).")
    parser.add_argument("--header_lines", type=int, default=0, help="Number of header lines to skip in the loci file.")
    parser.add_argument("--standardized_locus_id_key", default="LOCUS_ID", help="Column name for locus IDs in the loci file.")
    parser.add_argument("--ref_bfile", required=True, help="Path to the reference PLINK binary fileset (prefix).")
    parser.add_argument("--ld_panel", required=True, help="Name of the LD reference panel (e.g., 1000G_EUR).")
    parser.add_argument("--ld_output_dir", required=True, help="Directory to write output LD matrices and manifest.")
    parser.add_argument("--ld_manifest", default="ld_manifest.tsv", help="Name of the LD manifest file to be created (default: ld_manifest.tsv).")
    parser.add_argument("--standardized_chr_key", default="CHR", help="Column name for chromosome in the loci file.")
    parser.add_argument("--standardized_left_bound_key", default="LEFT_500KB", help="Column name for left boundary in the loci file.")
    parser.add_argument("--standardized_right_bound_key", default="RIGHT_500KB", help="Column name for right boundary in the loci file.")
    parser.add_argument("--maf_min", type=float, default=0.01, help="Minimum minor allele frequency for filtering.")
    parser.add_argument("--geno_max", type=float, default=0.05, help="Maximum missing genotype rate for filtering.")
    parser.add_argument("--ld_file_key", default="LD", help="Column name for LD file paths in the manifest (default: ld_file).")
    parser.add_argument("--bim_file_key", default="BIM", help="Column name for BIM file paths in the manifest (default: bim_file).")

    args = parser.parse_args()

    main(args)