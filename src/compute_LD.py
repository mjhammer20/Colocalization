# Imports
import shlex
import subprocess
import shutil
import argparse
from pathlib import Path
import pandas as pd
import numpy as np #type: ignore (silences pylance warning)
from helpers import load_input_data

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
    cmd = [str(x) for x in cmd]
    cmd_str = " ".join(shlex.quote(x) for x in cmd)

    # Keep this hard-coded
    full_cmd = f"ml plink/1.9; {cmd_str}"

    print("  $", full_cmd)

    res = subprocess.run(["bash", "-lc", full_cmd], capture_output=True, text=True)
    if res.returncode != 0:
        print("    [stderr]", res.stderr.strip()[-1000:])
        print("    [stdout]", res.stdout.strip()[-1000:])
    return res.returncode
        
# ------------------------- Class Definition -------------------------

def compute_ld_for_locus(bfile: str, ld_output_dir: str, tag: str, row: pd.Series, chr_col: str, left_bound_col: str, right_bound_col: str, maf_min: float, geno_max: float) -> int:
    """
    Extracts one locus from `bfile` and writes a signed-r square LD matrix, returning the number of variants in the locus. The .gwas.* naming is kept for downstream compatibility

    Args:
        bfile (str): Path to the reference PLINK binary fileset (prefix).
        ld_output_dir (str): Directory to write output LD matrices.
        tag (str): Filesystem-safe tag for the locus.
        row (pd.Series): A row from the loci DataFrame containing locus information.
        chr_col (str): Column name for chromosome in the loci DataFrame.
        left_bound_col (str): Column name for left boundary in the loci DataFrame.
        right_bound_col (str): Column name for right boundary in the loci DataFrame.
        maf_min (float): Minimum minor allele frequency for filtering.
        geno_max (float): Maximum missing genotype rate for filtering.

    Produces:
      {tag}.gwas.ld: signed r LD matrix from PLINK --r square
      {tag}.gwas.bim: variant order + alleles

    Returns:
        int: Number of variants in the locus, or -1 if an error occurred.
    """
    region = str(Path(ld_output_dir) / f"{tag}.gwas.region")
    ldpref = str(Path(ld_output_dir) / f"{tag}.gwas")

    # 1) Extract locus
    rc = _run([
        "plink",
        "--bfile", bfile,
        "--chr", row[chr_col],
        "--from-bp", row[left_bound_col],
        "--to-bp", row[right_bound_col],
        "--maf", str(maf_min),
        "--geno", str(geno_max),
        "--keep-allele-order",
        "--make-bed",
        "--out", region
    ])
    if rc != 0:
        return -1

    # 2) Compute signed LD matrix
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

# ------------------------- Main Function Definition -------------------------
def main(args: argparse.Namespace):

    loci = load_input_data(args.loci_file, header_lines=args.header_lines)
    manifest = []

    for _, row in loci.iterrows():
        tag = _safe_tag(str(row[args.locus_id_col]))
        print(f"\n[{row[args.locus_id_col]}]")

        n = compute_ld_for_locus(args.ref_bfile, args.ld_output_dir, tag, row, args.chr_col, args.left_bound_col, args.right_bound_col, args.maf_min, args.geno_max)

        manifest.append({
            args.locus_id_col: row[args.locus_id_col],
            "tag": tag,
            args.chr_col: row[args.chr_col],
            args.left_bound_col: row[args.left_bound_col],
            args.right_bound_col: row[args.right_bound_col],
            "ld_panel": args.ld_panel,
            "ld_bfile": args.ref_bfile,

            args.ld_file_col: f"{tag}.gwas.ld",
            args.bim_file_col: f"{tag}.gwas.bim",

            "n_gwas": n,
            "n_eqtl": n,
            "same_panel": True,
            "note": "Same EUR LD reference used for GWAS and GTEx/eQTL because GTEx genotypes are unavailable."
        })

    man = pd.DataFrame(manifest)
    man_path = Path(args.ld_output_dir) / args.ld_manifest
    man.to_csv(man_path, sep="\t", index=False)

    print(f"\nWrote LD manifest: {man_path}")


# ------------------------- Command-Line Interface -------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute LD matrices for loci.")
    parser.add_argument("--loci_file", required=True, help="Path to the input loci file (TSV).")
    parser.add_argument("--header_lines", type=int, default=0, help="Number of header lines to skip in the loci file.")
    parser.add_argument("--locus_id_col", required=True, help="Column name for locus IDs in the loci file.")
    parser.add_argument("--ref_bfile", required=True, help="Path to the reference PLINK binary fileset (prefix).")
    parser.add_argument("--ld_panel", required=True, help="Name of the LD reference panel (e.g., 1000G_EUR).")
    parser.add_argument("--ld_output_dir", required=True, help="Directory to write output LD matrices and manifest.")
    parser.add_argument("--ld_manifest", default="ld_manifest.tsv", help="Name of the LD manifest file to be created (default: ld_manifest.tsv).")
    parser.add_argument("--chr_col", required=True, help="Column name for chromosome in the loci file.")
    parser.add_argument("--left_bound_col", required=True, help="Column name for left boundary in the loci file.")
    parser.add_argument("--right_bound_col", required=True, help="Column name for right boundary in the loci file.")
    parser.add_argument("--maf_min", type=float, default=0.01, help="Minimum minor allele frequency for filtering.")
    parser.add_argument("--geno_max", type=float, default=0.05, help="Maximum missing genotype rate for filtering.")

    args = parser.parse_args()

    main(args)