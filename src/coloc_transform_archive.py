# def _safe_get(row, col: Optional[str]):
#     """
#     Safely get a value from a DataFrame row by column name.
#     """
#     return row[col] if (col is not None and col in row.index) else None

# def extract_maf_from_record(record: dict, alt_allele: str) -> float:
#     """
#     Extract the minor allele frequency (MAF) from the NCBI dbSNP record.
#     """

#     # Extract GLOBAL_MAFS field from the record
#     global_mafs = record.get('GLOBAL_MAFS', None)
#     if not global_mafs:
#         print(f"Error: GLOBAL_MAFS not found in record for rsID rs{record.get('@uid', 'Unknown')}")
#         return np.nan

#     # Extract MAFs list from the GLOBAL_MAFS field
#     mafs_list = global_mafs.get('MAF', None)
#     if not mafs_list:
#         print(f"Error: MAFs not found in GLOBAL_MAFS for rsID rs{record.get('@uid', 'Unknown')}")
#         return np.nan

#     # Convert the MAFs list to a dictionary for easier access
#     try:
#         mafs_dict = {x['STUDY']: x['FREQ'] for x in mafs_list}
#     except TypeError:
#         # Records with single MAF entry may be a dictionary instead of a list, so handle that case
#         try:
#             mafs_dict = {mafs_list['STUDY']: mafs_list['FREQ']}

#         # Handle cases where the MAFs list is not in the expected format
#         except Exception as e:
#             print(f"Error: MAFs field not in expected format for rsID rs{record.get('@uid', 'Unknown')}: {e}")
#             return np.nan

#     # Extract MAF for specific studies, prioritizing GnomAD_genomes, then 1000Genomes, and finally any available study
#     maf = mafs_dict.get('GnomAD_genomes') or mafs_dict.get('1000Genomes') or list(mafs_dict.values())[0]

#     # Parse MAF string to extract the frequency of the alternate allele
#     maf_parsed = maf.split("=")
#     ma = maf_parsed[0]
#     freqs = maf_parsed[1].split("/")

#     # Determine which allele is the minor allele and return its frequency
#     if ma == alt_allele:
#         return float(freqs[0])
#     else:
#         return float(f'0.{freqs[1]}')

# def check_allele_match(
#         ref: str,
#         alt: str,
#         non_effect_allele: str,
#         effect_allele: str
#     ) -> bool:
#     """
#     Check if the alleles from the summary statistics match the alleles from dbSNP.

#     Args:
#         ref (str): The reference allele from dbSNP.
#         alt (str): The alternate allele from dbSNP.
#         effect_allele (str): The effect allele from the summary statistics.
#         non_effect_allele (str): The non-effect allele from the summary statistics.

#     Returns:
#         bool: True if the alleles match, False otherwise.
#     """

#     if not ref == non_effect_allele and alt == effect_allele:
#         return False
#     return True

# def add_variant_id(
#         chr: str,
#         pos: int,
#         ref: str,
#         alt: str
#     ) -> str:
#     """
#     Create a marker ID in the format chrN:bp:ref:alt.

#     Args:
#         chr (str): The chromosome number or identifier.
#         pos (int): The base pair position on the chromosome.
#         ref (str): The reference allele.
#         alt (str): The alternate allele.

#     Returns:
#         str: The marker ID in the format chr:bp:ref:alt.
#     """
#     return f"{chr}:{pos}:{ref}:{alt}"

# def calculate_n(
#         mac: Optional[int] = None,
#         maf: Optional[float] = None
#     ) -> float:
#     """
#     Calculate the sample size (N) from minor allele count (MAC) and minor allele frequency (MAF).

#     Args:
#         mac (Optional[int]): The minor allele count.
#         maf (Optional[float]): The minor allele frequency.

#     Returns:
#         float: The calculated sample size or None if insufficient information is provided.
#     """

#     if mac is not None and maf is not None:
#         return float(mac / (2 * float(maf)))

#     return np.nan

# def calculate_se(
#         beta: Optional[float] = None,
#         t: Optional[float] = None,
#         z: Optional[float] = None,
#         p: Optional[float] = None,
#         n: Optional[int] = None,
#         maf: Optional[float] = None
#     ) -> float:
#     """
#     Calculate the standard error (SE) from available summary statistics

#     Args:
#         beta (Optional[float]): The beta coefficient.
#         t (Optional[float]): The t-statistic.
#         z (Optional[float]): The z-score.
#         p (Optional[float]): The p-value.
#         n (Optional[int]): The sample size.
#         maf (Optional[float]): The minor allele frequency.
#         mac (Optional[int]): The minor allele count.

#     Returns:
#         float: The calculated standard error or None if insufficient information is provided.
#     """

#     # Calculate SE from beta, if available, using t, z, or p-value
#     if beta is not None:

#         # Calculate SE from beta using t-statistic, if available
#         if t is not None:
#             return abs(float(beta) / float(t))

#         # If t-statistic is not available, calculate SE from beta using z-score, if available
#         elif z is not None:
#             return abs(float(beta) / float(z))

#         # If t-statistic and z-score are not available, calculate SE from beta using p-value, if available
#         elif p is not None:
#             z = float(norm.ppf(1 - float(p) / 2))
#             return abs(float(beta) / z)

#     # If beta is not available, calculate SE from MAF, if available, using sample size or minor allele count
#     elif maf is not None:

#         # Calculate SE from MAF using sample size, if available
#         if n is not None and pd.isna(n) == False:
#             return 1 / np.sqrt(2 * float(maf) * (1 - float(maf)) * int(n))

#     return np.nan

# def calculate_pvalue(
#         beta: Optional[float] = None,
#         se: Optional[float] = None,
#         t: Optional[float] = None,
#         z: Optional[float] = None
#     ) -> float:
#     """
#     Calculate the p-value from available summary statistics.

#     Args:
#         beta (Optional[float]): The beta coefficient.
#         se (Optional[float]): The standard error.
#         t (Optional[float]): The t-statistic.
#         z (Optional[float]): The z-score.

#     Returns:
#         float: The calculated p-value or None if insufficient information is provided.
#     """

#     # Calculate p-value from beta and SE, if available
#     if beta is not None and se is not None:
#         z = float(beta) / float(se)
#         return float(2 * norm.sf(abs(z)))

#     # If beta and SE are not available, calculate p-value from t-statistic, if available
#     elif t is not None:
#         return float(2 * norm.sf(abs(float(t))))

#     # If beta, SE, and t-statistic are not available, calculate p-value from z-score, if available
#     elif z is not None:
#         return float(2 * norm.sf(abs(float(z))))

#     return np.nan

# def calculate_beta(
#         se: Optional[float] = None,
#         t: Optional[float] = None,
#         z: Optional[float] = None,
#         p: Optional[float] = None,
#         n: Optional[int] = None,
#         maf: Optional[float] = None
#     ) -> float:
#     """
#     Calculate the beta coefficient from available summary statistics.

#     Args:
#         se (Optional[float]): The standard error.
#         t (Optional[float]): The t-statistic.
#         z (Optional[float]): The z-score.
#         p (Optional[float]): The p-value.
#         n (Optional[int]): The sample size.
#         maf (Optional[float]): The minor allele frequency.

#     Returns:
#         float: The calculated beta coefficient or None if insufficient information is provided.
#     """

#     # Calculate beta from SE, if available, using t-statistic, z-score, or p-value
#     if se is not None:

#         # Calculate beta from SE using t-statistic, if available
#         if t is not None:
#             return float(t) * float(se)

#         # If t-statistic is not available, calculate beta from SE using z-score, if available
#         elif z is not None:
#             print(f"Calculating beta from SE using z-score: z={z}, se={se}")
#             beta = float(z) * float(se)
#             print(f"Calculated beta: {beta}")
#             return beta

#         # If t-statistic and z-score are not available, calculate beta from SE using p-value, if available
#         elif p is not None:
#             print(f"Calculating beta from SE using p-value: p={p}, se={se}")
#             z_from_p = float(norm.ppf(1 - float(p) / 2))
#             return float(float(se) * z_from_p)

#     # If SE is not available, calculate beta from MAF, if available, using sample size or minor allele count, and t-statistic, z-score, or p-value
#     elif maf is not None and n is not None and pd.isna(n) == False:
            
#         # Calculate beta from MAF using t-statistic and sample size, if available
#         if t is not None:
#             return float(1 / np.sqrt(2 * float(maf) * (1 - float(maf)) * (int(n) + float(t)**2)))
        
#         # If t-statistic is not available, calculate beta from MAF using z-score and sample size, if available
#         elif z is not None:
#             return float(1 / np.sqrt(2 * float(maf) * (1 - float(maf)) * (int(n) + float(z)**2)))

#         # If t-statistic and z-score are not available, calculate beta from MAF using p-value and sample size, if available
#         elif p is not None:
#             z_from_p = float(norm.ppf(1 - float(p) / 2))
#             return float(1 / np.sqrt(2 * float(maf) * (1 - float(maf)) * (int(n) + float(z_from_p)**2)))

#     return np.nan

# Function to transform coloc results
# def annotate_sum_stats(
#         df: pd.DataFrame,
#         chr_col: str,
#         pos_col: str,  
#         rsid_col: Optional[str] = None,
#         non_effect_allele_col: Optional[str] = None,
#         effect_allele_col: Optional[str] = None,
#         pval_col: Optional[str] = None,
#         beta_col: Optional[str] = None,
#         se_col: Optional[str] = None,
#         n_col: Optional[str] = None,
#         t_col: Optional[str] = None,
#         z_col: Optional[str] = None, 
#         maf_col: Optional[str] = None,
#         mac_col: Optional[str] = None
#     ) -> pd.DataFrame:
#     """
#     Transform the coloc results DataFrame by selecting relevant columns and renaming them.

#     Args:
#         df (pd.DataFrame): The input DataFrame containing coloc results.
#         chr_col (str): The name of the column containing chromosome information.
#         pos_col (str): The name of the column containing position information.
#         rsid_col (Optional[str]): The name of the column containing rsID information.
#         non_effect_allele_col (Optional[str]): The name of the column containing non effect allele information.
#         effect_allele_col (Optional[str]): The name of the column containing effect allele information.
#         pval_col (Optional[str]): The name of the column containing p-value information.
#         beta_col (Optional[str]): The name of the column containing beta coefficient information.
#         se_col (Optional[str]): The name of the column containing standard error information.
#         n_col (Optional[str]): The name of the column containing sample size information.
#         t_col (Optional[str]): The name of the column containing t-statistic information.
#         z_col (Optional[str]): The name of the column containing z-score information.
#         maf_col (Optional[str]): The name of the column containing minor allele frequency information.
#         mac_col (Optional[str]): The name of the column containing minor allele count information.

#     Returns:
#         pd.DataFrame: The transformed DataFrame with selected and renamed columns.
#     """

#     # Normalize chromosome column
#     df[chr_col] = df[chr_col].astype(str).apply(normalize_chromosome)

#     # Subset the DataFrame to unique SNPs based on chromosome and position, and fetch rsID if not provided
#     if rsid_col is None:
#         rsid_col = 'rsID'
#         snps_df = df[[chr_col, pos_col]].drop_duplicates()
#         snps_df[rsid_col] = df.apply(lambda row: fetch_rsid(row[chr_col], row[pos_col]), axis=1)
    
#     else:
#         snps_df = df[[chr_col, pos_col, rsid_col]].drop_duplicates()

#     # Add dbSNP information (ref allele, alt allele, MAF) based on rsID, skip if rsID is NULL/empty
#     ref_allele_col_dbsnp = 'ref_dbSNP'
#     alt_allele_col_dbsnp = 'alt_dbSNP'
#     maf_col_dbsnp = 'maf_dbSNP'
#     snps_df[ref_allele_col_dbsnp], snps_df[alt_allele_col_dbsnp], snps_df[maf_col_dbsnp] = zip(*snps_df[rsid_col].apply(lambda rsid: add_dbSNP_info(rsid) if not rsid in ["NA", ".", ""] else ("NA", "NA", np.nan)))

#     # Merge the dbSNP information back into the original DataFrame
#     if rsid_col not in df.columns:
#         annotated_df = pd.merge(df, snps_df, on=[chr_col, pos_col], how='left')
#     else:
#         annotated_df = pd.merge(df, snps_df, on=[chr_col, pos_col, rsid_col], how='left')

#     # Drop MAF from dbSNP if maf_col is provided
#     if maf_col is None:
#         maf_col = maf_col_dbsnp

#     # Check provided effect/non-effect alleles from sum stats against dbSNP alleles, if available. If not available, set to NA.
#     allele_match_col = 'Allele_Match'
#     annotated_df[allele_match_col] = annotated_df.apply(lambda row: check_allele_match(
#         ref= row[ref_allele_col_dbsnp],
#         alt= row[alt_allele_col_dbsnp],
#         non_effect_allele= row[non_effect_allele_col],
#         effect_allele= row[effect_allele_col])
#         if all(_safe_get(row, col) for col in [non_effect_allele_col, effect_allele_col])
#         else "NA",
#         axis=1
#     )

#     # Add variant ID column in the format chr:bp:ref:alt. Use ref/alt from dbSNP if available, otherwise use non-effect/effect alleles if available, else set to NA.
#     variant_id_col = 'Variant_ID'
#     annotated_df[variant_id_col] = annotated_df.apply(lambda row: add_variant_id(
#         chr=row[chr_col],
#         pos=row[pos_col],
#         ref=row[ref_allele_col_dbsnp],
#         alt=row[alt_allele_col_dbsnp]) 
#         if row[ref_allele_col_dbsnp] != "NA" and row[alt_allele_col_dbsnp] != "NA" 
#         else add_variant_id(
#             chr=row[chr_col],
#             pos=row[pos_col],
#             ref=row[non_effect_allele_col],
#             alt=row[effect_allele_col])
#         if all(_safe_get(row, col) for col in [non_effect_allele_col, effect_allele_col])
#         else "NA",
#         axis=1
#     )

#     # Add se column if not provided
#     if se_col is None:
#         se_col = 'SE'
#         annotated_df[se_col] = annotated_df.apply(
#             lambda row: calculate_se(
#                 beta=_safe_get(row, beta_col),
#                 t=_safe_get(row, t_col),
#                 z=_safe_get(row, z_col),
#                 p=_safe_get(row, pval_col),
#                 n=_safe_get(row, n_col),
#                 maf=_safe_get(row, maf_col),
#                 mac=_safe_get(row, mac_col)),
#             axis=1
#         )

#     # Add n column if not provided
#     if n_col is None:
#         n_col = 'N'
#         annotated_df[n_col] = annotated_df.apply(
#             lambda row: calculate_n_from_mac_and_maf(
#                 mac=_safe_get(row, mac_col),
#                 maf=_safe_get(row, maf_col)),
#             axis=1
#         )

#     # Add beta column if not provided
#     if beta_col is None:
#         beta_col = 'Beta'
#         annotated_df[beta_col] = annotated_df.apply(
#             lambda row: calculate_beta(
#                 se=_safe_get(row, se_col),
#                 t=_safe_get(row, t_col),
#                 z=_safe_get(row, z_col),
#                 p=_safe_get(row, pval_col),
#                 n=_safe_get(row, n_col),
#                 maf=_safe_get(row, maf_col), 
#                 mac=_safe_get(row, mac_col)),
#             axis=1
#         )

#     # Add p-value column if not provided
#     if pval_col is None:
#         pval_col = 'P'
#         annotated_df[pval_col] = annotated_df.apply(
#             lambda row: calculate_pvalue(
#                 beta=_safe_get(row, beta_col),
#                 se=_safe_get(row, se_col),
#                 t=_safe_get(row, t_col),
#                 z=_safe_get(row, z_col)),
#             axis=1
#         )

#     return annotated_df

# def annotate_sum_stats(
#         df: pd.DataFrame,
#         chr_col: str,
#         pos_col: str,  
#         rsid_col: Optional[str] = None,
#         non_effect_allele_col: Optional[str] = None,
#         effect_allele_col: Optional[str] = None,
#         pval_col: Optional[str] = None,
#         beta_col: Optional[str] = None,
#         se_col: Optional[str] = None,
#         n_col: Optional[str] = None,
#         t_col: Optional[str] = None,
#         z_col: Optional[str] = None, 
#         maf_col: Optional[str] = None,
#         mac_col: Optional[str] = None
#     ) -> pd.DataFrame:
#     """
#     Transform the coloc results DataFrame by selecting relevant columns and renaming them.

#     Args:
#         df (pd.DataFrame): The input DataFrame containing coloc results.
#         chr_col (str): The name of the column containing chromosome information.
#         pos_col (str): The name of the column containing position information.
#         rsid_col (Optional[str]): The name of the column containing rsID information.
#         non_effect_allele_col (Optional[str]): The name of the column containing non effect allele information.
#         effect_allele_col (Optional[str]): The name of the column containing effect allele information.
#         pval_col (Optional[str]): The name of the column containing p-value information.
#         beta_col (Optional[str]): The name of the column containing beta coefficient information.
#         se_col (Optional[str]): The name of the column containing standard error information.
#         n_col (Optional[str]): The name of the column containing sample size information.
#         t_col (Optional[str]): The name of the column containing t-statistic information.
#         z_col (Optional[str]): The name of the column containing z-score information.
#         maf_col (Optional[str]): The name of the column containing minor allele frequency information.
#         mac_col (Optional[str]): The name of the column containing minor allele count information.

#     Returns:
#         pd.DataFrame: The transformed DataFrame with selected and renamed columns.
#     """

#     # Normalize chromosome column
#     df[chr_col] = df[chr_col].astype(str).apply(normalize_chromosome)

#     # Subset the DataFrame to unique SNPs based on chromosome and position, and fetch rsID if not provided
#     if rsid_col is None:
#         rsid_col = 'rsID'
#         snps_df = df[[chr_col, pos_col]].drop_duplicates()
#         snps_df[rsid_col] = df.apply(lambda row: fetch_rsid(row[chr_col], row[pos_col]), axis=1)
    
#     else:
#         snps_df = df[[chr_col, pos_col, rsid_col]].drop_duplicates()

#     # Add dbSNP information (ref allele, alt allele, MAF) based on rsID, skip if rsID is NULL/empty
#     if non_effect_allele_col is None or effect_allele_col is None or maf_col is None:
#         non_effect_allele_col = 'Non_Effect_Allele'
#         effect_allele_col = 'Effect_Allele'
#         maf_col = 'MAF'
#         snps_df[non_effect_allele_col], snps_df[effect_allele_col], snps_df[maf_col] = zip(*snps_df[rsid_col].apply(lambda rsid: add_dbSNP_info(rsid) if not rsid in ["NA", ".", ""] else ("NA", "NA", np.nan)))

#         # Merge the dbSNP information back into the original DataFrame
#         if rsid_col not in df.columns:
#             annotated_df = pd.merge(df, snps_df, on=[chr_col, pos_col], how='left')
#         else:
#             annotated_df = pd.merge(df, snps_df, on=[chr_col, pos_col, rsid_col], how='left')
    
#     else:
#         annotated_df = df.copy()

#     # Add variant ID column in the format chr:bp:ref:alt. Use ref/alt from dbSNP if available, otherwise use non-effect/effect alleles if available, else set to NA.
#     variant_id_col = 'Variant_ID'
#     annotated_df[variant_id_col] = annotated_df.apply(lambda row:  add_variant_id(
#         chr=row[chr_col],
#         pos=row[pos_col],
#         ref=row[non_effect_allele_col],
#         alt=row[effect_allele_col])
#         if all(_safe_get(row, col) for col in [non_effect_allele_col, effect_allele_col])
#         else "NA",
#         axis=1
#     )

#     # Add n column if not provided
#     if n_col is None:
#         n_col = 'N'
#         annotated_df[n_col] = annotated_df.apply(
#             lambda row: calculate_n(
#                 mac=_safe_get(row, mac_col),
#                 maf=_safe_get(row, maf_col)),
#             axis=1
#         )

#     # Add se column if not provided
#     if se_col is None:
#         se_col = 'SE'
#         annotated_df[se_col] = annotated_df.apply(
#             lambda row: calculate_se(
#                 beta=_safe_get(row, beta_col),
#                 t=_safe_get(row, t_col),
#                 z=_safe_get(row, z_col),
#                 p=_safe_get(row, pval_col),
#                 n=_safe_get(row, n_col),
#                 maf=_safe_get(row, maf_col)),
#             axis=1
#         )

#     # Add beta column if not provided
#     if beta_col is None:
#         beta_col = 'Beta'
#         annotated_df[beta_col] = annotated_df.apply(
#             lambda row: calculate_beta(
#                 se=_safe_get(row, se_col),
#                 t=_safe_get(row, t_col),
#                 z=_safe_get(row, z_col),
#                 p=_safe_get(row, pval_col),
#                 n=_safe_get(row, n_col),
#                 maf=_safe_get(row, maf_col)),
#             axis=1
#         )

#     # Add p-value column if not provided
#     if pval_col is None:
#         pval_col = 'P'
#         annotated_df[pval_col] = annotated_df.apply(
#             lambda row: calculate_pvalue(
#                 beta=_safe_get(row, beta_col),
#                 se=_safe_get(row, se_col),
#                 t=_safe_get(row, t_col),
#                 z=_safe_get(row, z_col)),
#             axis=1
#         )

#     return annotated_df

# def main(args: argparse.Namespace):
#     """
#     Main function to load input data, annotate summary statistics, and save the transformed data.

#     Args:
#         args (argparse.Namespace): Command line arguments.
#     """

#     # Define input file path
#     input_fp = Path(args.input_file)

#     # Load input data
#     df = load_input_data(input_fp, args.header_lines)

#     # Annotate summary statistics
#     annotated_df = annotate_sum_stats(
#         df=df,
#         chr_col=args.chr_col,
#         pos_col=args.pos_col,
#         rsid_col=args.rsid_col,
#         non_effect_allele_col=args.non_effect_allele_col,
#         effect_allele_col=args.effect_allele_col,
#         pval_col=args.pval_col,
#         beta_col=args.beta_col,
#         se_col=args.se_col,
#         n_col=args.n_col,
#         t_col=args.t_col,
#         z_col=args.z_col,
#         maf_col=args.maf_col,
#         mac_col=args.mac_col
#     )

#     # Filter out rows with mismatched alleles in the annotated DataFrame
#     matched_df = annotated_df[annotated_df['Allele_Match'] == True]

#     # Transform summary statistics
#     transformed_df = transform_sum_stats(
#         df=matched_df,
#         chr_col=args.chr_col,
#         pos_col=args.pos_col,
#         variant_id_col='Variant_ID',
#         rsid_col=args.rsid_col or 'rsID',
#         gene_id_col=args.gene_id_col,
#         non_effect_allele_col=args.non_effect_allele_col,
#         effect_allele_col=args.effect_allele_col,
#         pval_col=args.pval_col or 'P',
#         beta_col=args.beta_col or 'Beta',
#         se_col=args.se_col or 'SE',
#         maf_col=args.maf_col or 'maf_dbSNP',
#         n_col=args.n_col or 'N'
#     )

#     # Filter out rows with NA values in the transformed DataFrame
#     standardized_df = transformed_df.dropna()

#     # Create a DataFrame for non-standardized rows (mismatched alleles or NA values)
#     non_standardized_df = annotated_df[~annotated_df.index.isin(standardized_df.index)]

#     # Extract the base name of the input file for naming the output files
#     if input_fp.suffix == '.gz':
#         input_file_stem = input_fp.stem  # Remove .gz
#         input_file_base = Path(input_file_stem).stem  # Remove file extension (e.g., .csv, .tsv, .txt)
#     else:
#         input_file_base = input_fp.stem  # Remove file extension (e.g., .csv, .tsv, .txt)
    
#     # Save the transformed DataFrame to the output file
#     standardized_fp = Path(f"{args.output_dir}/{input_file_base}.coloc_ready.tsv")
#     standardized_df.to_csv(standardized_fp, sep='\t', index=False)

#     # Save the mismatched alleles DataFrame to a separate file
#     non_standardized_fp = Path(f"{args.output_dir}/{input_file_base}.not_coloc_ready.tsv")
#     non_standardized_df.to_csv(non_standardized_fp, sep='\t', index=False)