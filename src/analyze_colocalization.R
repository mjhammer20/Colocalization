# Imports

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(purrr)
  library(arrow)
  library(susieR)
  library(coloc)
  library(R6)
})

# ------------------------------ Helper Function Definitions -------------------------------

.log <- function(...) {
    """
    Logs a message with a timestamp to both the console and a log file.

    Args:
        ...: The message to log, which can be formatted using sprintf-style formatting.
    
    Returns:
        None. The function prints the message to the console and appends it to a log file.

    """
    # Format the message with a timestamp and write it to both the console and the log file
    msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), " | ", sprintf(...))
    cat(msg, "\n"); write(msg, file = log_file, append = TRUE)
}


.norm_chr <- function(x) {
    """
    Standardizes chromosome labels by converting them to uppercase and ensuring they start with 'chr'.

    Args:
        x: A character vector of chromosome labels (e.g., 'chr1', 'CHRX', '2').
    
    Returns:
        A character vector of standardized chromosome labels (e.g., 'chr1', 'chrX', 'chr2').

    """
    ifelse(grepl("^chr", x, ignore.case = TRUE), substr(x, 4, nchar(x)), x)

}


.add_variant_id <- function(chr, pos, a1, a2) {
    """
    Constructs a variant ID for a genetic variant based on chromosome, position, and alleles.

    Args:
        chr: Chromosome identifier (e.g., 'chr1', 'chrX').
        pos: Position of the variant on the chromosome (integer).
        a1: First allele (string).
        a2: Second allele (string).

    Returns:
        A string representing the canonical key for the variant in the format 'chr:pos:allele1_allele2'

    """
    # Construct the canonical key using the standardized chromosome, position, and ordered alleles
    paste0(.norm_chr(chr), ":", as.integer(pos), ":", a1, ":", a2)
}

.is_ambiguous <- function(a1, a2) {
    """
    Checks if a pair of alleles is ambiguous (i.e., A/T or C/G).

    Args:
        a1: First allele (string).
        a2: Second allele (string).

    Returns:
        TRUE if the allele pair is ambiguous, FALSE otherwise.

    """
    # Construct a string representing the allele pair in uppercase
    p <- paste0(toupper(a1), toupper(a2))

    # Check if the allele pair is one of the ambiguous pairs (A/T or C/G)
    p %in% c("AT", "TA", "CG", "GC")
}


.stream_write <- function(df, path) {
    """
    Writes a data frame to a file in a streaming manner. If the file already exists, it appends the data frame to the existing file without writing column names.

    Args:
        df: The data frame to be written.
        path: The file path to write the data frame to.

    Returns:
        None. The function writes the data frame to the specified file path.
        
    """
    # Use data.table's fwrite function to write the data frame to the specified path
    if (!file.exists(path)) data.table::fwrite(df, path, sep = "\t", quote = FALSE)
    else data.table::fwrite(df, path, sep = "\t", quote = FALSE, append = TRUE, col.names = FALSE)
}


.credset_for_row <- function(coloc_result, row_k) {
    """
    Extracts the credible set size and top SNP for a specific row in the coloc result.

    Args:
        coloc_result: The result object from the coloc analysis.
        row_k: The index of the row to extract credible set information from.

    Returns:
        A list containing:
            - size: The size of the credible set (number of SNPs).
            - top: The SNP with the highest posterior probability in the credible set.

    """
    # Initialize output with NA values
    out <- list(size = NA_integer_, top = NA_character_)

    # Attempt to convert the coloc result to a data frame and handle any errors
    rr <- tryCatch(as.data.frame(coloc_result$results), error = function(e) NULL)

    # Check if the result is valid and contains the necessary columns
    if (is.null(rr) || !nrow(rr) || !"snp" %in% names(rr)) return(out)

    # Find columns that match the pattern for SNP posterior probabilities for hypothesis H4
    h4cols <- grep("^SNP\\.PP\\.H4", names(rr), value = TRUE)

    # If there are no matching columns, return the output with NA values
    if (!length(h4cols)) return(out)

    # Determine the appropriate column for the specified row index
    col <- if (length(h4cols) >= row_k) {
        byname <- paste0("SNP.PP.H4.row", row_k)
        if (byname %in% h4cols) byname else h4cols[row_k]
    } else h4cols[1]

    # Extract the posterior probabilities for the specified column and handle any warnings
    v <- suppressWarnings(as.numeric(rr[[col]]))
    if (!any(is.finite(v))) return(out)
    
    # Order the posterior probabilities in decreasing order and calculate the cumulative sum to determine the credible set
    o <- order(v, decreasing = TRUE)
    credible_set <- cumsum(v[o]); w <- which(credible_set >= CRED_COVERAGE)[1]
    
    # If no credible set is found, set the size to the total number of SNPs
    if (is.na(w)) w <- length(o)
    out$size <- as.integer(w)
    out$top  <- as.character(rr$snp[o][1])
    out
}

# ------------------------------ Class Definition -------------------------------

ColocalizationAnalyzer <- R6Class("ColocalizationAnalyzer",
    public = list(

        # Initialize the class with paths and column keys
        initialize = function(args) {

            # Inputs/Outputs
            self$output_dir <- args$output_dir
            self$ld_dir <- paste0(self$output_dir, args$ld_dir)
            self$gwas <- args$gwas
            self$gwas_path <- paste0(self$output_dir, args$gwas)
            self$gwas_target_allele_key <- args$gwas_target_allele_key
            self$qtl <- args$qtl
            self$qtl_path <- paste0(self$output_dir, args$qtl)
            self$qtl_target_allele_key <- args$qtl_target_allele_key
            self$ld_manifest <- args$ld_manifest
            self$ld_manifest_path <- paste0(self$ld_dir, args$ld_manifest)

            # Parameters
            self$drop_ambiguous <- args$drop_ambiguous
            self$window_bp <- args$window_bp
            self$min_overlap <- args$min_overlap
            self$susie_min_snps <- args$susie_min_snps
            self$pp_h4_threshold <- args$pp_h4_threshold
            self$coloc_cred_coverage <- args$coloc_cred_coverage
            self$susie_cred_coverage <- args$susie_cred_coverage
            self$susie_max_iter <- args$susie_max_iter
            self$susie_l <- args$susie_l
            self$susie_repeat_until_converged <- args$susie_repeat_until_converged
            self$coloc_priors <- list(
                p1 = args$coloc_prior_p1,
                p2 = args$coloc_prior_p2,
                p12 = args$coloc_prior_p12
            )

            # Standardized Column Keys
            self.standardized_variant_id_key <- args$standardized_variant_id_key
            self.standardized_chr_key <- args$standardized_chr_key
            self.standardized_pos_key <- args$standardized_pos_key
            self.standardized_effect_allele_key <- args$standardized_effect_allele_key
            self.standardized_non_effect_allele_key <- args$standardized_non_effect_allele_key
            self.standardized_beta_key <- args$standardized_beta_key
            self.standardized_var_beta_key <- args$standardized_var_beta_key
            self.standardized_se_key <- args$standardized_se_key
            self.standardized_sdy_key <- args$standardized_sdy_key
        }

        read_ld <- function(ld_path, bim_path) {
            """
            Reads LD matrix and BIM file, aligns them, and returns a list containing the LD matrix and allele information.

            Args:
                ld_path: Path to the LD matrix file.
                bim_path: Path to the BIM file.

            Returns:
                A list containing:
                    - M: The LD matrix with canonical key dimnames.
                    - A1: A named vector of the allele the LD is counted on.
                    - A2: A named vector of the other allele.
                Returns NULL if the LD or BIM files are missing or if there are issues with reading them

            """
            # Check if the LD and BIM files exist
            if (is.na(ld_path) || is.na(bim_path) ||
                !file.exists(ld_path) || !file.exists(bim_path)) {
                .log("  LD missing: ld=%s bim=%s", as.character(ld_path), as.character(bim_path))
                return(NULL)
            }

            # Read the BIM file and handle any errors
            bim <- tryCatch(fread(bim_path, header = FALSE, data.table = FALSE),
                            error = function(e) NULL)

            # Check if the BIM file is missing or empty
            if (is.null(bim) || !nrow(bim)) {
                .log("  BIM missing/empty: %s", as.character(bim_path))
                return(NULL)
            }

            # Extract chromosome, position, and alleles from the BIM file and create canonical keys
            # PLINK .bim: V1 chr, V2 id, V3 cM, V4 bp, V5 A1, V6 A2 ; --r counts A1
            chr <- .norm_chr(bim$V1)
            pos <- as.integer(bim$V4)
            a1  <- toupper(bim$V5)
            a2  <- toupper(bim$V6)
            keys <- .add_variant_id(chr, pos, a1, a2)

            # Read the LD matrix and handle any errors
            M <- tryCatch(as.matrix(fread(ld_path, header = FALSE, data.table = FALSE)),
                            error = function(e) NULL)
            if (is.null(M) || nrow(M) != length(keys) || ncol(M) != length(keys)) {
                .log("  LD shape mismatch (%s): matrix %dx%d vs bim %d",
                    basename(ld_path), nrow(M %||% matrix(0)), ncol(M %||% matrix(0)), length(keys))
                return(NULL)
            }

            # Keep the first occurrence of duplicate variants
            dup <- duplicated(keys)
            if (any(dup)) {
                M <- M[!dup, !dup, drop = FALSE]
                a1 <- a1[!dup]
                a2 <- a2[!dup]
                keys <- keys[!dup]
            }

            # Set dimnames and storage mode for the LD matrix, and replace non-finite values with 0
            dimnames(M) <- list(keys, keys)
            storage.mode(M) <- "numeric"
            M[!is.finite(M)] <- 0

            # Return the LD matrix and allele information as a list
            list(M = M, A1 = setNames(a1, keys), A2 = setNames(a2, keys))
        }


        align_ld <- function(ld, snps, target_allele) {
            """
            Aligns the LD matrix to a target allele for each variant, subsets it to the specified SNPs, and returns the aligned LD matrix.

            Args:
                ld: A list containing the LD matrix and allele information (output from read_ld).
                snps: A character vector of SNPs to subset the LD matrix to.
                target_allele: A named vector of target alleles for each SNP (names should match the SNPs) in which the datasets beta is on (ie. GWAS effect allele / GTEx ALT).
            
            Returns:
                A numeric matrix representing the aligned LD matrix for the specified SNPs, with rows and columns corresponding to the SNPs in the order they appear in the `snps` vector.
                Returns NULL if there are fewer than 2 SNPs present or if the LD matrix is NULL.
            """
            # Check if the LD matrix is NULL and return NULL if it is
            if (is.null(ld)) return(NULL)

            # Subset the SNPs to those present in the LD matrix
            present <- snps[snps %in% rownames(ld$M)]

            # If there are fewer than 2 SNPs present, return NULL
            if (length(present) < 2) return(NULL)

            # Extract the LD matrix and allele information for the present SNPs
            M  <- ld$M[present, present, drop = FALSE]
            a1 <- ld$A1[present]
            a2 <- ld$A2[present]
            tgt <- target_allele[present]

            # Align the LD sign to the effect allele
            s <- rep(NA_real_, length(present))
            s[a1 == tgt] <- 1
            s[a2 == tgt] <- -1

            # Check for finite values in the sign vector
            ok <- is.finite(s)

            # If there are fewer than 2 SNPs with finite signs, return NULL
            if (sum(ok) < 2) return(NULL)

            # Subset the present SNPs, LD matrix, and sign vector to those with finite signs
            present <- present[ok]; M <- M[ok, ok, drop = FALSE]; s <- s[ok]

            # Align the LD matrix by multiplying it with the outer product of the sign vector
            M <- (s %o% s) * M 

            # Set the diagonal of the LD matrix to 1 and set the dimnames to the present SNPs
            diag(M) <- 1

            # Set the dimnames of the LD matrix to the present SNPs and return the aligned LD matrix
            dimnames(M) <- list(present, present)

            # Return the aligned LD matrix
            M
        }


        build_dataset <- function(tbl, ld, type, N, s = NULL, sdY = NULL) {
            """
            Builds a dataset for coloc analysis by standardizing column names, collapsing duplicate SNPs, filtering ambiguous SNPs, aligning the LD matrix, and returning a list containing the dataset and the number of SNPs.

            Args:
                tbl: A data frame containing the input data with standardized column names for SNPs, beta, variance of beta, MAF, effect allele, non-effect allele, and position.
                ld: A list containing the LD matrix and allele information (output from read_ld).
                type: A string indicating the type of dataset (e.g., 'quant' for quantitative traits).
                N: An integer representing the sample size.
                s: An optional numeric vector representing the standard errors of the beta estimates.
                sdY: An optional numeric value representing the standard deviation of the trait.

            Returns:
                A list containing:
                    - D: A list representing the dataset for coloc analysis, including beta, variance of beta, SNPs, positions, type, sample size, and aligned LD matrix.
                    - n: An integer representing the number of SNPs in the dataset.
                Returns NULL if there are fewer than 2 SNPs in the dataset or if the aligned LD matrix is NULL.

            """
            # Mutate the input table to create a new data frame D0 with standardized column names for SNPs, beta, variance of beta, MAF, effect allele, non-effect allele, and position
            D0 <- tbl %>% transmute(
                snp = self.standardized_variant_id_key,
                beta = self.standardized_beta_key,
                varbeta = self.standardized_var_beta_key,
                MAF = self.standardized_maf_key,
                EA = self.standardized_effect_allele_key,
                OA = self.standardized_non_effect_allele_key,
                POS = self.standardized_pos_key
            )

            # Collapse duplicate SNPs by keeping the one with the smallest p-value (largest absolute z-score)
            D0 <- .collapse_keys(D0)
            if (self$drop_ambiguous) D0 <- D0 %>% filter(!.is_ambiguous(EA, OA))
            
            # Check if there are fewer than 2 SNPs in the dataset and return NULL if so
            if (nrow(D0) < 2) return(NULL)

            # Align the LD matrix to the effect allele and subset it to the SNPs in D0
            LD <- align_ld(ld, D0$snp, setNames(D0$EA, D0$snp))
            
            # Check if the aligned LD matrix is NULL and return NULL if so
            if (is.null(LD)) return(NULL)

            # Subset D0 to only include SNPs present in the aligned LD matrix and reorder D0 to match the order of SNPs in the LD matrix
            keep <- D0$snp %in% rownames(LD)
            D0 <- D0[keep, , drop = FALSE]

            # Reorder D0 to match the order of SNPs in the LD matrix
            D0 <- D0[match(rownames(LD), D0$snp), , drop = FALSE]
            
            # Check if there are fewer than the minimum number of SNPs required for SuSiE and return a list indicating too few SNPs if so
            if (nrow(D0) < self$susie_min_snps) {
                return(list(too_few = TRUE, n = nrow(D0)))
            }

            # Create a list D containing the necessary information for the coloc analysis, including beta, variance of beta, SNPs, positions, type, sample size, and aligned LD matrix
            D <- list(
                beta = D0$beta,
                varbeta = D0$varbeta,
                snp = D0$snp,
                position = D0$POS,
                type = type,
                N = N,
                LD = LD
            )

            # Add optional parameters s and sdY to the list D if they are not NULL
            if (!is.null(s))   D$s   <- s
            if (!is.null(sdY)) D$sdY <- sdY

            # Set the MAF in the list D to 0.5 for any non-finite values in D0$MAF
            if (any(is.finite(D0$MAF))) D$MAF <- ifelse(is.finite(D0$MAF), D0$MAF, 0.5)
            
            # Return a list containing the dataset D and the number of SNPs in D0
            list(D = D, n = nrow(D0))
        }


        safe_runsusie <- function(D, label) {
            """
            Safely runs the SuSiE algorithm on a given dataset, handling errors and logging the results.

            Args:
                D: A list representing the dataset for SuSiE analysis, including beta, variance of beta, SNPs, positions, type, sample size, and aligned LD matrix.
                label: A string label for logging purposes.

            Returns:
                The result of the SuSiE analysis if successful, or NULL if there was an error

            """
            # Check if the dataset is valid for SuSiE analysis using the check_dataset function. If it fails, log the error and return NULL.
            if (!is.null(check_dataset(D, req = "LD"))) {     # NULL means OK
                .log("  [%s] check_dataset failed (LD) -> skip SuSiE", label); return(NULL)
            }

            # Attempt to run the SuSiE algorithm on the dataset, suppressing warnings and handling errors
            S <- tryCatch(
                suppressWarnings(runsusie(
                    D, coverage = self$susie_cred_coverage, max_iter = self$susie_max_iter, L = self$susie_l,
                    repeat_until_convergence = self$susie_repeat_until_converged
                )),
                error = function(e) { .log("  [%s] runsusie error: %s", label, conditionMessage(e)); NULL }
            )

            # Check if the SuSiE result is NULL and return NULL if it is
            if (is.null(S)) return(NULL)

            # Count the number of credible sets in the SuSiE result, handling any errors
            ncs <- tryCatch(length(S$sets$cs), error = function(e) 0L)
            .log("  [%s] SuSiE credible sets: %d", label, ncs %||% 0L)

            # Return NULL if there are no credible sets, indicating that there is nothing to colocalize
            if (is.null(ncs) || ncs < 1) return(NULL)

            # Return the SuSiE result
            S
        },

        run_susie_coloc <- function(self, gwas_susie_fit, qtl_susie_fit, gwas_stratum, qtl_stratum, locus, ld, gene_id, gene_name, n_overlap, n_gwas_ld, n_eqtl_ld, n_gwas_signals, n_eqtl_signals, top_gwas_variant, top_gwas_pval, lead_eqtl_id, lead_eqtl_p, full_path) {
            """
            Runs SuSiE colocalization analysis between GWAS and QTL datasets, logs the results, and writes the summary to a specified file.

            Args:
                gwas_susie_fit: The SuSiE fit result for the GWAS dataset.
                qtl_susie_fit: The SuSiE fit result for the QTL dataset.
                gwas_stratum: The stratum identifier for the GWAS dataset.
                qtl_stratum: The stratum identifier for the QTL dataset.
                locus: A list containing locus information (e.g., LOCUS_ID, CHR_STR, START, END).
                ld: A list containing LD information (e.g., status, e_same).
                gene_id: The gene identifier for the QTL dataset.
                gene_name: The gene name for the QTL dataset.
                n_overlap: The number of overlapping SNPs between the GWAS and QTL datasets.
                n_gwas_ld: The number of SNPs in LD for the GWAS dataset.
                n_eqtl_ld: The number of SNPs in LD for the QTL dataset.
                n_gwas_signals: The number of signals in the GWAS dataset.
                n_eqtl_signals: The number of signals in the QTL dataset.
                top_gwas_variant: The top variant in the GWAS dataset.
                top_gwas_pval: The p-value of the top variant in the GWAS dataset.
                lead_eqtl_id: The lead eQTL identifier.
                lead_eqtl_p: The p-value of the lead eQTL.
                full_path: The file path to write the colocalization summary results.

            Returns:
                None. The function logs the results and writes the colocalization summary to the specified file.

            """
            .log("Running SuSiE colocalization...")
            
            # Check if both GWAS and QTL SuSiE fits are not NULL before proceeding with colocalization analysis
            if (!is.null(gwas_susie_fit) && !is.null(qtl_susie_fit)) {
                
                # Run coloc.susie with the provided GWAS and QTL SuSiE fits, handling any errors
                coloc_result <- tryCatch(coloc.susie(gwas_susie_fit, qtl_susie_fit, p1 = self$coloc_priors$p1, p2 = self$coloc_priors$p2, p12 = self$coloc_priors$p12),
                                error = function(e) { .log("  coloc.susie error: %s", conditionMessage(e)); NULL })
                
                # Check if the coloc result is valid and contains a summary before proceeding to write the results
                if (!is.null(coloc_result) && !is.null(coloc_result$summary) && nrow(coloc_result$summary)) {
                    
                    # Convert the coloc summary to a data frame for easier manipulation
                    coloc_summary <- as.data.frame(coloc_result$summary)

                    # Loop through each row of the coloc summary and extract the credible set information for that row
                    for (result_index in seq_len(nrow(coloc_summary))) {

                        # Extract the credible set information for the current row of the coloc summary
                        credible_set <- .credset_for_row(coloc_result, result_index)
                        
                        # Write the colocalization summary results to the specified file path in a streaming manner, appending to the file if it already exists
                        .stream_write(tibble(
                            gwas_stratum = gwas_stratum, qtl_stratum = qtl_stratum, locus_id = locus$LOCUS_ID,
                            is_mhc = .is_mhc_locus(locus$CHR_STR, locus$START, locus$END),
                            ld_status = ld$status,
                            same_ld_panel_for_both = isTRUE(ld$e_same),
                            ld_panel_note = LD_PANEL_NOTE,
                            gene_id = gene_id, gene_name = gene_name, method = "susie",
                            idx1 = coloc_summary$idx1[result_index] %||% NA, idx2 = coloc_summary$idx2[result_index] %||% NA,
                            hit1 = as.character(coloc_summary$hit1[result_index]), hit2 = as.character(coloc_summary$hit2[result_index]),
                            nsnps = coloc_summary$nsnps[result_index] %||% NA_integer_,
                            PP.H0 = coloc_summary$PP.H0.abf[result_index], PP.H1 = coloc_summary$PP.H1.abf[result_index], PP.H2 = coloc_summary$PP.H2.abf[result_index],
                            PP.H3 = coloc_summary$PP.H3.abf[result_index], PP.H4 = coloc_summary$PP.H4.abf[result_index],
                            n_snps_overlap = n_overlap,
                            n_gwas_ld = n_gwas_ld, n_eqtl_ld = n_eqtl_ld,
                            n_gwas_signals = n_gwas_signals, n_eqtl_signals = n_eqtl_signals,
                            cred_set_size_95 = credible_set$size, top_snp_h4 = credible_set$top,
                            top_gwas_variant = top_gwas_variant, top_gwas_pval = top_gwas_pval,
                            lead_eqtl_id = lead_eqtl_id, lead_eqtl_p = lead_eqtl_p
                        ), full_path)

                        # Update the bucket counter for SuSiE pairs and set the wrote_any flag to TRUE
                        bucket["susie_pairs"] <- bucket["susie_pairs"] + 1L

                        # Set the wrote_any flag to TRUE to indicate that results have been written
                        wrote_any <- TRUE
                    }
                }
            }
            
        }
    )
)