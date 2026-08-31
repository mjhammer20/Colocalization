#!/usr/bin/env Rscript

# Imports

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(purrr)
  library(susieR)
  library(coloc)
  library(R6)
  library(argparse)
})

# ------------------------------ Helper Function Definitions -------------------------------

.load_table <- function(path){
    #
    # Loads a table from a specified file path, automatically detecting the file format based on the file extension.

    # Args:
    #     path: The file path to load the table from.

    # Returns:
    #     A data frame containing the loaded table.
    
    # Check if the file exists
    if (!file.exists(path)) {
        stop(sprintf("File not found: %s", path))
    }

    # Determine the file extension and load the table accordingly

    ext <- tools::file_ext(path)
    if (ext == "csv" || ext == "csv.gz") {
        return(readr::read_csv(path))
    } else if (ext == "tsv" || ext == "tsv.gz") {
        return(readr::read_tsv(path))
    } else if (ext == "txt" || ext == "txt.gz") {
        df <- try(readr::read_table(path))
        if (inherits(df, "try-error")) {
            df <- try(readr::read_csv(path))
            if (inherits(df, "try-error")) {
                stop(sprintf("Failed to read table from %s", path))
            }
        }
        return(df)
    } else {
        stop(sprintf("Unsupported file format: %s", ext))
    }
}


.log <- function(...) {
    
    # Logs a message with a timestamp to both the console and a log file.

    # Args:
    #     ...: The message to log, which can be formatted using sprintf-style formatting.
    
    # Returns:
    #     None. The function prints the message to the console and appends it to a log file.
    
    # Format the message with a timestamp and write it to both the console and the log file
    msg <- sprintf(paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), " | ", sprintf(...)))
    cat(msg, "\n", file = stderr())
}


.norm_chr <- function(x) {

    # Standardizes chromosome labels by converting them to uppercase and ensuring they self$manifest_left_bound_key with 'chr'.

    # Args:
    #     x: A character vector of chromosome labels (e.g., 'chr1', 'CHRX', '2').
    
    # Returns:
    #     A character vector of standardized chromosome labels (e.g., 'chr1', 'chrX', 'chr2').

    return(ifelse(grepl("^chr", x, ignore.case = TRUE), substr(x, 4, nchar(x)), x))
}


.add_variant_id <- function(chr, pos, a1, a2) {
    
    # Constructs a variant ID for a genetic variant based on chromosome, position, and alleles.

    # Args:
    #     chr: Chromosome identifier (e.g., 'chr1', 'chrX').
    #     pos: Position of the variant on the chromosome (integer).
    #     a1: First allele (string).
    #     a2: Second allele (string).

    # Returns:
    #     A string representing the canonical key for the variant in the format 'chr:pos:allele1_allele2'
    
    # Construct the canonical key using the standardized chromosome, position, and ordered alleles
    return(paste0(.norm_chr(chr), ":", as.integer(pos), ":", a1, ":", a2))
}


.is_ambiguous <- function(a1, a2) {
    
    # Checks if a pair of alleles is ambiguous (i.e., A/T or C/G).

    # Args:
    #     a1: First allele (string).
    #     a2: Second allele (string).

    # Returns:
    #     TRUE if the allele pair is ambiguous, FALSE otherwise.
    
    # Construct a string representing the allele pair in uppercase
    p <- paste0(toupper(a1), toupper(a2))

    # Check if the allele pair is one of the ambiguous pairs (A/T or C/G)
    return(p %in% c("AT", "TA", "CG", "GC"))
}


.stream_write <- function(df, path) {
    
    # Writes a data frame to a file in a streaming manner. If the file already exists, it appends the data frame to the existing file without writing column names.

    # Args:
    #     df: The data frame to be written.
    #     path: The file path to write the data frame to.

    # Returns:
    #     None. The function writes the data frame to the specified file path.    
    
    # Use data.table's fwrite function to write the data frame to the specified path
    if (!file.exists(path)) data.table::fwrite(df, path, sep = "\t", quote = FALSE)
    else data.table::fwrite(df, path, sep = "\t", quote = FALSE, append = TRUE, col.names = FALSE)
}


.collapse_keys <- function(D) {
    
    # Collapses duplicate SNPs in a data frame by keeping the one with the smallest p-value (largest absolute z-score).

    # Args:
    #     D: A data frame containing SNP information, including 'snp', 'beta', 'varbeta', 'MAF', 'EA', 'OA', and 'POS' columns.

    # Returns:
    #     A data frame with duplicate SNPs collapsed, keeping the one with the smallest p-value (largest absolute z-score).
    
    return(D %>%
        dplyr::filter(is.finite(beta), is.finite(varbeta), varbeta > 0) %>%
        dplyr::group_by(snp) %>%
        dplyr::summarise(
        beta    = sum(beta / varbeta) / sum(1 / varbeta),
        varbeta = 1 / sum(1 / varbeta),
        MAF     = suppressWarnings(mean(MAF, na.rm = TRUE)),
        EA      = dplyr::first(EA),     # effect/ALT allele defining beta sign
        OA      = dplyr::first(OA),
        POS     = dplyr::first(POS),
        .groups = "drop"
        ) %>%
        dplyr::filter(is.finite(beta), is.finite(varbeta), varbeta > 0)
    )
}


# ------------------------------ Class Definition -------------------------------

ColocalizationAnalyzer <- R6Class("ColocalizationAnalyzer",
    public = list(

        # Fields
        output_dir = NULL,
        ld_dir = NULL,
        qc_dir = NULL,
        gwas = NULL,
        gwas_path = NULL,
        gwas_target_allele_key = NULL,
        gwas_strata_key = NULL,
        gwas_sample_size = NULL,
        gwas_s = NULL,
        qtl = NULL,
        qtl_path = NULL,
        qtl_target_allele_key = NULL,
        qtl_strata_key = NULL,
        qtl_sample_size = NULL,
        qtl_sdY = NULL,
        ld_manifest = NULL,
        ld_manifest_path = NULL,
        drop_ambiguous = NULL,
        min_overlap = NULL,
        susie_min_snps = NULL,
        pp_h4_threshold = NULL,
        cred_coverage = NULL,
        susie_max_iter = NULL,
        susie_l = NULL,
        susie_repeat_until_converged = NULL,
        coloc_priors = NULL,
        manifest_locus_id_key = NULL,
        manifest_left_bound_key = NULL,
        manifest_right_bound_key = NULL,
        manifest_ld_key = NULL,
        manifest_bim_key = NULL,
        manifest_note_key = NULL,
        standardized_variant_id_key = NULL,
        standardized_chr_key = NULL,
        standardized_pos_key = NULL,
        standardized_gene_id_key = NULL,
        standardized_effect_allele_key = NULL,
        standardized_non_effect_allele_key = NULL,
        standardized_p_key = NULL,
        standardized_beta_key = NULL,
        standardized_var_beta_key = NULL,
        standardized_se_key = NULL,
        standardized_maf_key = NULL,
        bucket = NULL,

        # Initialize the class with paths and column keys
        initialize = function(args) {

            # Inputs/Outputs
            self$output_dir <- args$output_dir
            self$ld_dir <- file.path(args$ld_dir)
            self$qc_dir <- file.path(args$qc_dir)
            self$gwas_path <- file.path(args$gwas_fp)
            self$gwas_target_allele_key <- args$gwas_target_allele_key
            self$gwas_strata_key <- args$gwas_strata_key
            self$gwas_sample_size <- args$gwas_sample_size
            self$gwas_s <- args$gwas_case_fraction
            self$qtl_path <- file.path(args$qtl_fp)
            self$qtl_target_allele_key <- args$qtl_target_allele_key
            self$qtl_strata_key <- args$qtl_strata_key
            self$qtl_sample_size <- args$qtl_sample_size
            self$qtl_sdY <- args$qtl_sdY
            self$ld_manifest <- args$ld_manifest
            self$ld_manifest_path <- file.path(self$ld_dir, args$ld_manifest)

            # Parameters
            self$drop_ambiguous <- args$drop_ambiguous
            self$min_overlap <- args$min_overlap
            self$susie_min_snps <- args$susie_min_snps
            self$pp_h4_threshold <- args$pp_h4_threshold
            self$cred_coverage <- args$cred_coverage
            self$susie_max_iter <- args$susie_max_iter
            self$susie_l <- args$susie_l
            self$susie_repeat_until_converged <- args$susie_repeat_until_converged
            self$coloc_priors <- list(
                p1 = args$coloc_prior_p1,
                p2 = args$coloc_prior_p2,
                p12 = args$coloc_prior_p12
            )

            # Manifest Keys
            self$manifest_locus_id_key <- args$manifest_locus_id_key
            self$manifest_left_bound_key <- args$manifest_left_bound_key
            self$manifest_right_bound_key <- args$manifest_right_bound_key
            self$manifest_ld_key <- args$manifest_ld_key
            self$manifest_bim_key <- args$manifest_bim_key
            self$manifest_note_key <- args$manifest_note_key

            # Standardized Column Keys
            self$standardized_variant_id_key <- args$standardized_variant_id_key
            self$standardized_chr_key <- args$standardized_chr_key
            self$standardized_pos_key <- args$standardized_pos_key
            self$standardized_gene_id_key <- args$standardized_gene_id_key
            self$standardized_effect_allele_key <- args$standardized_effect_allele_key
            self$standardized_non_effect_allele_key <- args$standardized_non_effect_allele_key
            self$standardized_p_key <- args$standardized_p_key
            self$standardized_beta_key <- args$standardized_beta_key
            self$standardized_var_beta_key <- args$standardized_var_beta_key
            self$standardized_se_key <- args$standardized_se_key
            self$standardized_maf_key <- args$standardized_maf_key

            # Tracking
            self$bucket <- c(susie_pairs = 0L, abf_fallback = 0L, low = 0L)
        },

        read_ld = function(ld_path, bim_path, panel_note) {
            
            # Reads LD matrix and BIM file, aligns them, and returns a list containing the LD matrix and allele information.

            # Args:
            #     ld_path: Path to the LD matrix file.
            #     bim_path: Path to the BIM file.

            # Returns:
            #     A list containing:
            #         - M: The LD matrix with canonical key dimnames.
            #         - A1: A named vector of the allele the LD is counted on.
            #         - A2: A named vector of the other allele.
            #     Returns NULL if the LD or BIM files are missing or if there are issues with reading them
            
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

            # Create a list containing the LD matrix and allele information
            ld <- list(M = M, A1 = setNames(a1, keys), A2 = setNames(a2, keys), panel_note = panel_note)

            # Return the list containing the LD matrix and allele information
            return(ld)
        },


        align_ld = function(ld, snps, target_allele) {
            
            # Aligns the LD matrix to a target allele for each variant, subsets it to the specified SNPs, and returns the aligned LD matrix.

            # Args:
            #     ld: A list containing the LD matrix and allele information (output from read_ld).
            #     snps: A character vector of SNPs to subset the LD matrix to.
            #     target_allele: A named vector of target alleles for each SNP (names should match the SNPs) in which the datasets beta is on (ie. GWAS effect allele / GTEx ALT).
            
            # Returns:
            #     A numeric matrix representing the aligned LD matrix for the specified SNPs, with rows and columns corresponding to the SNPs in the order they appear in the `snps` vector.
            #     Returns NULL if there are fewer than 2 SNPs present or if the LD matrix is NULL.
            
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
            return(M)
        },


        build_dataset = function(tbl, ld, type, N, s = NULL, sdY = NULL, target_allele_key) {
            
            # Builds a dataset for coloc analysis by standardizing column names, collapsing duplicate SNPs, filtering ambiguous SNPs, aligning the LD matrix, and returning a list containing the dataset and the number of SNPs.

            # Args:
            #     tbl: A data frame containing the input data with standardized column names for SNPs, beta, variance of beta, MAF, effect allele, non-effect allele, and position.
            #     ld: A list containing the LD matrix and allele information (output from read_ld).
            #     type: A string indicating the type of dataset (e.g., 'quant' for quantitative traits).
            #     N: An integer representing the sample size.
            #     s: An optional numeric vector representing the standard errors of the beta estimates.
            #     sdY: An optional numeric value representing the standard deviation of the trait.

            # Returns:
            #     A list containing:
            #         - D: A list representing the dataset for coloc analysis, including beta, variance of beta, SNPs, positions, type, sample size, and aligned LD matrix.
            #         - n: An integer representing the number of SNPs in the dataset.
            #     Returns NULL if there are fewer than 2 SNPs in the dataset or if the aligned LD matrix is NULL.
            
            # Mutate the input table to create a new data frame D0 with standardized column names for SNPs, beta, variance of beta, MAF, effect allele, non-effect allele, and position
            D0 <- tbl %>% transmute(
                snp = .data[[self$standardized_variant_id_key]],
                beta = .data[[self$standardized_beta_key]],
                varbeta = .data[[self$standardized_var_beta_key]],
                MAF = .data[[self$standardized_maf_key]],
                EA = .data[[self$standardized_effect_allele_key]],
                OA = .data[[self$standardized_non_effect_allele_key]],
                POS = .data[[self$standardized_pos_key]]
            )

            # Collapse duplicate SNPs by keeping the one with the smallest p-value (largest absolute z-score)
            D0 <- .collapse_keys(D0)
            if (self$drop_ambiguous) D0 <- D0 %>% filter(!.is_ambiguous(EA, OA))
            
            # Check if there are fewer than 2 SNPs in the dataset and return NULL if so
            if (nrow(D0) < 2) return(NULL)

            # Align the LD matrix to the effect allele and subset it to the SNPs in D0
            if (target_allele_key == self$standardized_effect_allele_key) {
                LD <- self$align_ld(ld, D0$snp, setNames(D0$EA, D0$snp))
            }
            else if (target_allele_key == self$standardized_non_effect_allele_key) {
                LD <- self$align_ld(ld, D0$snp, setNames(D0$OA, D0$snp))
            }
            else {
                .log("  Unknown target allele key: %s", target_allele_key)
                return(NULL)
            }

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
            if (!is.null(s)) D$s <- s
            if (!is.null(sdY)) D$sdY <- sdY

            # Set the MAF in the list D to 0.5 for any non-finite values in D0$MAF
            if (any(is.finite(D0$MAF))) D$MAF <- ifelse(is.finite(D0$MAF), D0$MAF, 0.5)
            
            # Create a list containing the dataset D and the number of SNPs in D0
            dataset <- list(D = D, n = nrow(D0))

            # Return the list containing the dataset and the number of SNPs
            return(dataset)
        },


        credset_for_row = function(coloc_result, row_k) {
            
            # Extracts the credible set size and top SNP for a specific row in the coloc result.

            # Args:
            #     coloc_result: The result object from the coloc analysis.
            #     row_k: The index of the row to extract credible set information from.

            # Returns:
            #     A list containing:
            #         - size: The size of the credible set (number of SNPs).
            #         - top: The SNP with the highest posterior probability in the credible set.
            
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
            credible_set <- cumsum(v[o]); w <- which(credible_set >= self$cred_coverage)[1]
            
            # If no credible set is found, set the size to the total number of SNPs
            if (is.na(w)) w <- length(o)
            out$size <- as.integer(w)
            out$top  <- as.character(rr$snp[o][1])
            
            # Return the output containing the credible set size and top SNP
            return(out)
        },


        safe_runsusie = function(D, label) {
            
            # Safely runs the SuSiE algorithm on a given dataset, handling errors and logging the results.

            # Args:
            #     D: A list representing the dataset for SuSiE analysis, including beta, variance of beta, SNPs, positions, type, sample size, and aligned LD matrix.
            #     label: A string label for logging purposes.

            # Returns:
            #     The result of the SuSiE analysis if successful, or NULL if there was an error
            
            # Check if the dataset is valid for SuSiE analysis using the check_dataset function. If it fails, log the error and return NULL.
            if (!is.null(check_dataset(D, req = "LD"))) {     # NULL means OK
                .log("  [%s] check_dataset failed (LD) -> skip SuSiE", label); return(NULL)
            }

            # Attempt to run the SuSiE algorithm on the dataset, suppressing warnings and handling errors
            S <- tryCatch(
                suppressWarnings(runsusie(
                    D, coverage = self$cred_coverage, max_iter = self$susie_max_iter, L = self$susie_l,
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
            if (is.null(ncs) || ncs < 1) {
                .log("  [%s] SuSiE has no credible sets -> skip coloc", label)
                return(NULL)
            }

            # Return the SuSiE result
            return(S)
        },


        run_susie_coloc = function(gwas_susie_fit, qtl_susie_fit, gwas_stratum, qtl_stratum, locus, ld, gene_id, n_overlap, n_gwas_ld, n_qtl_ld, n_gwas_signals, n_qtl_signals, top_gwas_variant, top_gwas_pval, lead_eqtl_id, lead_eqtl_p, coloc_full_fp) {
            
            # Runs SuSiE colocalization analysis between GWAS and QTL datasets, logs the results, and writes the summary to a specified file.

            # Args:
            #     gwas_susie_fit: The SuSiE fit result for the GWAS dataset.
            #     qtl_susie_fit: The SuSiE fit result for the QTL dataset.
            #     gwas_stratum: The stratum identifier for the GWAS dataset.
            #     qtl_stratum: The stratum identifier for the QTL dataset.
            #     locus: A list containing locus information (e.g., LOCUS_ID, self$standardized_chr_key, self$manifest_left_bound_key, self$manifest_right_bound_key).
            #     ld: A list containing LD information.
            #     gene_id: The gene identifier for the QTL dataset.
            #     n_overlap: The number of overlapping SNPs between the GWAS and QTL datasets.
            #     n_gwas_ld: The number of SNPs in LD for the GWAS dataset.
            #     n_qtl_ld: The number of SNPs in LD for the QTL dataset.
            #     n_gwas_signals: The number of signals in the GWAS dataset.
            #     n_qtl_signals: The number of signals in the QTL dataset.
            #     top_gwas_variant: The top variant in the GWAS dataset.
            #     top_gwas_pval: The p-value of the top variant in the GWAS dataset.
            #     lead_eqtl_id: The lead QTL identifier.
            #     lead_eqtl_p: The p-value of the lead eQTL.
            #     coloc_full_fp: The file path to write the colocalization summary results.

            # Returns:
            #     A logical value indicating whether any results were written to the specified file (TRUE if results were written, FALSE otherwise).

            # Initialize a flag to track whether any results were written to the specified file
            wrote_any <- FALSE 
            
            # Check if both GWAS and QTL SuSiE fits are not NULL before proceeding with colocalization analysis
            if (!is.null(gwas_susie_fit) && !is.null(qtl_susie_fit)) {

                .log("Running SuSiE colocalization...")
                
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
                        credible_set <- self$credset_for_row(coloc_result, result_index)
                        
                        # Write the colocalization summary results to the specified file path in a streaming manner, appending to the file if it already exists
                        .stream_write(tibble(
                            !!self$gwas_strata_key := gwas_stratum,
                            !!self$qtl_strata_key := qtl_stratum,
                            !!self$manifest_locus_id_key := locus[[self$manifest_locus_id_key]],
                            ld_panel_note = locus[[self$manifest_note_key]],
                            gene_id = gene_id,
                            method = "susie",
                            idx1 = coloc_summary$idx1[result_index] %||% NA,
                            idx2 = coloc_summary$idx2[result_index] %||% NA,
                            hit1 = as.character(coloc_summary$hit1[result_index]),
                            hit2 = as.character(coloc_summary$hit2[result_index]),
                            nsnps = coloc_summary$nsnps[result_index] %||% NA_integer_,
                            PP.H0 = coloc_summary$PP.H0.abf[result_index],
                            PP.H1 = coloc_summary$PP.H1.abf[result_index],
                            PP.H2 = coloc_summary$PP.H2.abf[result_index],
                            PP.H3 = coloc_summary$PP.H3.abf[result_index],
                            PP.H4 = coloc_summary$PP.H4.abf[result_index],
                            n_snps_overlap = n_overlap,
                            n_gwas_ld = n_gwas_ld,
                            n_qtl_ld = n_qtl_ld,
                            n_gwas_signals = n_gwas_signals,
                            n_qtl_signals = n_qtl_signals,
                            cred_set_size_95 = credible_set$size,
                            top_snp_h4 = credible_set$top,
                            top_gwas_variant = top_gwas_variant,
                            top_gwas_pval = top_gwas_pval,
                            lead_eqtl_id = lead_eqtl_id,
                            lead_eqtl_p = lead_eqtl_p
                        ), coloc_full_fp)

                        # Update the bucket counter for SuSiE pairs and set the wrote_any flag to TRUE
                        self$bucket["susie_pairs"] <- self$bucket["susie_pairs"] + 1L

                        # Set the wrote_any flag to TRUE to indicate that results have been written
                        wrote_any <- TRUE
                    }
                }
            }

                # Return the wrote_any flag to indicate whether any results were written to the specified file
                return(wrote_any)
        },


        run_coloc_abf_fallback = function(gwas_table, qtl_table, gwas_sample_size, gwas_s, qtl_sample_size, qtl_sdY) {
            
            # Runs coloc.abf analysis as a fallback method for colocalization between GWAS and QTL datasets, filtering for shared SNPs and handling errors.

            # Args:
            #     gwas_table: A data frame containing the GWAS dataset with standardized column names for SNPs, beta, variance of beta, and MAF.
            #     qtl_table: A data frame containing the QTL dataset with standardized column names for SNPs, beta, variance of beta, and MAF.
            #     gwas_sample_size: An integer representing the sample size for the GWAS dataset.
            #     gwas_s: An optional numeric vector representing the standard errors of the beta estimates for the GWAS dataset.
            #     qtl_sample_size: An integer representing the sample size for the QTL dataset.
            #     qtl_sdY: An optional numeric value representing the standard deviation of the trait for the QTL dataset.
            #     coloc_full_fp: The file path to write the colocalization summary results.

            # Returns:
            #     A list containing:
            #         - n: The number of shared SNPs between the GWAS and QTL datasets.
            #         - PP: A numeric vector of posterior probabilities for hypotheses H0 to H4.
            #         - top: The SNP with the highest posterior probability for hypothesis H4.
            #     Returns NULL if there are fewer than the minimum overlap required for ABF analysis or if there are issues with the coloc.abf analysis.
            
            .log("Running coloc.abf fallback...")

            # Join the GWAS and QTL tables on the SNP column, filtering for finite beta and variance values
            shared <- inner_join(
                gwas_table %>% transmute(snp = .data[[self$standardized_variant_id_key]], gwas_beta = .data[[self$standardized_beta_key]], gwas_variance = .data[[self$standardized_se_key]]^2, gwas_maf = .data[[self$standardized_maf_key]]),
                qtl_table %>% transmute(snp = .data[[self$standardized_variant_id_key]], qtl_beta = .data[[self$standardized_beta_key]], qtl_variance = .data[[self$standardized_se_key]]^2, qtl_maf = .data[[self$standardized_maf_key]]),
                by = "snp"
            ) %>% filter(is.finite(gwas_beta), gwas_variance > 0, is.finite(qtl_beta), qtl_variance > 0)

            # Check if the number of shared SNPs is less than the minimum overlap required for ABF analysis, and return NULL if so
            if (nrow(shared) < self$min_overlap) return(NULL)

            # Collapse duplicate SNPs by calculating the weighted average of beta and variance for both GWAS and QTL datasets
            if (any(duplicated(shared$snp))) {
                shared <- shared %>% group_by(snp) %>% summarise(
                gwas_beta = sum(gwas_beta/gwas_variance)/sum(1/gwas_variance), gwas_variance = 1/sum(1/gwas_variance),
                qtl_beta = sum(qtl_beta/qtl_variance)/sum(1/qtl_variance), qtl_variance = 1/sum(1/qtl_variance),
                gwas_maf = mean(gwas_maf, na.rm = TRUE), qtl_maf = mean(qtl_maf, na.rm = TRUE), .groups = "drop")
            }

            # Prepare the datasets for coloc.abf analysis, including beta, variance of beta, sample size, case fraction, type, and MAF for both GWAS and QTL datasets
            d1 <- list(snp = shared$snp, beta = shared$gwas_beta, varbeta = shared$gwas_variance,
                        N = gwas_sample_size, s = gwas_s, type = "cc",
                        MAF = ifelse(is.finite(shared$gwas_maf), shared$gwas_maf, 0.5))
            d2 <- list(snp = shared$snp, beta = shared$qtl_beta, varbeta = shared$qtl_variance,
                        N = qtl_sample_size, type = "quant", sdY = qtl_sdY,
                        MAF = ifelse(is.finite(shared$qtl_maf), shared$qtl_maf, 0.5))

            # Run coloc.abf analysis with the prepared datasets and specified priors, handling any errors
            coloc_result <- tryCatch(coloc.abf(d1, d2, p1 = self$coloc_priors$p1, p2 = self$coloc_priors$p2, p12 = self$coloc_priors$p12),
                            error = function(e) NULL)

            # Check if the coloc result is NULL and return NULL if it is
            if (is.null(coloc_result)) return(NULL)

            # Extract the summary from the coloc result and initialize the top SNP variable
            s <- coloc_result$summary

            # Determine the top SNP with the highest posterior probability for hypothesis H4, if available
            top <- NA_character_
            if (!is.null(coloc_result$results) && "SNP.PP.H4" %in% names(coloc_result$results)) {
                tr <- coloc_result$results %>% arrange(desc(SNP.PP.H4)) %>% slice(1)
                if (nrow(tr)) top <- as.character(tr$snp)
            }

            # Create a list containing the number of shared SNPs, posterior probabilities for hypotheses H0 to H4, and the top SNP
            coloc_abf_result <- list(n = nrow(shared),
                PP = c(s["PP.H0.abf"], s["PP.H1.abf"], s["PP.H2.abf"], s["PP.H3.abf"], s["PP.H4.abf"]),
                top = top)

            # Set wrote_any to TRUE to indicate that results have been written
            wrote_any <- TRUE

            # Return the wrote_any flag to indicate whether any results were written to the specified file
            return(coloc_abf_result) 
        },


        process_locus = function(locus, gwas_sub, gwas_stratum, gwas_sample_size, gwas_s, qtl_sub, qtl_stratum, qtl_sample_size, qtl_sdY, susie_low_overlap_fp, coloc_full_fp){
            
            # Processes a single locus for colocalization analysis between GWAS and QTL datasets, filtering variants, running SuSiE, and recording results.

            # Args:
            #     locus: A list containing locus information (e.g., LOCUS_ID, self$standardized_chr_key, self$manifest_left_bound_key, self$manifest_right_bound_key).
            #     gwas_sub: A data frame containing the GWAS dataset for the current stratum.
            #     gwas_stratum: The stratum identifier for the GWAS dataset.
            #     gwas_sample_size: An integer representing the sample size for the GWAS dataset.
            #     gwas_s: An optional numeric vector representing the standard errors of the beta estimates for the GWAS dataset.
            #     qtl_sub: A data frame containing the QTL dataset for the current stratum.
            #     qtl_stratum: The stratum identifier for the QTL dataset.
            #     qtl_sample_size: An integer representing the sample size for the QTL dataset.
            #     qtl_sdY: An optional numeric value representing the standard deviation of the trait for the QTL dataset.
            #     susie_low_overlap_fp: The file path to write results with low overlap between GWAS and QTL datasets.
            #     coloc_full_fp: The file path to write full colocalization results.

            # Returns:
            #     locus_result: A tibble containing the results of the colocalization analysis for the current locus, including stratum identifiers, locus ID, gene ID, number of overlapping SNPs, number of signals, and posterior probabilities for hypotheses H0 to H4.
            
            .log("Processing locus %s for GWAS stratum %s and QTL stratum %s", locus[[self$manifest_locus_id_key]], gwas_stratum, qtl_stratum)

            # Filter the GWAS and QTL datasets to include only variants within the specified locus boundaries (chromosome and position range)
            locus_gwas <- gwas_sub %>% filter(
                .data[[self$standardized_chr_key]] == locus[[self$standardized_chr_key]],
                .data[[self$standardized_pos_key]] >= locus[[self$manifest_left_bound_key]],
                .data[[self$standardized_pos_key]] <= locus[[self$manifest_right_bound_key]]
            )
            locus_qtl <- qtl_sub %>% filter(
                .data[[self$standardized_chr_key]] == locus[[self$standardized_chr_key]],
                .data[[self$standardized_pos_key]] >= locus[[self$manifest_left_bound_key]],
                .data[[self$standardized_pos_key]] <= locus[[self$manifest_right_bound_key]]
            )

            # Count the number of variants in the filtered GWAS and QTL datasets, skipping the locus if either dataset has no variants
            n_gwas = nrow(locus_gwas)
            n_qtl = nrow(locus_qtl)
            .log("%s | %s | %s: GWAS=%s QTL=%s", gwas_stratum, qtl_stratum, locus[[self$manifest_locus_id_key]], format(n_gwas, big.mark = ","), format(n_qtl, big.mark = ","))
            if (!n_gwas || !n_qtl) {
                .log("Skipping locus %s: insufficient variants (GWAS=%s, QTL=%s)", 
                    locus[[self$manifest_locus_id_key]], n_gwas, n_qtl)
                return(NULL)
            }

            # Load the LD matrix for the locus using the read_ld function
            ld <- self$read_ld(
                ld_path  = file.path(self$ld_dir, locus[[self$manifest_ld_key]]),
                bim_path = file.path(self$ld_dir, locus[[self$manifest_bim_key]]),
                panel_note = locus[[self$manifest_note_key]])

            # Build GWAS dataset for SuSiE analysis
            gwas_dataset <- self$build_dataset(locus_gwas, ld, type = "cc", N = gwas_sample_size, s = gwas_s, target_allele_key = self$gwas_target_allele_key)

            # Initialize variables for SuSiE fits, LD counts, and signal counts
            gwas_susie_fit <- NULL
            n_gwas_ld <- NA_integer_
            n_gwas_signals <- 0L

            # Run SuSiE for GWAS dataset
            if (!is.null(gwas_dataset) && is.null(gwas_dataset$too_few) && !is.null(gwas_dataset$D)) {
                n_gwas_ld <- gwas_dataset$n
                gwas_susie_fit <- self$safe_runsusie(gwas_dataset$D, sprintf("%s GWAS %s", gwas_stratum, locus[[self$manifest_locus_id_key]]))
                n_gwas_signals <- tryCatch(length(gwas_susie_fit$sets$cs), error = function(e) 0L) %||% 0L
            }

            # Record the lead GWAS variant
            top_gwas <- locus_gwas %>% arrange(.data[[self$standardized_p_key]]) %>% slice(1)
            top_gwas_variant <- top_gwas[[self$standardized_variant_id_key]] %||% NA_character_
            top_gwas_pval <- top_gwas[[self$standardized_p_key]] %||% NA_real_

            best_overlap_locus <- 0L

            # Extract gene ids from the filtered QTL dataset
            gene_ids <- unique(locus_qtl[[self$standardized_gene_id_key]])

            # Loop through each gene id and perform susie colocalization analysis
            for (gene_id in gene_ids) {
                
                # Filter the QTL dataset for the current gene id
                gene_qtl  <- locus_qtl %>% filter(.data[[self$standardized_gene_id_key]] == gene_id)

                # Check GWAS and QTL variant overlap
                n_overlap <- length(intersect(gwas_sub[[self$standardized_variant_id_key]], gene_qtl[[self$standardized_variant_id_key]]))
                if (n_overlap < self$min_overlap) {
                    .stream_write(tibble(
                        !!self$gwas_strata_key := gwas_stratum,
                        !!self$qtl_strata_key := qtl_stratum,
                        !!self$manifest_locus_id_key := locus[[self$manifest_locus_id_key]],
                        ld_panel_note = locus[[self$manifest_note_key]],
                        gene_id = gene_id,
                        n_overlap = n_overlap,
                        reason = paste0("overlap<", self$min_overlap)
                    ), susie_low_overlap_fp)
                    self$bucket["low"] <- self$bucket["low"] + 1L
                    next
                }

                # Record the lead QTL variant
                lead_qtl    <- gene_qtl %>% filter(is.finite(.data[[self$standardized_p_key]])) %>% arrange(.data[[self$standardized_p_key]]) %>% slice(1)
                lead_qtl_id <- lead_qtl[[self$standardized_variant_id_key]][1] %||% NA_character_
                lead_qtl_p  <- lead_qtl[[self$standardized_p_key]][1] %||% NA_real_

                # Build QTL dataset for SuSiE analysis
                gene_qtl_dataset <- self$build_dataset(gene_qtl, ld, type = "quant", N = qtl_sample_size, sdY = qtl_sdY, target_allele_key = self$qtl_target_allele_key)

                # Initialize variables for QTL SuSiE fit, LD counts, and signal counts
                qtl_susie_fit <- NULL
                n_qtl_ld <- NA_integer_
                n_qtl_signals <- 0L

                if (!is.null(gene_qtl_dataset) && is.null(gene_qtl_dataset$too_few) && !is.null(gene_qtl_dataset$D)) {
                    n_qtl_ld <- gene_qtl_dataset$n
                    qtl_susie_fit <- self$safe_runsusie(gene_qtl_dataset$D, sprintf("%s QTL %s/%s", gwas_stratum, locus[[self$manifest_locus_id_key]], gene_id))
                    n_qtl_signals <- tryCatch(length(qtl_susie_fit$sets$cs), error = function(e) 0L) %||% 0L
                }

                # Initialize wrote_any flag to FALSE before running SuSiE colocalization
                wrote_any <- FALSE

                # Run SuSiE colocalization analysis and write results to the specified file path
                wrote_any <- self$run_susie_coloc(
                    gwas_susie_fit,
                    qtl_susie_fit,
                    gwas_stratum,
                    qtl_stratum,
                    locus,
                    ld,
                    gene_id,
                    n_overlap,
                    n_gwas_ld,
                    n_qtl_ld,
                    n_gwas_signals,
                    n_qtl_signals,
                    top_gwas_variant,
                    top_gwas_pval,
                    lead_qtl_id,
                    lead_qtl_p,
                    coloc_full_fp
                )

                # If no results were written, run coloc.abf fallback and write results to the specified file path
                if (!wrote_any) {
                    coloc_abf_result <- self$run_coloc_abf_fallback(gwas_sub, gene_qtl, gwas_sample_size, gwas_s, qtl_sample_size, qtl_sdY)
                    
                    if (!is.null(coloc_abf_result)) {
                        # Write the coloc.abf fallback results to the specified file path in a streaming manner, appending to the file if it already exists
                        .stream_write(tibble(
                                !!self$gwas_strata_key := gwas_stratum,
                                !!self$qtl_strata_key := qtl_stratum,
                                !!self$manifest_locus_id_key := locus[[self$manifest_locus_id_key]],
                                ld_panel_note = locus[[self$manifest_note_key]],
                                gene_id = gene_id,
                                method = "coloc.abf fallback",
                                idx1 = NA_integer_,
                                idx2 = NA_integer_,
                                hit1 = NA_character_,
                                hit2 = NA_character_,
                                nsnps = coloc_abf_result$n,
                                PP.H0 = coloc_abf_result$PP[1],
                                PP.H1 = coloc_abf_result$PP[2],
                                PP.H2 = coloc_abf_result$PP[3],
                                PP.H3 = coloc_abf_result$PP[4],
                                PP.H4 = coloc_abf_result$PP[5],
                                n_snps_overlap = n_overlap,
                                n_gwas_ld = n_gwas_ld,
                                n_qtl_ld = n_qtl_ld,
                                n_gwas_signals = n_gwas_signals,
                                n_qtl_signals = n_qtl_signals,
                                cred_set_size_95 = NA_integer_,
                                top_snp_h4 = coloc_abf_result$top,
                                top_gwas_variant = top_gwas_variant,
                                top_gwas_pval = top_gwas_pval,
                                lead_qtl_id = lead_qtl_id,
                                lead_qtl_p = lead_qtl_p,
                            ), coloc_full_fp)

                            # Set wrote_any to TRUE to indicate that results have been written for coloc.abf fallback
                            wrote_any <- TRUE

                            # Update the bucket counter for coloc.abf fallback and log the successful writing of results
                            self$bucket["abf_fallback"] <- self$bucket["abf_fallback"] + 1L

                            # Log the successful writing of coloc.abf fallback results for the current gene and locus
                            .log("Coloc.abf fallback results written for gene %s in locus %s", gene_id, locus[[self$manifest_locus_id_key]])
                    } else {
                        # Log that no colocalization results were written for the current gene and locus
                        .log("No colocalization results written for gene %s in locus %s", gene_id, locus[[self$manifest_locus_id_key]])
                    }
                } else {
                    .log("SuSiE colocalization results written for gene %s in locus %s", gene_id, locus[[self$manifest_locus_id_key]])
                }

                # Update the best_overlap_locus variable if the current n_overlap is greater than the previous best_overlap_locus
                if (n_overlap > best_overlap_locus) best_overlap_locus <- n_overlap
            }

            # Create a tibble summarizing the locus results, including locus ID, number of GWAS and QTL signals, and best overlap
            locus_result <- tibble(
                !!self$manifest_locus_id_key := locus[[self$manifest_locus_id_key]],
                ld_panel_note = locus[[self$manifest_note_key]],
                n_gwas = n_gwas,
                n_qtl = n_qtl,
                n_gwas_ld = n_gwas_ld %||% NA_integer_,
                n_gwas_signals = n_gwas_signals,
                best_overlap = best_overlap_locus
            )
            
            # Return the locus result tibble summarizing the results of the colocalization analysis for the current locus
            .log("Finished processing locus %s for GWAS stratum %s and QTL stratum %s", locus[[self$manifest_locus_id_key]], gwas_stratum, qtl_stratum)
            return(locus_result)
        },


        run_coloc_analysis = function() {

            # Runs colocalization analysis between GWAS and QTL datasets for specified loci, processing each locus and recording results.

            # Args:
            #     gwas_data: A data frame containing the GWAS dataset with standardized column names for SNPs, beta, variance of beta, and MAF.
            #     qtl_data: A data frame containing the QTL dataset with standardized column names for SNPs, beta, variance of beta, and MAF.
            #     loci_data: A data frame containing locus information, including locus ID, chromosome, left bound, right bound, and LD file paths.

            # Returns:
            #     A tibble summarizing the results of the colocalization analysis across all loci and strata.

            # Load GWAS, QTL, and loci data from the specified file paths
            .log("Loading GWAS Summary Statistics...")
            gwas_data <- .load_table(self$gwas_path)
            .log("Loading QTL Summary Statistics...")
            qtl_data <- .load_table(self$qtl_path)
            .log("Loading LD Manifest...")
            loci_data <- .load_table(self$ld_manifest_path)

            # Add strata columns and fill values if necessary for GWAS and QTL datasets
            if (is.null(self$gwas_strata_key)) {
                self$gwas_strata_key <- "GWAS_Strata"
                gwas_data[[self$gwas_strata_key]] <- "Full"
            } else if (!self$gwas_strata_key %in% colnames(gwas_data)) {
                gwas_data[[self$gwas_strata_key]] <- "Full"
            } else if (gwas_data[[self$gwas_strata_key]] %>% is.na() %>% any()) {
                gwas_data[[self$gwas_strata_key]] <- gwas_data[[self$gwas_strata_key]] %>% replace_na("Full")
            }

            if (is.null(self$qtl_strata_key)) {
                self$qtl_strata_key <- "QTL_Strata"
                qtl_data[[self$qtl_strata_key]] <- "Bulk"
            } else if (!self$qtl_strata_key %in% colnames(qtl_data)) {
                qtl_data[[self$qtl_strata_key]] <- "Bulk"
            } else if (qtl_data[[self$qtl_strata_key]] %>% is.na() %>% any()) {
                qtl_data[[self$qtl_strata_key]] <- qtl_data[[self$qtl_strata_key]] %>% replace_na("Bulk")
            }

            # Initialize an empty list to store locus results
            locus_results <- list()            

            # Extract unique strata identifiers from the GWAS and QTL datasets
            gwas_strata <- unique(gwas_data[[self$gwas_strata_key]])
            qtl_strata <- unique(qtl_data[[self$qtl_strata_key]])

            # Loop through each combination of GWAS and QTL strata
            for (gwas_stratum in gwas_strata) {
                susie_low_overlap_fp <- file.path(self$qc_dir, sprintf("%s_coloc_susie_low_overlap.tsv", gwas_stratum))
                coloc_full_fp <- file.path(self$output_dir, sprintf("%s_coloc_full.tsv", gwas_stratum))
                for (qtl_stratum in qtl_strata) {
                    .log("Running colocalization analysis for GWAS stratum %s and QTL stratum %s", gwas_stratum, qtl_stratum)

                    # Filter the GWAS and QTL datasets for the current strata
                    gwas_sub <- gwas_data %>% filter(!!sym(self$gwas_strata_key) == gwas_stratum)
                    qtl_sub <- qtl_data %>% filter(!!sym(self$qtl_strata_key) == qtl_stratum)
                    
                    # Loop through each locus and process it
                    for (i in seq_len(nrow(loci_data))) {
                        locus <- loci_data[i, ]
                        locus_result <- self$process_locus(locus, gwas_sub, gwas_stratum, self$gwas_sample_size, self$gwas_s,
                                                            qtl_sub, qtl_stratum, self$qtl_sample_size, self$qtl_sdY, susie_low_overlap_fp, coloc_full_fp)
                        if (!is.null(locus_result)) {
                            locus_results <- append(locus_results, list(locus_result))
                        }
                    }
                }

                # Write summaries
                full_df <- if (file.exists(coloc_full_fp)) fread(coloc_full_fp) else data.table()
                strong  <- if (nrow(full_df)) as_tibble(full_df) %>%
                    filter(!is.na(PP.H4), PP.H4 >= self$pp_h4_threshold, n_snps_overlap >= self$min_overlap) else tibble()
                write_tsv(strong, file.path(self$qc_dir, paste0(tolower(gwas_stratum), "_coloc_susie_strong.tsv")))

                locus_tbl <- if (length(locus_results)) bind_rows(locus_results) else tibble()
                write_tsv(locus_tbl %>% arrange(self$manifest_locus_id_key),
                            file.path(self$qc_dir, paste0(tolower(gwas_stratum), "_locus_overlap_summary.tsv")))

                report_path <- file.path(self$output_dir, paste0(tolower(gwas_stratum), "_summary_report.txt"))
                report <- c(
                    sprintf("coloc.susie summary report — %s", gwas_stratum),
                    sprintf("Generated: %s", format(Sys.time())),
                    sprintf("%s: %s", self$gwas_strata_key, gwas_stratum),
                    "",
                    sprintf("Signal pairs written (susie): %d", self$bucket["susie_pairs"]),
                    sprintf("Loci processed via abf fallback: %d", self$bucket["abf_fallback"]),
                    sprintf("Loci w/ overlap below %d: %d", self$min_overlap, self$bucket["low"]),
                    sprintf("Strong (PP.H4 >= %.2f): %d rows", self$pp_h4_threshold, nrow(strong)),
                    "",
                    sprintf("Priors: p1 = %g, p2 = %g, p12 = %g", self$coloc_priors$p1, self$coloc_priors$p2, self$coloc_priors$p12),
                    sprintf("SuSiE: coverage = %.2f, max_iter = %d, L = %d", self$cred_coverage, self$susie_max_iter, self$susie_l),
                    sprintf("MIN OVERLAP = %d, MIN SNPS SUSIE = %d, DROP AMBIGUOUS = %s", self$min_overlap, self$susie_min_snps, self$drop_ambiguous),
                    sprintf("QTL sdY = %g", self$qtl_sdY),
                    sprintf("LD panel note: %s", unique(loci_data[[self$manifest_note_key]])[1])
                )
                writeLines(report, report_path)
                .log("%s: Completed. susie_pairs = %d, abf = %d, low = %d, strong = %d",
                    gwas_stratum, self$bucket["susie_pairs"], self$bucket["abf_fallback"], self$bucket["low"], nrow(strong))
                invisible(list(full = full_df, strong = strong, locus = locus_tbl, report = report_path))
            }
        }
    )
) 

# ------------------------------------- Command Line Interface -------------------------------------

if (!interactive()) {
    
    parser <- ArgumentParser(description = "Colocalization analysis using SuSiE and coloc.abf")
    parser$add_argument("--output_dir", required = TRUE, help = "Path to the directory to save output files")
    parser$add_argument("--ld_dir", required = TRUE, help = "Path to the directory containing LD matrices and BIM files.")
    parser$add_argument("--qc_dir", required = TRUE, help = "Path to the directory containing QC files.")
    parser$add_argument("--gwas_fp", required = TRUE, help = "Path to the standardized GWAS summary statistics file")
    parser$add_argument("--gwas_target_allele_key", default = "EFFECT", help = "Column name for the target allele in the GWAS summary statistics file (optional). Default is 'EFFECT'.")
    parser$add_argument("--gwas_strata_key", default = NULL, help = "Column name for the GWAS strata in the GWAS summary statistics file (optional). Default is 'STRATUM'.")
    parser$add_argument("--gwas_sample_size", type = "numeric", required = TRUE, help = "Sample size for the GWAS dataset")
    parser$add_argument("--gwas_case_fraction", type = "numeric", required = TRUE, help = "Case fraction for the GWAS dataset.")
    parser$add_argument("--qtl_fp", required = TRUE, help = "Path to the standardized QTL summary statistics file")
    parser$add_argument("--qtl_target_allele_key", default = "EFFECT", help = "Column name for the target allele in the QTL summary statistics file (optional). Default is 'EFFECT'.")
    parser$add_argument("--qtl_strata_key", default = NULL, help = "Column name for the QTL strata in the QTL summary statistics file (optional). Default is 'STRATUM'.")
    parser$add_argument("--qtl_sample_size", type = "numeric", required = TRUE, help = "Sample size for the QTL dataset")
    parser$add_argument("--qtl_sdY", type = "numeric", default = 1, help = "Standard deviation of the trait for the QTL dataset (optional). Default is 1.")
    parser$add_argument("--ld_manifest", default = "ld_manifest.tsv", help = "File name of the ld manifest file containing locus information (optional). Default is 'ld_manifest.tsv'.")
    parser$add_argument("--drop_ambiguous", default = TRUE, help = "Drop ambiguous SNPs (A/T or C/G) from the analysis (optional). Default is TRUE.")
    parser$add_argument("--min_overlap", type = "numeric", default = 50, help = "Minimum number of overlapping SNPs required for colocalization analysis (optional). Default is 50.")
    parser$add_argument("--cred_coverage", type = "numeric", default = 0.95, help = "Credible set coverage for SuSiE analysis (optional). Default is 0.95.")
    parser$add_argument("--pp_h4_threshold", type = "numeric", default = 0.80, help = "Posterior probability threshold for strong colocalization (optional). Default is 0.80.")
    parser$add_argument("--susie_min_snps", type = "numeric", default = 50, help = "Minimum number of SNPs required for SuSiE analysis (optional). Default is 50.")
    parser$add_argument("--susie_max_iter", type = "numeric", default = 100, help = "Maximum number of iterations for SuSiE analysis (optional). Default is 100.")
    parser$add_argument("--susie_l", type = "numeric", default = 10, help = "Maximum number of causal variants for SuSiE analysis (optional). Default is 10.")
    parser$add_argument("--susie_repeat_until_converged", default = TRUE, help = "Repeat SuSiE until convergence (optional). Default is TRUE.")
    parser$add_argument("--coloc_prior_p1", type = "numeric", default = 1e-4, help = "Prior probability for hypothesis H1 (optional). Default is 1e-4.")
    parser$add_argument("--coloc_prior_p2", type = "numeric", default = 1e-4, help = "Prior probability for hypothesis H2 (optional). Default is 1e-4.")
    parser$add_argument("--coloc_prior_p12", type = "numeric", default = 1e-5, help = "Prior probability for hypothesis H4 (optional). Default is 1e-5.")
    parser$add_argument("--manifest_locus_id_key", default = "LOCUS_ID", help = "Column name for the locus ID in the LD manifest file (optional). Default is 'LOCUS_ID'.")
    parser$add_argument("--manifest_left_bound_key", default = "LEFT_500KB", help = "Column name for the left bound of the locus in the LD manifest file (optional). Default is 'LEFT_500KB'.")
    parser$add_argument("--manifest_right_bound_key", default = "RIGHT_500KB", help = "Column name for the right bound of the locus in the LD manifest file (optional). Default is 'RIGHT_500KB'.")
    parser$add_argument("--manifest_ld_key", default = "LD", help = "Column name for the LD matrix file path in the LD manifest file (optional). Default is 'LD'.")
    parser$add_argument("--manifest_bim_key", default = "BIM", help = "Column name for the BIM file path in the LD manifest file (optional). Default is 'BIM'.")
    parser$add_argument("--manifest_note_key", default = "note", help = "Column name for the LD panel note in the LD manifest file (optional). Default is 'PANEL_NOTE'.")
    parser$add_argument("--standardized_variant_id_key", default = "VAR", help = "Column name for the standardized variant ID in the GWAS and QTL summary statistics files (optional). Default is 'VAR'.")
    parser$add_argument("--standardized_chr_key", default = "CHR", help = "Column name for the standardized chromosome in the GWAS and QTL summary statistics files (optional). Default is 'CHR'.")
    parser$add_argument("--standardized_pos_key", default = "BP", help = "Column name for the standardized position in the GWAS and QTL summary statistics files (optional). Default is 'POS'.")
    parser$add_argument("--standardized_gene_id_key", default = "GENE", help = "Column name for the standardized gene ID in the QTL summary statistics file (optional). Default is 'GENE_ID'.")
    parser$add_argument("--standardized_non_effect_allele_key", default = "NON_EFFECT", help = "Column name for the standardized reference allele in the GWAS and QTL summary statistics files (optional). Default is 'NON_EFFECT'.")
    parser$add_argument("--standardized_effect_allele_key", default = "EFFECT", help = "Column name for the standardized effect allele in the GWAS and QTL summary statistics files (optional). Default is 'EFFECT'.")
    parser$add_argument("--standardized_p_key", default = "P", help = "Column name for the standardized p-value in the GWAS and QTL summary statistics files (optional). Default is 'P'.")
    parser$add_argument("--standardized_beta_key", default = "BETA", help = "Column name for the standardized beta in the GWAS and QTL summary statistics files (optional). Default is 'BETA'.")
    parser$add_argument("--standardized_var_beta_key", default = "VARBETA", help = "Column name for the standardized variance of beta in the GWAS and QTL summary statistics files (optional). Default is 'VAR_BETA'.")
    parser$add_argument("--standardized_se_key", default = "SE", help = "Column name for the standardized standard error in the GWAS and QTL summary statistics files (optional). Default is 'SE'.")
    parser$add_argument("--standardized_maf_key", default = "MAF", help = "Column name for the standardized minor allele frequency in the GWAS and QTL summary statistics files (optional). Default is 'MAF'.")


    # Create an instance of the ColocalizationAnalyzer class and run the colocalization analysis
    args <- parser$parse_args()
    analyzer <- ColocalizationAnalyzer$new(args)
    analyzer$run_coloc_analysis()
}