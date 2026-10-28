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
  library(ggplot2)
  library(R.utils)
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
    msg <- sprintf(paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), " | ", "%s"), sprintf(...))
    cat(msg, "\n", file = stderr())
    flush(stderr())
    flush(stdout())
}


.norm_chr <- function(x) {

    # Standardizes chromosome labels by converting them to uppercase and ensuring they self$manifest_left_bound_key with 'chr'.

    # Args:
    #     x: A character vector of chromosome labels (e.g., 'chr1', 'CHRX', '2').
    
    # Returns:
    #     A character vector of standardized chromosome labels (e.g., 'chr1', 'chrX', 'chr2').

    return(ifelse(grepl("^chr", x, ignore.case = TRUE), substr(x, 4, nchar(x)), x))
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
        beta = sum(beta / varbeta) / sum(1 / varbeta),
        varbeta = 1 / sum(1 / varbeta),
        pval = 2 * pnorm(-abs(beta / sqrt(varbeta))),
        MAF = suppressWarnings(mean(MAF, na.rm = TRUE)),
        EA = dplyr::first(EA),     # effect/ALT allele defining beta sign
        OA = dplyr::first(OA),
        POS = dplyr::first(POS),
        .groups = "drop"
        ) %>%
        dplyr::filter(is.finite(beta), is.finite(varbeta), varbeta > 0)
    )
}

.read_ld = function(ld_path, bim_path, panel_note) {
    
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
    .log("  Reading LD matrix: %s", basename(ld_path))
    bim <- tryCatch(fread(bim_path, header = FALSE, data.table = FALSE),
                    error = function(e) NULL)

    # Check if the BIM file is missing or empty
    if (is.null(bim) || !nrow(bim)) {
        .log("  BIM missing/empty: %s", as.character(bim_path))
        return(NULL)
    }

    .log("  BIM file loaded: %d variants", nrow(bim))

    # Extract chromosome, position, and alleles from the BIM file and create canonical keys
    # PLINK .bim: V1 chr, V2 id, V3 cM, V4 bp, V5 A1, V6 A2 ; --r counts A1
    chr <- .norm_chr(bim$V1)
    pos <- as.integer(bim$V4)
    a1  <- toupper(bim$V5)
    a2  <- toupper(bim$V6)
    keys <- bim$V2

    # Read the LD matrix and handle any errors
    M <- tryCatch(as.matrix(fread(ld_path, header = FALSE, data.table = FALSE)),
                    error = function(e) {
                        .log("  ERROR reading LD matrix: %s", conditionMessage(e))
                        NULL
                    })
    if (is.null(M) || nrow(M) != length(keys) || ncol(M) != length(keys)) {
        .log("  LD shape mismatch (%s): matrix %dx%d vs bim %d",
            basename(ld_path), nrow(M %||% matrix(0)), ncol(M %||% matrix(0)), length(keys))
        return(NULL)
    }

    .log("  LD matrix loaded: %dx%d", nrow(M), ncol(M))


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
}

.align_ld = function(ld, snps, target_allele) {
    
    # Aligns the LD matrix to a target allele for each variant, subsets it to the specified SNPs, and returns the aligned LD matrix.

    # Args:
    #     ld: A list containing the LD matrix and allele information (output from .read_ld).
    #     snps: A character vector of SNPs to subset the LD matrix to.
    #     target_allele: A named vector of target alleles for each SNP (names should match the SNPs) in which the datasets beta is on (ie. GWAS effect allele / GTEx ALT).
    
    # Returns:
    #     A numeric matrix representing the aligned LD matrix for the specified SNPs, with rows and columns corresponding to the SNPs in the order they appear in the `snps` vector.
    #     Returns NULL if there are fewer than 2 SNPs present or if the LD matrix is NULL.
    
    # Check if the LD matrix is NULL and return NULL if it is
    if (is.null(ld)) return(NULL)

    # Subset the SNPs to those present in the LD matrix
    present <- snps[snps %in% rownames(ld$M)]

    # Log SNP match percentage
    pct_match <- 100 * length(present) / length(snps)
    if (pct_match < 80) {
        .log("    WARNING .align_ld: only %.1f%% SNP match (present=%d, input=%d)", 
            pct_match, length(present), length(snps))
    }

    # If there are fewer than 2 SNPs present, return NULL
    if (length(present) < 2) return(NULL)

    # Extract the LD matrix and allele information for the present SNPs
    M  <- ld$M[present, present, drop = FALSE]
    a1 <- ld$A1[present]
    a2 <- ld$A2[present]
    tgt <- target_allele[present]

    # Align the LD sign to the effect allele
    # If target allele is the same as the LD panels counted allele (a1, ref), flip is needed
    # If target allele is the same as the LD panels non-counted allele (a2, alt), no flip is needed
    s <- rep(NA_real_, length(present))
    s[a1 == tgt] <- -1
    s[a2 == tgt] <- 1

    # Log allele mismatches. If more than 10% of SNPs have mismatched alleles, return NULL to force fallback to coloc.abf
    unmatched <- sum(is.na(s))
    if (unmatched > 0) {
        pct_unmatched <- 100 * unmatched / length(present)
        .log("    WARNING .align_ld: %.1f%% allele mismatch (%d/%d SNPs have target not in A1/A2)", 
            pct_unmatched, unmatched, length(present))
        if (unmatched >= 0.1 * length(present)) {
            .log("      This may cause SuSiE convergence issues or 'unreasonably large prior variance' errors")
            return(NULL)
        }
    }

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
}


# ------------------------------ Class Definition -------------------------------

Finemapper <- R6Class("Finemapper",
    public = list(

        # Fields
        ld_dir = NULL,
        ld_manifest = NULL,
        ld_manifest_path = NULL,
        ss = NULL,
        ss_path = NULL,
        ss_target_allele_key = NULL,
        ss_sample_size = NULL,
        ss_type = NULL,
        ss_s = NULL,
        ss_sdY = NULL,
        drop_ambiguous = NULL,
        min_overlap = NULL,
        susie_min_snps = NULL,
        cred_coverage = NULL,
        susie_max_iter = NULL,
        susie_l = NULL,
        susie_repeat_until_converged = NULL,
        manifest_locus_id_key = NULL,
        manifest_left_bound_key = NULL,
        manifest_right_bound_key = NULL,
        manifest_ld_key = NULL,
        manifest_bim_key = NULL,
        manifest_note_key = NULL,
        standardized_variant_id_key = NULL,
        standardized_chr_key = NULL,
        standardized_pos_key = NULL,
        standardized_effect_allele_key = NULL,
        standardized_non_effect_allele_key = NULL,
        standardized_p_key = NULL,
        standardized_beta_key = NULL,
        standardized_var_beta_key = NULL,
        standardized_se_key = NULL,
        standardized_maf_key = NULL,
        skip_MHC = NULL,
        MHC_chr = NULL,
        MHC_start = NULL,
        MHC_end = NULL,
        susie_results_dir = NULL,
        ss_label = NULL,
        variant_info_dt = NULL,
        result_details_dt = NULL,

        # Initialize the class with paths and column keys
        initialize = function(args) {

            # Inputs/Outputs
            self$ld_dir <- file.path(args$ld_dir)
            self$ld_manifest <- args$ld_manifest
            self$ld_manifest_path <- file.path(self$ld_dir, args$ld_manifest)
            self$ss_path <- file.path(args$ss_fp)
            self$ss <- .load_table(self$ss_path)
            self$ss_target_allele_key <- args$ss_target_allele_key
            self$ss_sample_size <- args$ss_sample_size
            self$ss_type <- args$ss_type
            self$ss_s <- args$ss_case_fraction
            self$ss_sdY <- args$ss_sdY
            self$susie_results_dir <- args$susie_results_dir
            self$ss_label <- args$ss_label
            self$variant_info_dt <- data.table(
                        locus = character(),
                        variant_idx = integer(),
                        variant_id = character(),
                        pip = numeric(),
                        z = numeric(),
                        cs_id = character()
                    )
            self$result_details_dt <- data.table(
                        locus = character(),
                        n_variants = integer(),
                        lead_variant = character(),
                        lead_variant_pval = numeric(),
                        n_variants_ld = integer(),
                        n_credible_sets = integer(),
                        credible_set_lengths = numeric(),
                        lead_in_credible_sets = logical(),
                        susie_status = character()
                    )

            # Parameters
            self$drop_ambiguous <- args$drop_ambiguous
            self$min_overlap <- args$min_overlap
            self$susie_min_snps <- args$susie_min_snps
            self$cred_coverage <- args$cred_coverage
            self$susie_max_iter <- args$susie_max_iter
            self$susie_l <- args$susie_l
            self$susie_repeat_until_converged <- args$susie_repeat_until_converged

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
            self$standardized_effect_allele_key <- args$standardized_effect_allele_key
            self$standardized_non_effect_allele_key <- args$standardized_non_effect_allele_key
            self$standardized_p_key <- args$standardized_p_key
            self$standardized_beta_key <- args$standardized_beta_key
            self$standardized_var_beta_key <- args$standardized_var_beta_key
            self$standardized_se_key <- args$standardized_se_key
            self$standardized_maf_key <- args$standardized_maf_key

            # MHC Region
            self$skip_MHC <- args$skip_MHC
            self$MHC_chr <- 6
            self$MHC_start <- 28000000
            self$MHC_end <- 34000000
        },

        susie_qc = function(D, label) {
            # Performs quality control checks on the SuSiE analysis by checking for non-finite values, calculating Z-scores, estimating s rss, and running kriging rss to detect possible allele flips in the LD matrix.

            # Args:
            #     D: A list representing the dataset for coloc analysis, including beta, variance of beta, SNPs, positions, type, sample size, and aligned LD matrix.
            #     label: A string label for the dataset, used for logging and file naming.

            # Returns:
            #     A list containing:
            #         - s_rss: A numeric value representing the estimated proportion of variance explained by the SNPs (s rss) from the susie_qc function.
            #         - flagged_snps: A character vector of SNPs that were flagged for possible allele flips in the LD matrix
            #     Returns NULL if there are non-finite values in beta or varbeta, or if the susie_rss_estimate_s or susie_rss_kriging_rss functions fail, indicating that the SuSiE analysis should be skipped.

            .log("  Running susie_qc for %s", label)
            # Check for non-finite values. If any non-finite values are found, log a warning and return FALSE to indicate that the SuSiE analysis should be skipped.
            if (!all(is.finite(D$beta)) || !all(is.finite(D$varbeta)) || any(D$varbeta <= 0)) {
                .log("  [%s] Non-finite beta/varbeta detected -> skip SuSiE", label)
                return(list(NULL, NULL))
            }

            # Estimate s rss (the proportion of variance explained by the SNPs) using the susie_rss_estimate_s function, handling any errors
            s_rss <- tryCatch(estimate_s_rss(D$z, D$LD, D$N), error = function(e) NULL)
            if (is.null(s_rss)) {
                .log("  [%s] susie_rss_estimate_s failed -> skip SuSiE", label)
                return(list(NULL, NULL))
            }

            # Run kriging rss to detect possible allele flips in the LD matrix, handling any errors
            krig <- tryCatch(kriging_rss(D$z, D$LD, D$N, 1e-08, s_rss), error = function(e) NULL)
            if (is.null(krig)) {
                .log("  [%s] susie_rss_kriging_rss failed", label)
                flag <- NULL
            } else {

                # Flag SNPs with absolute standardized Z-score difference greater than 3 and log likelihood ratio greater than 2, indicating possible allele flips
                flag <- abs(krig$conditional_dist$z_std_diff) > 3 & krig$conditional_dist$logLR > 2
            }

            # Return s_rss and the flagged SNPs indicating possible allele flips
            return(list(s_rss, flag))
        },

        build_dataset = function(tbl, ld, type, N, s = NULL, sdY = NULL, target_allele_key, label, lead_variant) {
            
            # Builds a dataset for coloc analysis by standardizing column names, collapsing duplicate SNPs, filtering ambiguous SNPs, aligning the LD matrix, and returning a list containing the dataset and the number of SNPs.

            # Args:
            #     tbl: A data frame containing the input data with standardized column names for SNPs, beta, variance of beta, MAF, effect allele, non-effect allele, and position.
            #     ld: A list containing the LD matrix and allele information (output from .read_ld).
            #     type: A string indicating the type of dataset (e.g., 'quant' for quantitative traits).
            #     N: An integer representing the sample size.
            #     s: An optional numeric vector representing the standard errors of the beta estimates.
            #     sdY: An optional numeric value representing the standard deviation of the trait.
            #     target_allele_key: A string indicating the key for the target allele in the input data (e.g., 'EA' for effect allele).
            #     label: A string label for the dataset, used for logging and file naming.
            #     lead_variant: A string representing the lead variant for the dataset, used for logging and file naming.

            # Returns:
            #     A list containing:
            #         - D: A list representing the dataset for coloc analysis, including beta, variance of beta, SNPs, positions, type, sample size, and aligned LD matrix.
            #         - n: An integer representing the number of SNPs in the dataset.
            #         - flagged_snps: A character vector of SNPs that were flagged for possible allele flips in the LD matrix.
            #         - s_rss: A numeric value representing the estimated proportion of variance explained by the SNPs (s rss) from the susie_qc function.
            #         - qc_note: A string note indicating any QC issues detected during the susie_qc function (e.g., "Possible allele flips detected").
            #         - too_few: A logical value indicating whether there were too few SNPs in the dataset for SuSiE analysis (TRUE if fewer than the minimum number of SNPs required, FALSE otherwise).
            #     Returns NULL if there are fewer than 2 SNPs in the dataset or if the aligned LD matrix is NULL.
            
            .log("  Building dataset for %s", label)

            # Initialize qc_note to NA_character_
            qc_note <- NA_character_

            # Mutate the input table to create a new data frame D0 with standardized column names for SNPs, beta, variance of beta, MAF, effect allele, non-effect allele, and position
            D0 <- tbl %>% transmute(
                snp = .data[[self$standardized_variant_id_key]],
                beta = .data[[self$standardized_beta_key]],
                varbeta = .data[[self$standardized_var_beta_key]],
                pval = .data[[self$standardized_p_key]],
                MAF = .data[[self$standardized_maf_key]],
                EA = .data[[self$standardized_effect_allele_key]],
                OA = .data[[self$standardized_non_effect_allele_key]],
                POS = .data[[self$standardized_pos_key]]
            )

            # Collapse duplicate SNPs by keeping the one with the smallest p-value (largest absolute z-score)
            D0 <- .collapse_keys(D0)
            if (self$drop_ambiguous) D0 <- D0 %>% filter(!.is_ambiguous(EA, OA))

            # Remove SNPs with zero or unreliable variance/SE
            D0 <- D0 %>% 
                filter(is.finite(beta), is.finite(varbeta), varbeta > 0, !is.na(beta), beta != 0)
            
            # Check if there are fewer than 2 SNPs in the dataset and return NULL if so
            if (nrow(D0) < 2) return(NULL)

            # Align the LD matrix to the effect allele and subset it to the SNPs in D0
            if (target_allele_key == self$standardized_effect_allele_key) {
                LD <- .align_ld(ld, D0$snp, setNames(D0$EA, D0$snp))
            }
            else if (target_allele_key == self$standardized_non_effect_allele_key) {
                LD <- .align_ld(ld, D0$snp, setNames(D0$OA, D0$snp))
            }
            else {
                .log("  [%s] Unknown target allele key: %s", label, target_allele_key)
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

            # Create a list D containing the necessary information for the coloc analysis, including beta, variance of beta, SNPs, positions, type, sample size, aligned LD matrix, and Z-scores
            D <- list(
                beta = D0$beta,
                varbeta = D0$varbeta,
                pval = D0$pval,
                snp = D0$snp,
                position = D0$POS,
                type = type,
                N = N,
                LD = LD,
                MAF = ifelse(is.finite(D0$MAF), D0$MAF, 0.5),
                z = D0$beta / sqrt(D0$varbeta)
            )

            # Run the susie_qc function to check for non-finite values, estimate s_rss, and detect possible allele flips in the LD matrix
            qc_result <- self$susie_qc(D, label)
            s_rss <- qc_result[[1]]
            flag <- qc_result[[2]]

            # If any possible allele flips are detected, log a message indicating the flagged SNPs
            if (!is.null(flag) && any(flag)) {
                flagged_snps <- D$snp[flag]
                .log("  [%s] Possible allele flips detected for SNPs: %s", label, paste(flagged_snps, collapse = ", "))
            } else {
                flagged_snps <- NULL
            }

            # If s_rss is not NULL and is greater than 0.1, log a warning message indicating that s_rss is large and may indicate mismatched LD or allele flips. If any flagged SNPs are present, log a message indicating that the flagged SNPs will be checked to see if they explain the large s_rss.
            if (!is.null(s_rss)) {
                if (s_rss > 0.1) {
                    .log("  [%s] WARNING: s_rss is large (%.3f) -> may indicate mismatched LD or allele flips", label, s_rss)
                    if (!is.null(flagged_snps)) {

                        # Extract indices for flagged SNPs
                        idx_flagged <- D$snp %in% flagged_snps
                        
                        # Check if flipping the sign of the flagged SNPs resolves the large s_rss. Create a new dataset D_flipped by flipping the sign of the beta and z-scores for the flagged SNPs and run QC on the modified dataset. If s_rss is reduced to below 0.1, log a message indicating that flipping the flagged SNPs explains the large s_rss and update D and s_rss accordingly. Otherwise, log a message indicating that flipping the flagged SNPs does not explain the large s_rss and keep the flagged SNPs in the analysis.
                        .log("  [%s] Checking if allele flips or strand mismatches explain large s_rss", label)
                        
                        # Create a new dataset D_flipped by flipping the sign of the beta and z-scores for the flagged SNPs
                        D_flipped <- list(
                            beta = D$beta,
                            varbeta = D$varbeta,
                            pval = D$pval,
                            snp = D$snp,
                            position = D$position,
                            type = D$type,
                            N = D$N,
                            LD = D$LD,
                            MAF = D$MAF,
                            z = D$z
                        )
                        D_flipped$beta[idx_flagged] <- -D_flipped$beta[idx_flagged]
                        D_flipped$z[idx_flagged] <- -D_flipped$z[idx_flagged]

                        # Run QC on the modified dataset
                        label_flipped <- paste0(label, "_flipped")
                        qc_result_flipped <- self$susie_qc(D_flipped, label_flipped)
                        s_rss_flipped <- qc_result_flipped[[1]]
                        if (!is.null(s_rss_flipped) && s_rss_flipped < 0.1) {
                            qc_note <- "Allele flips or strand mismatches explain large s_rss -> flipping sign of flagged SNPs in analysis"
                            .log("  %s", qc_note)
                            D <- D_flipped
                            s_rss <- s_rss_flipped
                        } else {
                            .log("  [%s] Allele flips or strand mismatches do not explain large s_rss -> check removing flagged SNPs", label)

                            # Check if removing flagged SNPs resolves the large s_rss. If the lead variant is not among the flagged SNPs, create a new dataset D_no_flagged by removing the flagged SNPs and run QC on the modified dataset. If s_rss is reduced to below 0.1, log a message indicating that the flagged SNPs explain the large s_rss and update D and s_rss accordingly. Otherwise, log a message indicating that the flagged SNPs do not explain the large s_rss and keep the flagged SNPs in the analysis.
                            if (!lead_variant %in% flagged_snps){
                                .log("  [%s] Checking if removing flagged SNPs resolves large s_rss", label)

                                # Create D_no_flagged by removing flagged SNPs from all list components
                                D_no_flagged <- list(
                                    beta = D$beta[!idx_flagged],
                                    varbeta = D$varbeta[!idx_flagged],
                                    pval = D$pval[!idx_flagged],
                                    snp = D$snp[!idx_flagged],
                                    position = D$position[!idx_flagged],
                                    type = D$type,
                                    N = D$N,
                                    LD = D$LD[!idx_flagged, !idx_flagged, drop = FALSE],
                                    MAF = D$MAF[!idx_flagged],
                                    z = D$z[!idx_flagged]
                                )

                                # Run QC on the modified dataset
                                label_no_flagged <- paste0(label, "_no_flagged")
                                qc_result <- self$susie_qc(D_no_flagged, label_no_flagged)
                                s_rss_no_flagged <- qc_result[[1]]
                                if (!is.null(s_rss_no_flagged) && s_rss_no_flagged < 0.1) {
                                    qc_note <- "Flagged SNPs explain large s_rss -> removing flagged SNPs from analysis"
                                    .log("  %s", qc_note)
                                    D <- D_no_flagged
                                    s_rss <- s_rss_no_flagged
                                } else {
                                    qc_note <- "Flagged SNPs do not explain large s_rss -> keeping flagged SNPs in analysis but results may still be unreliable"
                                    .log("  %s", qc_note)
                                }
                            } else {
                                qc_note <- "Lead variant is flagged -> keeping flagged SNPs in analysis but results may still be unreliable"
                                .log("  %s", qc_note)
                            }
                        }
                    } else {
                        qc_note <- "s_rss is large but no flagged SNPs -> results may be unreliable"
                        .log("  %s", qc_note)
                    }
                } else {
                    qc_note <- "s_rss is small -> no issues detected"
                    .log("  %s", qc_note)
                }
            } else {
                qc_note <- "s_rss could not be estimated -> results may be unreliable"
                .log("  %s", qc_note)
            }

            # Add optional parameters s and sdY to the list D if they are not NULL
            if (D$type == "cc" && !is.null(s)) D$s <- s
            if (D$type == "quant" && !is.null(sdY)) D$sdY <- sdY
            
            # Create a list containing the dataset D and the number of SNPs in D0
            dataset <- list(D = D, n = length(D$snp), flagged_snps = paste0(flagged_snps, collapse = ", ") %||% NA_character_, s_rss = s_rss %||% NA_real_, qc_note = qc_note %||% NA_character_, too_few = FALSE)

            # Return the list containing the dataset and the number of SNPs
            return(dataset)
        },


        safe_runsusie = function(D, label) {
            
            # Safely runs the SuSiE algorithm on a given dataset, handling errors and logging the results.

            # Args:
            #     D: A list representing the dataset for SuSiE analysis, including beta, variance of beta, SNPs, positions, type, sample size, and aligned LD matrix.
            #     label: A string label for logging purposes.

            # Returns:
            #     A list containing:
            #         - fit: The result of the SuSiE analysis if successful, or NULL if there was an error.
            #         - reason: A string indicating the reason for failure, or NULL if the analysis was successful.
            
            # Check if the dataset is valid for SuSiE analysis using the check_dataset function. If it fails, log the error and return "dataset_check_failed".
            dataset_check <- tryCatch(check_dataset(D, req = "LD"), error = function(e) conditionMessage(e))
            if (!is.null(dataset_check)) {     # NULL means OK
                ran_susie_reason <- paste0("check_dataset failed")
                .log("  [%s] %s with the following error: %s -> skip SuSiE", label, ran_susie_reason, dataset_check)
                return(list(fit = NULL, reason = ran_susie_reason))
            }

            # Attempt to run the SuSiE algorithm on the dataset, suppressing warnings and handling errors
            S <- tryCatch(
                withTimeout({
                    fit <- suppressWarnings(runsusie(
                        D,
                        coverage = self$cred_coverage,
                        max_iter = self$susie_max_iter,
                        L = self$susie_l,
                        repeat_until_convergence = self$susie_repeat_until_converged,
                        compute_univariate_zscore = TRUE
                    ))
                    list(fit = fit, reason = NULL)
                }, timeout = 600, cpu = Inf, elapsed = 600),
                error = function(e) { 
                    if (inherits(e, "TimeoutException")) {
                        ran_susie_reason <- "SuSiE timeout"
                        .log("  [%s] %s -> fall back to SuSiE w/ refined parameters", label, ran_susie_reason)
                    }
                    else {
                        ran_susie_reason <- paste0("SuSiE error: ", gsub("[\r\n\t ]+", " ", conditionMessage(e)))
                        .log("  [%s] %s -> fall back to SuSiE w/ refined parameters", label, ran_susie_reason)
                    }
                    list(fit = NULL, reason = ran_susie_reason)
                }
            )

            # Check if the SuSiE result is NULL and return NULL if it is
            if (is.null(S$fit)) {
                return(S)
            }

            # Count the number of credible sets in the SuSiE result, handling any errors
            ncs <- tryCatch(length(S$fit$sets$cs), error = function(e) 0L)
            .log("  [%s] SuSiE credible sets: %d", label, ncs %||% 0L)

            # Return NULL if there are no credible sets, indicating that there is nothing to colocalize
            if (ncs < 1) {
                ran_susie_reason <- "SuSiE found no credible sets"
                .log("  [%s] %s -> fall back to SuSiE w/ refined parameters", label, ran_susie_reason)
                return(list(fit = NULL, reason = ran_susie_reason))
            }

            # Return the SuSiE result
            return(S)
        },

        run_susie_fallback = function(D, label) {

            # Runs the SuSiE fallback w/ refined parameters on a given dataset. Specifically, explicitly estimate prior variance and reduce max credible sets.

            # Args:
            #     D: A list representing the dataset for SuSiE analysis, including beta, variance of beta, SNPs, positions, type, sample size, and aligned LD matrix.
            #     label: A string label for logging purposes.

            # Returns:
            #     The result of the SuSiE analysis if successful, or NULL if there was an error

            # Attempt to run the SuSiE w/ refined parameters on the dataset, suppressing warnings and handling errors
            S <- tryCatch(
                withTimeout({
                    fit <- suppressWarnings(runsusie(
                        D,
                        coverage = self$cred_coverage, 
                        max_iter = self$susie_max_iter, 
                        L = 5,  # Limit the number of components to 5 for stability
                        estimate_prior_variance = TRUE,  # Enable prior variance estimation for regularization
                        estimate_prior_method = "EM",     # Use EM algorithm for stable estimation
                        scaled_prior_variance = 0.01,     # Regularization strength (smaller = more regularization)
                        repeat_until_convergence = self$susie_repeat_until_converged,
                        compute_univariate_zscore = TRUE
                    ))
                    list(fit = fit, reason = NULL)
                }, timeout = 600, cpu = Inf, elapsed = 600),
                error = function(e) { 
                    if (inherits(e, "TimeoutException")) {
                        ran_susie_reason <- "SuSiE fallback timeout"
                        .log("  [%s] %s -> skip susie.coloc", label, ran_susie_reason)
                    } else {
                        ran_susie_reason <- paste0("SuSiE fallback error: ", gsub("[\r\n\t ]+", " ", conditionMessage(e)))
                        .log("  [%s] %s -> skip susie.coloc", label, ran_susie_reason)
                    }
                    list(fit = NULL, reason = ran_susie_reason)
                }
            )

            # Check if the SuSiE result is NULL and return NULL if it is
            if (is.null(S$fit)) {
                return(S)
            }

            # Count the number of credible sets in the SuSiE result, handling any errors
            ncs <- tryCatch(length(S$fit$sets$cs), error = function(e) 0L)
            .log("  [%s] SuSiE fallback credible sets: %d", label, ncs %||% 0L)

            # Return NULL if there are no credible sets, indicating that there is nothing to colocalize
            if (ncs < 1) {
                ran_susie_reason <- "SuSiE fallback found no credible sets"
                .log("  [%s] %s -> skip susie.coloc", label, ran_susie_reason)
                return(list(fit = NULL, reason = ran_susie_reason))
            }

            # Return the SuSiE result
            return(S)

        },

        process_locus = function(locus) {

            # Run SuSiE finemapping for a single locus

            # Args:
            #     locus: A list containing locus information (e.g., LOCUS_ID, self$standardized_chr_key, self$manifest_left_bound_key, self$manifest_right_bound_key).

            # Returns:
            #     None. The function updates the internal state of the object with SuSiE results and logs relevant information.

            locus_label <- paste0(self$ss_label, "_", locus[[self$manifest_locus_id_key]])

            .log("Processing locus: %s", locus_label)

            # Initialize file paths for SuSiE results and status
            susie_fit_fp <- file.path(self$susie_results_dir, paste0(locus_label, "_susie_fit.rds"))

            # Initialize variables
            susie_fit <- NULL
            susie_status <- NA_character_
            n_signals_susie <- 0L
            n_credible_sets <- 0L
            credible_set_lengths <- integer(0)
            lead_in_credible_sets <- FALSE
            
            # Filter the summary stats to include only variants within the specified locus boundaries (chromosome and position range)
            locus_ss <- self$ss %>% filter(
                .data[[self$standardized_chr_key]] == locus[[self$standardized_chr_key]],
                .data[[self$standardized_pos_key]] >= locus[[self$manifest_left_bound_key]],
                .data[[self$standardized_pos_key]] <= locus[[self$manifest_right_bound_key]]
            )

            # Count the number of variants in the locus and log the information
            n_variants = nrow(locus_ss)
            .log("  [%s] Variants in locus: %d", locus_label, n_variants)

            # Record the lead variant
            top_variant <- locus_ss %>% arrange(.data[[self$standardized_p_key]]) %>% slice(1)
            lead_variant <- top_variant[[self$standardized_variant_id_key]] %||% NA_character_
            lead_variant_pval <- top_variant[[self$standardized_p_key]] %||% NA_real_

            # Check if locus is in MHC region and skip SuSiE analysis if specified in the configuration. Log a message indicating that the locus is being skipped due to being in the MHC region.
            if (self$skip_MHC & 
            locus[[self$standardized_chr_key]] == self$MHC_chr &
            ((locus[[self$manifest_left_bound_key]] >= self$MHC_start & locus[[self$manifest_left_bound_key]] <= self$MHC_end) |
            (locus[[self$manifest_right_bound_key]] >= self$MHC_start & locus[[self$manifest_right_bound_key]] <= self$MHC_end))) {
                .log("  [%s] Skipping SuSiE for locus in MHC region", locus_label)
                ld <- NULL
            } else{
                # Load the LD matrix for the locus using the .read_ld function
                ld <- .read_ld(
                    ld_path  = file.path(self$ld_dir, locus[[self$manifest_ld_key]]),
                    bim_path = file.path(self$ld_dir, locus[[self$manifest_bim_key]]),
                    panel_note = locus[[self$manifest_note_key]])
            }

            # Build summary stats dataset for SuSiE analysis
            dataset <- self$build_dataset(
                tbl = locus_ss,
                ld = ld, type = self$ss_type,
                N = self$ss_sample_size,
                s = self$ss_s,
                sdY = self$ss_sdY,
                target_allele_key = self$ss_target_allele_key,
                label = locus_label,
                lead_variant = lead_variant
            )

            # Log the number of variants in the dataset (= number of SNPs in the LD matrix)
            n_variants_ld <- dataset$n %||% 0L
            .log("  [%s] Variants in dataset (LD): %d", locus_label, n_variants_ld)

            # Skip SuSiE analysis if the dataset has too few variants, logging a message indicating that SuSiE is being skipped due to insufficient variants.
            if (!is.null(dataset) && isFALSE(dataset$too_few)) {
                .log("  [%s] Passed initial checks for dataset", locus_label)

                if (!is.null(dataset$D)) {
                    .log("[%s] Running SuSiE...", locus_label)
                    susie_result <- self$safe_runsusie(
                        D = dataset$D,
                        label = locus_label
                    )

                    # Extract the SuSiE fit and status from the result
                    susie_fit <- susie_result$fit
                    susie_status <- susie_result$reason

                    # Only attempt fallback if we didn't get a "check_dataset failed" sentinel (indicating LD check failure)
                    if (is.null(susie_fit) && !identical(susie_status, "check_dataset failed")) {
                        # Attempt fallback only if safe_runsusie returned NULL (not "check_dataset failed")
                        susie_fallback_result <- self$run_susie_fallback(
                            D = dataset$D,
                            label = locus_label
                        )
                        # Extract the SuSiE fit and status from the fallback result
                        susie_fit <- susie_fallback_result$fit
                        susie_status <- susie_fallback_result$reason
                    }
                }

                # Save SuSiE fit and set ran_susie_gwas flag to TRUE
                if (!is.null(susie_fit)) {

                    # Plot -log10 P-values (Manhattan plot)
                    tryCatch({
                        png(file.path(self$susie_results_dir, paste0(locus_label, "_manhattan_plot.png")), width = 1200, height = 600, res = 150)
                        plot(
                            x = seq_along(dataset$D$pval),
                            y = -log10(dataset$D$pval),
                            xlab = "variable",
                            ylab = expression(-log[10](P)),
                            pch = 16, cex = 0.5
                        )
                    }, error = function(e) .log("  [%s] Manhattan plot failed: %s", locus_label, conditionMessage(e)),
                    finally = if (dev.cur() > 1) dev.off())

                    # Plot PIP Scores
                    tryCatch({
                        png(file.path(self$susie_results_dir, paste0(locus_label, "_PIP_plot.png")), width = 1200, height = 600, res = 150)
                        susie_plot(susie_fit, y = "PIP", add_legend = TRUE)
                    }, error = function(e) .log("  [%s] PIP plot failed: %s", locus_label, conditionMessage(e)),
                    finally = if (dev.cur() > 1) dev.off())

                    # Extract credible sets
                    credible_sets <- susie_fit$sets$cs
                    n_credible_sets <- length(credible_sets)

                    # Determine credible set lengths
                    credible_set_lengths <- sapply(credible_sets, length)

                    # Create a data.table to store locus variants information, including locus name, variant index, variant ID, PIP scores, Z-scores, and credible set IDs (if available)
                    locus_variants_dt <- data.table(
                        locus = locus_label,
                        variant_idx = seq_along(susie_fit$pip),
                        variant_id = names(susie_fit$pip),
                        pip = as.numeric(susie_fit$pip),
                        z = if (is.null(susie_fit$z)) NA_real_ else as.numeric(susie_fit$z),
                        cs_id = NA_character_
                    )

                    # Fill in the credible set IDs for each variant based on the credible sets extracted from the SuSiE fit
                    for (cs_name in names(credible_sets)) {
                        locus_variants_dt[variant_idx %in% credible_sets[[cs_name]], cs_id := cs_name]
                        if (lead_variant %in% names(credible_sets[[cs_name]])) {
                            lead_in_credible_sets <- TRUE
                        }
                    }

                    # Append the locus variants information to the overall susie_variants_dt data.table
                    self$variant_info_dt <- rbindlist(list(self$variant_info_dt, locus_variants_dt))

                    # Save the SuSiE fit for the dataset to an RDS file and log the action
                    .log("[%s] Saving SuSiE fit to %s", locus_label, susie_fit_fp)
                    saveRDS(list(fit = susie_fit, reason = susie_status), susie_fit_fp)
                    susie_status <- "Succeeded and saved"

                }
            } else {
                .log("[%s] Skipping SuSiE for locus due to too few variants: n_variants=%d < min_variants=%d", locus_label, n_variants_ld, self$susie_min_snps)
                susie_status <- sprintf("Skipped due to too few variants: n_variants=%d < min_variants=%d", n_variants_ld, self$susie_min_snps)
            }

            # Append this locus' details as a new row in the result_details_dt data.table
            self$result_details_dt <- rbindlist(list(
                self$result_details_dt,
                data.table(
                    locus = locus_label,
                    n_variants = n_variants,
                    lead_variant = lead_variant,
                    lead_variant_pval = lead_variant_pval,
                    n_variants_ld = n_variants_ld,
                    n_credible_sets = n_credible_sets,
                    credible_set_lengths = paste(credible_set_lengths, collapse = ","),
                    lead_in_credible_sets = lead_in_credible_sets,
                    susie_status = susie_status %||% NA_character_
                )
            ))

            # Clean up memory by removing large objects and running garbage collection
            rm(list = intersect(c("ld", "dataset", "susie_fit", "locus_ss"), ls()))
            gc(verbose = FALSE)
        },

        run_finemapping = function() {

            # Runs fine-mapping analysis for all loci in provided SS datasets using SuSiE, logging results.
            .log("Loading LD Manifest...")
            loci_data <- .load_table(self$ld_manifest_path)
            .log("LD Manifest loaded: %d rows", nrow(loci_data))

            # Define file paths for variant information and result details output files
            variant_info_fp <- file.path(self$susie_results_dir, sprintf("%s_susie_variant_info.tsv", self$ss_label))
            result_details_fp <- file.path(self$susie_results_dir, sprintf("%s_susie_result_details.tsv", self$ss_label))

            # Loop through each locus and process it
            for (i in seq_len(nrow(loci_data))) {
                locus <- loci_data[i, ]
                self$process_locus(locus)
            }

            # Save the variant information and result details data.tables to TSV files
            .log("Saving variant information to %s", variant_info_fp)
            fwrite(self$variant_info_dt, variant_info_fp, sep = "\t", na = "NA", quote = FALSE)
            .log("Saving result details to %s", result_details_fp)
            fwrite(self$result_details_dt, result_details_fp, sep = "\t", na = "NA", quote = FALSE) 
        }
    )
) 

# ------------------------------------- Command Line Interface -------------------------------------

if (!interactive()) {
    
    parser <- ArgumentParser(description = "Colocalization analysis using SuSiE and coloc.abf")
    parser$add_argument("--ld_dir", required = TRUE, help = "Path to the directory containing LD matrices and BIM files.")
    parser$add_argument("--ld_manifest", default = "ld_manifest.tsv", help = "File name of the ld manifest file containing locus information (optional). Default is 'ld_manifest.tsv'.")
    parser$add_argument("--ss_fp", required = TRUE, help = "Path to the standardized summary statistics file")
    parser$add_argument("--ss_target_allele_key", default = "EFFECT", help = "Column name for the target allele in the summary statistics file (optional). Default is 'EFFECT'.")
    parser$add_argument("--ss_sample_size", type = "numeric", required = TRUE, help = "Sample size for the dataset")
    parser$add_argument("--ss_type", choices = c("cc", "quant"), required = TRUE, help = "Type of the trait in the dataset: 'cc' for case-control or 'quant' for quantitative.")
    parser$add_argument("--ss_case_fraction", type = "numeric", required = TRUE, help = "Case fraction for the dataset.")
    parser$add_argument("--ss_sdY", type = "numeric", default = 1, help = "Standard deviation of the trait for the dataset (optional). Default is 1.")
    parser$add_argument("--drop_ambiguous", default = TRUE, type = "logical", help = "Drop ambiguous SNPs (A/T or C/G) from the analysis (optional). Default is TRUE.")
    parser$add_argument("--min_overlap", type = "numeric", default = 50, help = "Minimum number of overlapping SNPs required for colocalization analysis (optional). Default is 50.")
    parser$add_argument("--cred_coverage", type = "numeric", default = 0.95, help = "Credible set coverage for SuSiE analysis (optional). Default is 0.95.")
    parser$add_argument("--susie_min_snps", type = "numeric", default = 5, help = "Minimum number of SNPs required for SuSiE analysis (optional). Default is 5.")
    parser$add_argument("--susie_max_iter", type = "numeric", default = 1000, help = "Maximum number of iterations for SuSiE analysis (optional). Default is 1000.")
    parser$add_argument("--susie_l", type = "numeric", default = 10, help = "Maximum number of causal variants for SuSiE analysis (optional). Default is 10.")
    parser$add_argument("--susie_repeat_until_converged", default = FALSE, type = "logical", help = "Repeat SuSiE until convergence (optional). Default is FALSE.")
    parser$add_argument("--manifest_locus_id_key", default = "LOCUS_ID", help = "Column name for the locus ID in the LD manifest file (optional). Default is 'LOCUS_ID'.")
    parser$add_argument("--manifest_left_bound_key", default = "LEFT_500KB", help = "Column name for the left bound of the locus in the LD manifest file (optional). Default is 'LEFT_500KB'.")
    parser$add_argument("--manifest_right_bound_key", default = "RIGHT_500KB", help = "Column name for the right bound of the locus in the LD manifest file (optional). Default is 'RIGHT_500KB'.")
    parser$add_argument("--manifest_ld_key", default = "LD", help = "Column name for the LD matrix file path in the LD manifest file (optional). Default is 'LD'.")
    parser$add_argument("--manifest_bim_key", default = "BIM", help = "Column name for the BIM file path in the LD manifest file (optional). Default is 'BIM'.")
    parser$add_argument("--manifest_note_key", default = "note", help = "Column name for the LD panel note in the LD manifest file (optional). Default is 'PANEL_NOTE'.")
    parser$add_argument("--standardized_variant_id_key", default = "VAR", help = "Column name for the standardized variant ID in the summary statistics files (optional). Default is 'VAR'.")
    parser$add_argument("--standardized_chr_key", default = "CHR", help = "Column name for the standardized chromosome in the summary statistics files (optional). Default is 'CHR'.")
    parser$add_argument("--standardized_pos_key", default = "BP", help = "Column name for the standardized position in the summary statistics files (optional). Default is 'POS'.")
    parser$add_argument("--standardized_non_effect_allele_key", default = "NON_EFFECT", help = "Column name for the standardized reference allele in the summary statistics files (optional). Default is 'NON_EFFECT'.")
    parser$add_argument("--standardized_effect_allele_key", default = "EFFECT", help = "Column name for the standardized effect allele in the summary statistics files (optional). Default is 'EFFECT'.")
    parser$add_argument("--standardized_p_key", default = "P", help = "Column name for the standardized p-value in the summary statistics files (optional). Default is 'P'.")
    parser$add_argument("--standardized_beta_key", default = "BETA", help = "Column name for the standardized beta in the summary statistics files (optional). Default is 'BETA'.")
    parser$add_argument("--standardized_var_beta_key", default = "VARBETA", help = "Column name for the standardized variance of beta in the summary statistics files (optional). Default is 'VAR_BETA'.")
    parser$add_argument("--standardized_se_key", default = "SE", help = "Column name for the standardized standard error in the summary statistics files (optional). Default is 'SE'.")
    parser$add_argument("--standardized_maf_key", default = "MAF", help = "Column name for the standardized minor allele frequency in the summary statistics files (optional). Default is 'MAF'.")
    parser$add_argument("--skip_MHC", default = TRUE, type = "logical", help = "Skip loci in the MHC region (optional). Default is TRUE.")
    parser$add_argument("--susie_results_dir", required = TRUE, help = "Directory to save SuSiE results (optional). If not provided, SuSiE results will be saved in the output directory.")
    parser$add_argument("--ss_label", default = "", help = "Label for the summary statistics (optional). If not provided, the default label will be used.")

    # Parse the command line arguments
    args <- parser$parse_args()

    # Print all arguments (including defaults)
    .log("=== Command Line Arguments ===")
    for (arg_name in names(args)) {
        arg_value <- args[[arg_name]]
        if (is.character(arg_value)) {
            .log("  %s = '%s'", arg_name, arg_value)
        } else if (is.logical(arg_value)) {
            .log("  %s = %s", arg_name, arg_value)
        } else if (is.numeric(arg_value)) {
            .log("  %s = %g", arg_name, arg_value)
        } else {
            .log("  %s = %s", arg_name, paste(arg_value, collapse = ", "))
        }
    }
    .log("===============================")

    # Create an instance of the Finemapper class
    finemapper <- Finemapper$new(args)

    # Run the colocalization analysis
    finemapper$run_finemapping()
}