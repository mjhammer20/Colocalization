#!/usr/bin/env Rscript

# Colocalization analysis (SuSiE coloc + coloc.abf) that inherits all fine-mapping logic
# (LD loading/alignment, dataset building, QC, cascading rescue, SuSiE + fallback) from
# finemap_susie_cascading_rescue.R. ColocalizationAnalyzer is a subclass of Finemapper.

# ------------------------------ Source parent script -------------------------------

.this_dir <- function() {

    # Returns the directory containing this script (works for Rscript and source()).

    f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
    if (length(f)) return(dirname(normalizePath(f[1])))
    of <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
    if (!is.null(of)) return(dirname(normalizePath(of)))
    return(".")
}

# Provides the libraries, the LD/dataset helper functions (.load_table, .log, .norm_chr, .is_ambiguous, .collapse_keys, .read_ld, .reconcile_ld_ids, .align_ld) and the Finemapper class
source(file.path(.this_dir(), "finemap_susie_cascading_rescue.R"))

# ------------------------------ Helper Function Definitions -------------------------------

.stream_write <- function(df, path) {

    # Writes a data frame to a file in a streaming manner. If the file already exists, it appends the data frame to the existing file without writing column names.

    # Args:
    #     df: The data frame to be written.
    #     path: The file path to write the data frame to.

    # Returns:
    #     None. The function writes the data frame to the specified file path.

    if (!file.exists(path)) data.table::fwrite(df, path, sep = "\t", quote = FALSE)
    else data.table::fwrite(df, path, sep = "\t", quote = FALSE, append = TRUE, col.names = FALSE)
}


.collapse_flagged <- function(dataset) {

    # Collapses the flagged SNP vector of a dataset into a single comma separated string (NA if none).

    flagged <- dataset$flagged_snps
    if (is.null(flagged) || length(flagged) == 0) return(NA_character_)
    return(paste(flagged, collapse = ", "))
}


# ------------------------------ Class Definition -------------------------------

ColocalizationAnalyzer <- R6Class("ColocalizationAnalyzer",
    inherit = Finemapper,
    public = list(

        # Fields (everything else is inherited from Finemapper)
        output_dir = NULL,
        qc_dir = NULL,
        gwas_path = NULL,
        gwas_target_allele_key = NULL,
        gwas_strata_key = NULL,
        gwas_sample_size = NULL,
        gwas_s = NULL,
        qtl_path = NULL,
        qtl_target_allele_key = NULL,
        qtl_strata_key = NULL,
        qtl_sample_size = NULL,
        qtl_sdY = NULL,
        pp_h4_threshold = NULL,
        max_s_rss = NULL,
        coloc_priors = NULL,
        standardized_gene_id_key = NULL,
        bucket = NULL,
        gwas_label = NULL,
        qtl_label = NULL,
        gwas_susie_results_dir = NULL,
        qtl_susie_results_dir = NULL,
        gwas_variant_info_dt = NULL,
        gwas_result_details_dt = NULL,
        qtl_variant_info_dt = NULL,
        qtl_result_details_dt = NULL,

        # Initialize the class with paths and column keys
        initialize = function(args) {

            # Initialize the inherited Finemapper state. The GWAS acts as the "ss" for the parent; the GWAS/QTL tables are loaded in run_coloc_analysis, so no ss_fp is passed.
            parent_args <- args
            parent_args$ss_fp <- NULL
            parent_args$ss_target_allele_key <- args$gwas_target_allele_key
            parent_args$ss_sample_size <- args$gwas_sample_size
            parent_args$ss_type <- "cc"
            parent_args$ss_case_fraction <- args$gwas_case_fraction
            parent_args$ss_sdY <- NULL
            parent_args$ss_label <- args$gwas_label
            parent_args$susie_results_dir <- NULL   # GWAS and QTL use separate results directories (set below)
            super$initialize(parent_args)

            # Inputs/Outputs
            self$output_dir <- args$output_dir

            # Separate SuSiE results directories and result tables for GWAS and QTL (copies of the empty tables defined by the parent)
            self$gwas_susie_results_dir <- args$gwas_susie_results_dir %||% file.path(args$output_dir, "gwas_susie_results")
            self$qtl_susie_results_dir <- args$qtl_susie_results_dir %||% file.path(args$output_dir, "qtl_susie_results")
            dir.create(self$gwas_susie_results_dir, recursive = TRUE, showWarnings = FALSE)
            dir.create(self$qtl_susie_results_dir, recursive = TRUE, showWarnings = FALSE)
            self$gwas_variant_info_dt <- copy(self$variant_info_dt)
            self$gwas_result_details_dt <- copy(self$result_details_dt)
            self$qtl_variant_info_dt <- copy(self$variant_info_dt)
            self$qtl_result_details_dt <- copy(self$result_details_dt)
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
            self$gwas_label <- args$gwas_label
            self$qtl_label <- args$qtl_label

            # Colocalization parameters
            self$pp_h4_threshold <- args$pp_h4_threshold
            self$max_s_rss <- args$max_s_rss %||% 0.3
            self$coloc_priors <- list(
                p1 = args$coloc_prior_p1,
                p2 = args$coloc_prior_p2,
                p12 = args$coloc_prior_p12
            )
            self$standardized_gene_id_key <- args$standardized_gene_id_key

            # Tracking
            self$bucket <- c(susie_pairs = 0L, coloc_abf_pairs = 0L, low_overlap = 0L)
        },


        susie_dir_for = function(source) {

            # Returns the SuSiE results directory for the given source ("gwas" or "qtl").

            if (identical(source, "gwas")) return(self$gwas_susie_results_dir)
            if (identical(source, "qtl")) return(self$qtl_susie_results_dir)
            stop(sprintf("Unknown source '%s' (expected 'gwas' or 'qtl')", source))
        },


        run_finemapping_for_coloc = function(tbl, ld, type, N, s = NULL, sdY = NULL, target_allele_key, label, locus_in_mhc, source) {

            # Builds a dataset from a summary stats table and fine-maps it with the inherited cascading rescue, for use in colocalization.
            # Reuses a precomputed SuSiE fit if one exists. A fit with no credible sets is treated as a failure since there is nothing to colocalize.

            # Args:
            #     tbl: Summary stats for the locus (GWAS table or single gene QTL table).
            #     ld: LD list from load_locus_ld (NULL if unavailable or skipped).
            #     type: "cc" or "quant".
            #     N, s, sdY: Sample size, case fraction, and trait sd passed to build_dataset.
            #     target_allele_key: Standardized allele key the beta is reported on.
            #     label: Label used for logging.
            #     locus_in_mhc: TRUE if the locus is skipped because it is in the MHC region.
            #     source: "gwas" or "qtl". Selects the results directory used to cache the SuSiE fit (<label>_susie_fit.rds, written by record_locus_results).

            # Returns:
            #     A list containing:
            #         - dataset: The dataset from build_dataset (the rescued dataset if rescue succeeded; may be NULL or a too_few list).
            #         - n_ld: Number of variants in the dataset after LD alignment.
            #         - fit: The SuSiE fit usable for colocalization (NULL if unavailable).
            #         - n_signals: Number of credible sets in the fit.
            #         - ran: TRUE if a usable SuSiE fit is available.
            #         - status: Status string describing what happened.
            #         - attempted: TRUE if SuSiE was run (or a precomputed fit was loaded) for this dataset.
            #         - susie_result: List with fit, reason, and mode for record_locus_results.
            #         - n_variants, n_rare_variants, lead_variant, lead_variant_pval: Locus summary values for record_locus_results.

            # Locus summary values (needed by record_locus_results)
            fit_fp <- file.path(self$susie_dir_for(source), paste0(label, "_susie_fit.rds"))
            n_variants <- nrow(tbl)
            n_rare_variants <- tbl %>% filter(.data[[self$standardized_eaf_key]] < 0.01 | .data[[self$standardized_eaf_key]] > 0.99) %>% nrow()
            top_variant <- tbl %>% filter(is.finite(.data[[self$standardized_p_key]])) %>% arrange(.data[[self$standardized_p_key]]) %>% slice(1)
            lead_variant <- if (nrow(top_variant) > 0) top_variant[[self$standardized_variant_id_key]][1] else NA_character_
            lead_variant_pval <- if (nrow(top_variant) > 0) top_variant[[self$standardized_p_key]][1] else NA_real_

            # Build dataset (inherited)
            self$skip_reason <- if (locus_in_mhc) "MHC region" else NULL
            dataset <- self$build_dataset(
                tbl = tbl, ld = ld, type = type, N = N, s = s, sdY = sdY,
                target_allele_key = target_allele_key, label = label
            )

            n_ld <- dataset$n %||% 0L
            .log("  [%s] Variants in dataset (LD): %d", label, n_ld)

            out <- list(
                dataset = dataset, n_ld = n_ld, fit = NULL, n_signals = 0L, ran = FALSE, status = NA_character_,
                attempted = FALSE, susie_result = list(fit = NULL, reason = NULL, mode = NULL),
                n_variants = n_variants, n_rare_variants = n_rare_variants,
                lead_variant = lead_variant, lead_variant_pval = lead_variant_pval
            )

            # Skip if the dataset could not be built or has too few variants
            if (is.null(dataset) || !isFALSE(dataset$too_few)) {
                reason <- self$skip_reason %||% sprintf("too few variants: n_variants=%d < min_variants=%d", n_ld, self$susie_min_snps)
                .log("[%s] Skipping SuSiE: %s", label, reason)
                out$status <- sprintf("Skipped due to %s", reason)
                return(out)
            }

            .log("  [%s] Passed initial checks for dataset", label)

            # Load a precomputed fit if available
            if (file.exists(fit_fp)) {
                .log("[%s] Loading precomputed SuSiE fit from %s", label, fit_fp)
                cached <- readRDS(fit_fp)
                cached_fit <- cached$fit
                n_signals <- tryCatch(length(cached_fit$sets$cs), error = function(e) 0L) %||% 0L
                .log("  [%s] Precomputed SuSiE fit credible sets: %d", label, n_signals)

                # Restore the rescue information so the results tables match the original run
                out$dataset$rescue_mode <- cached$rescue_mode
                out$dataset$s_rss_rescued <- cached$s_rss_rescued
                out$susie_result <- list(fit = cached_fit, reason = NULL, mode = cached$mode)
                out$attempted <- TRUE
                out$n_signals <- n_signals

                if (is.null(cached_fit) || n_signals < 1) {
                    out$status <- "SuSiE found no credible sets (precomputed)"
                } else {
                    out$fit <- cached_fit
                    out$ran <- TRUE
                    out$status <- "Loaded precomputed"
                }
                return(out)
            }

            # Skip datasets with unreliable LD that cannot be rescued (no flagged SNPs)
            if (is.na(dataset$s_rss) || (dataset$s_rss >= self$max_s_rss && is.null(dataset$flagged_snps))) {
                .log("[%s] Skipping SuSiE due to unreliable LD: estimated s_rss=%.3f >= %.2f and no flagged SNPs to rescue", label, dataset$s_rss, self$max_s_rss)
                out$status <- sprintf("Skipped due to unreliable LD: estimated s_rss=%.3f >= %.2f", dataset$s_rss, self$max_s_rss)
                return(out)
            }

            # Run fine-mapping with cascading rescue (inherited)
            fm <- self$run_finemapping_w_rescue(dataset, lead_variant, label)
            susie_result <- fm$susie_result
            dataset <- fm$dataset
            out$dataset <- dataset
            out$susie_result <- susie_result
            out$attempted <- TRUE

            # Check the s_rss of the dataset that was actually fit (rescued value if a rescue was used)
            effective_s_rss <- dataset$s_rss_rescued %||% dataset$s_rss
            if (!is.null(susie_result$fit) && !is.na(effective_s_rss) && effective_s_rss >= self$max_s_rss) {
                .log("[%s] Discarding SuSiE fit due to unreliable LD: s_rss=%.3f >= %.2f", label, effective_s_rss, self$max_s_rss)
                out$status <- sprintf("Skipped due to unreliable LD: estimated s_rss=%.3f >= %.2f", effective_s_rss, self$max_s_rss)
                # Do not record or cache the discarded fit
                out$susie_result <- list(fit = NULL, reason = out$status, mode = susie_result$mode)
                return(out)
            }

            n_signals <- tryCatch(length(susie_result$fit$sets$cs), error = function(e) 0L) %||% 0L
            out$n_signals <- n_signals
            if (is.null(susie_result$fit)) {
                out$status <- susie_result$reason %||% "Failed, check log for details"
                return(out)
            }
            if (n_signals < 1) {
                .log("  [%s] SuSiE found no credible sets -> nothing to colocalize", label)
                out$status <- "SuSiE found no credible sets"
                return(out)
            }

            # Report success (the fit is saved to fit_fp by record_locus_results)
            out$fit <- susie_result$fit
            out$ran <- TRUE
            out$status <- "Succeeded and saved"
            return(out)
        },


        record_locus_results = function(locus_label, ran_susie, susie_result, susie_status, dataset, n_variants, n_rare_variants, lead_variant, lead_variant_pval, n_variants_ld, source) {

            # Plots results for a single GWAS or QTL fine-mapping run, saves the SuSiE fit, and updates the variant info and result details tables for that source.
            # Overrides Finemapper$record_locus_results so GWAS and QTL results go to separate directories and tables. Each locus label is recorded only once.

            # Args:
            #     locus_label: Label for the dataset (used for logging and file names).
            #     ran_susie: TRUE if SuSiE was attempted (or a precomputed fit loaded), FALSE if skipped.
            #     susie_result: List with fit, reason, and mode from the final SuSiE run.
            #     susie_status: Status string ("Succeeded and saved", failure reason, or "Skipped ...").
            #     dataset: The dataset used for the final SuSiE fit (may be NULL or a too_few list for skipped datasets).
            #     n_variants: Number of variants in the summary statistics for the locus.
            #     n_rare_variants: Number of rare variants in the locus.
            #     lead_variant: ID of the lead variant (smallest p-value).
            #     lead_variant_pval: P-value of the lead variant.
            #     n_variants_ld: Number of variants in the final dataset (after LD alignment).
            #     source: "gwas" or "qtl". Selects the results directory and tables.

            # Returns:
            #     None. Updates the source specific tables and writes plots/RDS files.

            results_dir <- self$susie_dir_for(source)
            variant_info_field <- paste0(source, "_variant_info_dt")
            result_details_field <- paste0(source, "_result_details_dt")

            # GWAS results are repeated for every QTL stratum, so only record each label once
            if (locus_label %in% self[[result_details_field]]$locus) return(invisible(NULL))

            # Initialize results variables
            susie_fit <- susie_result$fit
            fn_locus_label <- gsub("[ /]", "_", locus_label)
            susie_fit_fp <- file.path(results_dir, paste0(fn_locus_label, "_susie_fit.rds"))
            n_credible_sets <- 0L
            credible_set_lengths <- integer(0)
            lead_in_credible_sets <- FALSE

            if (ran_susie) {

                # Plot -log10 P-values (Manhattan plot)
                tryCatch({
                    png(file.path(results_dir, paste0(fn_locus_label, "_manhattan_plot.png")), width = 1200, height = 600, res = 150)
                    plot(
                        x = seq_along(dataset$D$pval),
                        y = -log10(dataset$D$pval),
                        xlab = "variable",
                        ylab = expression(-log[10](P)),
                        pch = 16, cex = 0.5
                    )
                }, error = function(e) .log("  [%s] Manhattan plot failed: %s", locus_label, conditionMessage(e)),
                finally = if (dev.cur() > 1) dev.off())

                # Plot PIP scores, append locus variant info to the variant info table, and save SuSiE fit
                if (!is.null(susie_fit)) {

                    # Plot PIP Scores
                    tryCatch({
                        png(file.path(results_dir, paste0(fn_locus_label, "_PIP_plot.png")), width = 1200, height = 600, res = 150)
                        susie_plot(susie_fit, y = "PIP", add_legend = TRUE)
                    }, error = function(e) .log("  [%s] PIP plot failed: %s", locus_label, conditionMessage(e)),
                    finally = if (dev.cur() > 1) dev.off())

                    # Extract credible sets and their lengths
                    credible_sets <- susie_fit$sets$cs
                    n_credible_sets <- length(credible_sets)
                    credible_set_lengths <- sapply(credible_sets, length)

                    # Create a data.table of variant level results (z matched by variant ID since removed variants change the dataset order)
                    locus_variants_dt <- data.table(
                        locus = locus_label,
                        variant_idx = seq_along(susie_fit$pip),
                        variant_id = names(susie_fit$pip),
                        pip = as.numeric(susie_fit$pip),
                        z = dataset$D$z[match(names(susie_fit$pip), dataset$D$snp)],
                        cs_id = NA_character_
                    )

                    # Fill in the credible set IDs for each variant
                    for (cs_name in names(credible_sets)) {
                        locus_variants_dt[variant_idx %in% credible_sets[[cs_name]], cs_id := cs_name]
                        if (lead_variant %in% names(credible_sets[[cs_name]])) {
                            lead_in_credible_sets <- TRUE
                        }
                    }

                    self[[variant_info_field]] <- rbindlist(list(self[[variant_info_field]], locus_variants_dt))

                    # Save the SuSiE fit (also used as the precomputed fit cache), unless one already exists
                    if (!file.exists(susie_fit_fp)) {
                        .log("[%s] Saving SuSiE fit to %s", locus_label, susie_fit_fp)
                        saveRDS(list(
                            fit = susie_fit, reason = susie_status, mode = susie_result$mode,
                            rescue_mode = dataset$rescue_mode, s_rss_rescued = dataset$s_rss_rescued
                        ), susie_fit_fp)
                    }
                }
            }

            # Append this dataset's details as a new row in the result details table
            self[[result_details_field]] <- rbindlist(list(
                self[[result_details_field]],
                data.table(
                    locus = locus_label,
                    n_variants = n_variants,
                    n_rare_variants = n_rare_variants,
                    lead_variant = lead_variant,
                    lead_variant_pval = lead_variant_pval,
                    n_variants_ld = n_variants_ld,
                    n_credible_sets = n_credible_sets,
                    credible_set_lengths = paste(credible_set_lengths, collapse = ","),
                    lead_in_credible_sets = lead_in_credible_sets,
                    susie_status = susie_status %||% NA_character_,
                    susie_mode = susie_result$mode %||% NA_character_,
                    rescue_mode = dataset$rescue_mode %||% NA_character_,
                    susie_estimated_s_rss = dataset$s_rss %||% NA_real_,
                    susie_estimated_s_rss_rescued = dataset$s_rss_rescued %||% NA_real_
                )
            ))

            # Clean up memory by running garbage collection
            gc(verbose = FALSE)
        },


        record_finemapping = function(fm, label, source) {

            # Convenience wrapper that passes a run_finemapping_for_coloc result to record_locus_results.

            self$record_locus_results(
                locus_label = label,
                ran_susie = fm$attempted,
                susie_result = fm$susie_result,
                susie_status = fm$status,
                dataset = fm$dataset,
                n_variants = fm$n_variants,
                n_rare_variants = fm$n_rare_variants,
                lead_variant = fm$lead_variant,
                lead_variant_pval = fm$lead_variant_pval,
                n_variants_ld = fm$n_ld,
                source = source
            )
        },


        extract_credset_for_row = function(coloc_result, row_k) {

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


        run_susie_coloc = function(gwas_susie_fit, qtl_susie_fit, gwas_stratum, qtl_stratum, locus, gene_id, susie_coloc_full_results_fp) {

            # Runs SuSiE colocalization analysis between GWAS and QTL datasets, logs the results, and writes the summary to a specified file.

            # Args:
            #     gwas_susie_fit: The SuSiE fit result for the GWAS dataset.
            #     qtl_susie_fit: The SuSiE fit result for the QTL dataset.
            #     gwas_stratum: The stratum identifier for the GWAS dataset.
            #     qtl_stratum: The stratum identifier for the QTL dataset.
            #     locus: A list containing locus information (e.g., LOCUS_ID, self$standardized_chr_key, self$manifest_left_bound_key, self$manifest_right_bound_key).
            #     gene_id: The gene identifier for the QTL dataset.
            #     susie_coloc_full_results_fp: The file path to write the colocalization summary results.

            # Returns:
            #     A logical value indicating whether any results were written to the specified file (TRUE if results were written, FALSE otherwise).

            .log("Running SuSiE colocalization...")

            # Initialize a flag to track whether any results were written to the specified file
            wrote_any <- FALSE

            # Run coloc.susie with the provided GWAS and QTL SuSiE fits, handling any errors
            coloc_result <- tryCatch(coloc.susie(
                gwas_susie_fit, qtl_susie_fit,
                p1 = self$coloc_priors$p1,
                p2 = self$coloc_priors$p2,
                p12 = self$coloc_priors$p12
            ), error = function(e) { .log("  coloc.susie error: %s", conditionMessage(e)); NULL })

            # Check if the coloc result is valid and contains a summary before proceeding to write the results
            if (!is.null(coloc_result) && !is.null(coloc_result$summary) && nrow(coloc_result$summary)) {

                # Convert the coloc summary to a data frame for easier manipulation
                coloc_summary <- as.data.frame(coloc_result$summary)

                # Determine which PP columns exist and whether any row has a finite PP value
                pp_cols <- intersect(c("PP.H0.abf","PP.H1.abf","PP.H2.abf","PP.H3.abf","PP.H4.abf"), names(coloc_summary))
                has_finite_pp <- FALSE
                if (length(pp_cols) == 5) {
                    # Coerce to numeric matrix then ensure every element is finite
                    pp_mat <- suppressWarnings(as.matrix(sapply(coloc_summary[, pp_cols, drop = FALSE], as.numeric)))
                    has_finite_pp <- all(is.finite(pp_mat))
                } else {
                    .log("  coloc.susie summary missing PP columns (%s)", paste(pp_cols, collapse = ","))
                }
                # Log a message if coloc.susie produced a summary but no finite posterior probabilities
                if (!has_finite_pp) {
                    .log("  coloc.susie produced summary but no finite posterior probabilities (likely due to credible set below coverage threshold)")
                } else {

                    # Loop through each row of the coloc summary and extract the credible set information for that row
                    for (result_index in seq_len(nrow(coloc_summary))) {

                        # Extract posterior probabilities for the current result
                        pp_h0 <- suppressWarnings(as.numeric(coloc_summary$PP.H0.abf[result_index]))
                        pp_h1 <- suppressWarnings(as.numeric(coloc_summary$PP.H1.abf[result_index]))
                        pp_h2 <- suppressWarnings(as.numeric(coloc_summary$PP.H2.abf[result_index]))
                        pp_h3 <- suppressWarnings(as.numeric(coloc_summary$PP.H3.abf[result_index]))
                        pp_h4 <- suppressWarnings(as.numeric(coloc_summary$PP.H4.abf[result_index]))

                        # Log the posterior probabilities
                        .log("  PP.H0.susie PP.H1.susie PP.H2.susie PP.H3.susie PP.H4.susie")
                        .log("  %g %g %g %g %g", pp_h0, pp_h1, pp_h2, pp_h3, pp_h4)
                        .log("  [SuSiE colocalization for %s] PP.H4 = %.2f%%", gene_id, pp_h4 * 100)

                        # Extract the credible set information for the current row of the coloc summary
                        credible_set <- self$extract_credset_for_row(coloc_result, result_index)

                        # Write the colocalization summary results to the specified file path in a streaming manner, appending to the file if it already exists
                        .stream_write(tibble(
                            !!self$gwas_strata_key := gwas_stratum,
                            !!self$qtl_strata_key := qtl_stratum,
                            !!self$manifest_locus_id_key := locus[[self$manifest_locus_id_key]],
                            ld_panel_note = locus[[self$manifest_note_key]],
                            gene_id = gene_id,
                            method = "susie.coloc",
                            hit1 = as.character(coloc_summary$hit1[result_index]),
                            hit2 = as.character(coloc_summary$hit2[result_index]),
                            nsnps = coloc_summary$nsnps[result_index] %||% NA_integer_,
                            PP.H0 = coloc_summary$PP.H0.abf[result_index],
                            PP.H1 = coloc_summary$PP.H1.abf[result_index],
                            PP.H2 = coloc_summary$PP.H2.abf[result_index],
                            PP.H3 = coloc_summary$PP.H3.abf[result_index],
                            PP.H4 = coloc_summary$PP.H4.abf[result_index],
                            cred_set_size_95 = credible_set$size,
                            top_snp_h4 = credible_set$top,
                        ), susie_coloc_full_results_fp)

                        # Update the bucket counter for SuSiE pairs
                        self$bucket["susie_pairs"] <- self$bucket["susie_pairs"] + 1L
                    }
                    # Set the wrote_any flag to TRUE to indicate that results have been written
                    wrote_any <- TRUE

                    .log("SuSiE colocalization results written for gene %s in locus %s", gene_id, locus[[self$manifest_locus_id_key]])
                }
            }

            # Return the wrote_any flag to indicate whether any results were written to the specified file
            return(wrote_any)
        },


        run_coloc_abf = function(gwas_table, gwas_stratum, qtl_table, qtl_stratum, gwas_sample_size, gwas_s, qtl_sample_size, qtl_sdY, locus, gene_id, coloc_abf_full_results_fp) {

            # Runs coloc.abf analysis between GWAS and QTL datasets, filtering for shared SNPs and handling errors.

            # Args:
            #     gwas_table: A data frame containing the GWAS dataset with standardized column names for SNPs, beta, variance of beta, and EAF.
            #     gwas_stratum: A string representing the stratum identifier for the GWAS dataset.
            #     qtl_table: A data frame containing the QTL dataset with standardized column names for SNPs, beta, variance of beta, and EAF.
            #     qtl_stratum: A string representing the stratum identifier for the QTL dataset.
            #     gwas_sample_size: An integer representing the sample size for the GWAS dataset.
            #     gwas_s: The case fraction for the GWAS dataset.
            #     qtl_sample_size: An integer representing the sample size for the QTL dataset.
            #     qtl_sdY: An optional numeric value representing the standard deviation of the trait for the QTL dataset.
            #     locus: A list containing locus information (e.g., LOCUS_ID, self$standardized_chr_key, self$manifest_left_bound_key, self$manifest_right_bound_key).
            #     gene_id: The gene identifier for the QTL dataset.
            #     coloc_abf_full_results_fp: The file path to write the colocalization summary results.

            # Returns:
            #     A logical value indicating whether any results were written to the specified file (TRUE if results were written, FALSE otherwise).

            .log("Running coloc.abf...")

            # Initialize write_any flag to FALSE to track if any results were written
            wrote_any <- FALSE

            # Prepare the GWAS and QTL datasets by selecting relevant columns and renaming them for compatibility with coloc.abf analysis
            t1 <- gwas_table %>% transmute(
                snp = .data[[self$standardized_variant_id_key]],
                gwas_beta = .data[[self$standardized_beta_key]],
                gwas_variance = .data[[self$standardized_var_beta_key]],
                gwas_eaf = .data[[self$standardized_eaf_key]]
            )

            t2 <- qtl_table %>% transmute(
                snp = .data[[self$standardized_variant_id_key]],
                qtl_beta = .data[[self$standardized_beta_key]],
                qtl_variance = .data[[self$standardized_var_beta_key]],
                qtl_eaf = .data[[self$standardized_eaf_key]]
            )

            # Identify shared SNPs between the GWAS and QTL datasets, filtering for finite beta and variance values
            raw_shared <- inner_join(t1 %>% select(snp), t2 %>% select(snp), by = "snp")

            shared <- inner_join(t1, t2, by = "snp") %>%
                filter(is.finite(gwas_beta), is.finite(gwas_variance), gwas_variance > 0,
                       is.finite(qtl_beta), is.finite(qtl_variance), qtl_variance > 0)

            .log("  coloc.abf: shared raw=%d filtered=%d (min_overlap=%d)", nrow(raw_shared), nrow(shared), self$min_overlap)

            # Check if the number of shared SNPs is less than the minimum overlap required for ABF analysis
            if (nrow(shared) < self$min_overlap) {
                .log("  coloc.abf: too few shared SNPs (%d < %d) -> skip", nrow(shared), self$min_overlap)
                return(FALSE)
            }

            # Collapse duplicate SNPs by calculating the weighted average of beta and variance for both GWAS and QTL datasets
            if (any(duplicated(shared$snp))) {
                shared <- shared %>% group_by(snp) %>% summarise(
                gwas_beta = sum(gwas_beta/gwas_variance)/sum(1/gwas_variance), gwas_variance = 1/sum(1/gwas_variance),
                qtl_beta = sum(qtl_beta/qtl_variance)/sum(1/qtl_variance), qtl_variance = 1/sum(1/qtl_variance),
                gwas_eaf = mean(gwas_eaf, na.rm = TRUE), qtl_eaf = mean(qtl_eaf, na.rm = TRUE), .groups = "drop")
            }

            # Prepare the datasets for coloc.abf analysis
            d1 <- list(snp = shared$snp, beta = shared$gwas_beta, varbeta = shared$gwas_variance,
                        N = gwas_sample_size, s = gwas_s, type = "cc",
                        EAF = ifelse(is.finite(shared$gwas_eaf), shared$gwas_eaf, 0.5))

            d2 <- list(snp = shared$snp, beta = shared$qtl_beta, varbeta = shared$qtl_variance,
                        N = qtl_sample_size, type = "quant", sdY = qtl_sdY,
                        EAF = ifelse(is.finite(shared$qtl_eaf), shared$qtl_eaf, 0.5))

            # Run coloc.abf analysis with the prepared datasets and specified priors, handling any errors
            coloc_result <- tryCatch(coloc.abf(d1, d2, p1 = self$coloc_priors$p1, p2 = self$coloc_priors$p2, p12 = self$coloc_priors$p12),
                            error = function(e) NULL)

            # Check if the coloc result is NULL
            if (is.null(coloc_result)) {
                .log("  coloc.abf failed -> skip")
                return(FALSE)
            }

            # Extract the summary from the coloc result
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

            # Check if the coloc.abf result contains finite posterior probabilities before writing the results
            if (!is.null(coloc_abf_result$PP) && all(is.finite(coloc_abf_result$PP))) {

                .stream_write(tibble(
                    !!self$gwas_strata_key := gwas_stratum,
                    !!self$qtl_strata_key := qtl_stratum,
                    !!self$manifest_locus_id_key := locus[[self$manifest_locus_id_key]],
                    ld_panel_note = locus[[self$manifest_note_key]],
                    gene_id = gene_id,
                    method = "coloc.abf",
                    nsnps = coloc_abf_result$n,
                    PP.H0 = coloc_abf_result$PP[1],
                    PP.H1 = coloc_abf_result$PP[2],
                    PP.H2 = coloc_abf_result$PP[3],
                    PP.H3 = coloc_abf_result$PP[4],
                    PP.H4 = coloc_abf_result$PP[5],
                    top_snp_h4 = coloc_abf_result$top,
                ), coloc_abf_full_results_fp)

                wrote_any <- TRUE
                self$bucket["coloc_abf_pairs"] <- self$bucket["coloc_abf_pairs"] + 1L
                .log("Coloc.abf results written for gene %s in locus %s", gene_id, locus[[self$manifest_locus_id_key]])
            }

            # Return the wrote_any flag to indicate whether any results were written to the specified file
            return(wrote_any)
        },


        process_coloc_locus = function(locus, gwas_sub, gwas_stratum, gwas_sample_size, gwas_s, qtl_sub, qtl_stratum, qtl_sample_size, qtl_sdY, coloc_summary_fp, susie_coloc_full_results_fp, coloc_abf_full_results_fp) {

            # Processes a single locus for colocalization between GWAS and QTL datasets. Fine-mapping (dataset building, QC, cascading rescue, SuSiE) is inherited from Finemapper.
            # (Named process_coloc_locus so it does not override the inherited Finemapper$process_locus, which has a different signature.)

            # Args:
            #     locus: A manifest row containing locus information.
            #     gwas_sub, gwas_stratum, gwas_sample_size, gwas_s: GWAS table for the current stratum and its parameters.
            #     qtl_sub, qtl_stratum, qtl_sample_size, qtl_sdY: QTL table for the current stratum and its parameters.
            #     coloc_summary_fp: File path for the per gene/locus summary.
            #     susie_coloc_full_results_fp: File path for the full SuSiE colocalization results.
            #     coloc_abf_full_results_fp: File path for the full coloc.abf results.

            # Returns:
            #     None. The function writes results to the specified file paths and updates internal state.

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

            # Define locus_gwas_label
            locus_gwas_label <- sprintf("%s_%s_%s", self$gwas_label, gwas_stratum, locus[[self$manifest_locus_id_key]])

            # Count the number of variants in the GWAS SNPs in locus
            .log("  [%s] Variants in GWAS locus: %d", locus_gwas_label, nrow(locus_gwas))

            # Load the LD matrix (inherited), skipping the MHC region if configured
            locus_in_mhc <- self$in_mhc(locus)
            if (locus_in_mhc) {
                .log("Skipping SuSiE for locus %s: in MHC region", locus[[self$manifest_locus_id_key]])
                ld <- NULL
            } else {
                ld <- self$load_locus_ld(locus)
            }

            # Build the GWAS dataset, fine-map it with the inherited cascading rescue, and record the results in the GWAS tables
            gwas_fm <- self$run_finemapping_for_coloc(
                tbl = locus_gwas, ld = ld, type = "cc", N = gwas_sample_size, s = gwas_s,
                target_allele_key = self$gwas_target_allele_key, label = locus_gwas_label, locus_in_mhc = locus_in_mhc,
                source = "gwas"
            )
            self$record_finemapping(gwas_fm, locus_gwas_label, "gwas")
            gwas_susie_fit <- gwas_fm$fit

            # Extract gene ids from the filtered QTL dataset
            gene_ids <- unique(locus_qtl[[self$standardized_gene_id_key]])

            # Loop through each gene id and perform colocalization analysis
            for (gene_id in gene_ids) {

                # Reset coloc flags for the current gene
                ran_susie_coloc <- FALSE
                ran_coloc_abf <- FALSE

                # Filter the QTL dataset for the current gene id
                gene_qtl <- locus_qtl %>% filter(.data[[self$standardized_gene_id_key]] == gene_id)

                # Define gene_qtl_label
                gene_qtl_label <- sprintf("%s_%s_%s_%s", self$qtl_label, qtl_stratum, locus[[self$manifest_locus_id_key]], gene_id)

                # Count the number of variants in the QTL SNPs in locus
                .log("  [%s] Variants in QTL locus: %d", gene_qtl_label, nrow(gene_qtl))

                # Check GWAS and QTL variant overlap
                n_overlap <- length(intersect(locus_gwas[[self$standardized_variant_id_key]], gene_qtl[[self$standardized_variant_id_key]]))
                if (n_overlap < self$min_overlap) {

                    # Write a summary of the low overlap gene locus pair (fine-mapping details are in the SuSiE results tables)
                    .stream_write(tibble(
                        !!self$manifest_locus_id_key := locus[[self$manifest_locus_id_key]],
                        gwas_locus = locus_gwas_label,
                        qtl_locus = gene_qtl_label,
                        n_overlap = n_overlap,
                        ld_panel_note = locus[[self$manifest_note_key]],
                        ran_susie_coloc = ran_susie_coloc,
                        ran_susie_coloc_reason = paste0("overlap<", self$min_overlap),
                        ran_coloc_abf = ran_coloc_abf,
                        ran_coloc_reason = paste0("overlap<", self$min_overlap)
                    ), coloc_summary_fp)

                    # Update the bucket counter for low overlap and skip this gene
                    self$bucket["low_overlap"] <- self$bucket["low_overlap"] + 1L
                    .log("Skipping coloc gene %s in locus %s: insufficient overlap (n_overlap=%s)", gene_id, locus[[self$manifest_locus_id_key]], n_overlap)
                    next
                }

                # Build the QTL dataset, fine-map it with the inherited cascading rescue, and record the results in the QTL tables
                qtl_fm <- self$run_finemapping_for_coloc(
                    tbl = gene_qtl, ld = ld, type = "quant", N = qtl_sample_size, sdY = qtl_sdY,
                    target_allele_key = self$qtl_target_allele_key, label = gene_qtl_label, locus_in_mhc = locus_in_mhc,
                    source = "qtl"
                )
                self$record_finemapping(qtl_fm, gene_qtl_label, "qtl")
                qtl_susie_fit <- qtl_fm$fit

                # Check if both GWAS and QTL SuSiE fits are available before running SuSiE colocalization analysis
                if (!is.null(gwas_susie_fit) && !is.null(qtl_susie_fit)) {

                    # Run SuSiE colocalization analysis and write results to the specified file path
                    wrote_any <- self$run_susie_coloc(
                        gwas_susie_fit = gwas_susie_fit,
                        qtl_susie_fit = qtl_susie_fit,
                        gwas_stratum = gwas_stratum,
                        qtl_stratum = qtl_stratum,
                        locus = locus,
                        gene_id = gene_id,
                        susie_coloc_full_results_fp = susie_coloc_full_results_fp
                    )

                    if (wrote_any) {
                        ran_susie_coloc <- TRUE
                        ran_susie_coloc_reason <- NA_character_
                    } else {
                        ran_susie_coloc_reason <- "SuSiE colocalization produced summary but no finite posterior probabilities (likely due to insufficient overlap or other issues)"
                        .log("No SuSiE colocalization results written for gene %s in locus %s: %s", gene_id, locus[[self$manifest_locus_id_key]], ran_susie_coloc_reason)
                    }
                } else {
                    ran_susie_coloc_reason <- sprintf("Missing SuSiE fits (GWAS: %s, QTL: %s)", !is.null(gwas_susie_fit), !is.null(qtl_susie_fit))
                    .log("Skipping SuSiE colocalization for gene %s in locus %s: %s", gene_id, locus[[self$manifest_locus_id_key]], ran_susie_coloc_reason)
                }

                # Run coloc.abf analysis
                wrote_any <- self$run_coloc_abf(
                    gwas_table = locus_gwas,
                    gwas_stratum = gwas_stratum,
                    qtl_table = gene_qtl,
                    qtl_stratum = qtl_stratum,
                    gwas_sample_size = gwas_sample_size,
                    gwas_s = gwas_s,
                    qtl_sample_size = qtl_sample_size,
                    qtl_sdY = qtl_sdY,
                    locus = locus,
                    gene_id = gene_id,
                    coloc_abf_full_results_fp = coloc_abf_full_results_fp
                )

                if (wrote_any) {
                    ran_coloc_abf <- TRUE
                    ran_coloc_abf_reason <- NA_character_
                } else {
                    ran_coloc_abf_reason <- "coloc.abf produced summary but no finite posterior probabilities (likely due to insufficient overlap or other issues)"
                    .log("No coloc.abf results written for gene %s in locus %s: %s", gene_id, locus[[self$manifest_locus_id_key]], ran_coloc_abf_reason)
                }

                # Write a summary of the gene locus pair (fine-mapping details are in the SuSiE results tables)
                .stream_write(tibble(
                    !!self$manifest_locus_id_key := locus[[self$manifest_locus_id_key]],
                    gwas_locus = locus_gwas_label,
                    qtl_locus = gene_qtl_label,
                    n_overlap = n_overlap,
                    ld_panel_note = locus[[self$manifest_note_key]],
                    ran_susie_coloc = ran_susie_coloc,
                    ran_susie_coloc_reason = ran_susie_coloc_reason %||% NA_character_,
                    ran_coloc_abf = ran_coloc_abf,
                    ran_coloc_reason = ran_coloc_abf_reason %||% NA_character_
                ), coloc_summary_fp)
            }

            # Log that processing of the locus for the current GWAS and QTL strata has finished
            .log("Finished processing locus %s for GWAS stratum %s and QTL stratum %s", locus[[self$manifest_locus_id_key]], gwas_stratum, qtl_stratum)

            # Cleanup: remove large objects to free memory
            rm(list = intersect(c("ld", "locus_gwas", "locus_qtl", "gwas_susie_fit", "qtl_susie_fit", "gwas_fm", "qtl_fm"), ls()))
            gc(verbose = FALSE, full = TRUE)

            # Flush the log to ensure all messages are written to the log file
            flush(stderr())
            flush(stdout())
        },


        run_coloc_analysis = function() {

            # Runs colocalization analysis between GWAS and QTL datasets for the loci in the LD manifest, processing each locus and recording results.

            # Returns:
            #     None. Results are written to the output directory.

            # Load GWAS, QTL, and loci data from the specified file paths
            .log("Loading GWAS Summary Statistics...")
            gwas_data <- .load_table(self$gwas_path)
            .log("GWAS Summary Statistics loaded: %d rows", nrow(gwas_data))
            .log("Loading QTL Summary Statistics...")
            qtl_data <- .load_table(self$qtl_path)
            .log("QTL Summary Statistics loaded: %d rows", nrow(qtl_data))
            .log("Loading LD Manifest...")
            loci_data <- .load_table(self$ld_manifest_path)
            .log("LD Manifest loaded: %d rows", nrow(loci_data))

            # Normalize chromosome representation across all tables (e.g., "chr1" vs "1")
            gwas_data[[self$standardized_chr_key]] <- .norm_chr(gwas_data[[self$standardized_chr_key]])
            qtl_data[[self$standardized_chr_key]] <- .norm_chr(qtl_data[[self$standardized_chr_key]])
            loci_data[[self$standardized_chr_key]] <- .norm_chr(loci_data[[self$standardized_chr_key]])

            # Add strata columns and fill values if necessary for GWAS and QTL datasets
            if (is.null(self$gwas_strata_key)) {
                self$gwas_strata_key <- "GWAS_Strata"
                gwas_data[[self$gwas_strata_key]] <- "Full"
            } else if (!self$gwas_strata_key %in% colnames(gwas_data)) {
                gwas_data[[self$gwas_strata_key]] <- "Full"
            } else if (gwas_data[[self$gwas_strata_key]] %>% is.na() %>% any()) {
                gwas_data[[self$gwas_strata_key]] <- as.character(gwas_data[[self$gwas_strata_key]])
                gwas_data[[self$gwas_strata_key]] <- gwas_data[[self$gwas_strata_key]] %>% replace_na("Full")
            }

            if (is.null(self$qtl_strata_key)) {
                self$qtl_strata_key <- "QTL_Strata"
                qtl_data[[self$qtl_strata_key]] <- "Bulk"
            } else if (!self$qtl_strata_key %in% colnames(qtl_data)) {
                qtl_data[[self$qtl_strata_key]] <- "Bulk"
            } else if (qtl_data[[self$qtl_strata_key]] %>% is.na() %>% any()) {
                qtl_data[[self$qtl_strata_key]] <- as.character(qtl_data[[self$qtl_strata_key]])
                qtl_data[[self$qtl_strata_key]] <- qtl_data[[self$qtl_strata_key]] %>% replace_na("Bulk")
            }

            # Extract unique strata identifiers from the GWAS and QTL datasets
            gwas_strata <- unique(gwas_data[[self$gwas_strata_key]])
            qtl_strata <- unique(qtl_data[[self$qtl_strata_key]])

            # Loop through each combination of GWAS and QTL strata
            for (gwas_stratum in gwas_strata) {

                # Define file paths for low overlap and full colocalization results based on the current GWAS stratum
                coloc_summary_fp <- file.path(self$qc_dir, sprintf("%s_%s_%s_coloc_summary.tsv", tolower(self$gwas_label), tolower(gwas_stratum), tolower(self$qtl_label)))
                susie_coloc_full_results_fp <- file.path(self$output_dir, sprintf("%s_%s_%s_susie_coloc_all_results.tsv", tolower(self$gwas_label), tolower(gwas_stratum), tolower(self$qtl_label)))
                coloc_abf_full_results_fp <- file.path(self$output_dir, sprintf("%s_%s_%s_coloc_abf_all_results.tsv", tolower(self$gwas_label), tolower(gwas_stratum), tolower(self$qtl_label)))

                # Ensure fresh outputs for this run: remove prior files so .stream_write will write headers
                if (file.exists(coloc_summary_fp)) file.remove(coloc_summary_fp)
                if (file.exists(susie_coloc_full_results_fp)) file.remove(susie_coloc_full_results_fp)
                if (file.exists(coloc_abf_full_results_fp)) file.remove(coloc_abf_full_results_fp)

                for (qtl_stratum in qtl_strata) {
                    .log("Running colocalization analysis for GWAS stratum %s and QTL stratum %s", gwas_stratum, qtl_stratum)

                    # Filter the GWAS and QTL datasets for the current strata
                    gwas_sub <- gwas_data %>% filter(!!sym(self$gwas_strata_key) == gwas_stratum)
                    qtl_sub <- qtl_data %>% filter(!!sym(self$qtl_strata_key) == qtl_stratum)

                    # Loop through each locus and process it
                    for (i in seq_len(nrow(loci_data))) {
                        locus <- loci_data[i, ]
                        self$process_coloc_locus(locus, gwas_sub, gwas_stratum, self$gwas_sample_size, self$gwas_s,
                            qtl_sub, qtl_stratum, self$qtl_sample_size, self$qtl_sdY, coloc_summary_fp, susie_coloc_full_results_fp, coloc_abf_full_results_fp)
                    }
                }

                # Write Report
                report_path <- file.path(self$output_dir, sprintf("%s_%s_%s_summary_report.txt", tolower(self$gwas_label), tolower(gwas_stratum), tolower(self$qtl_label)))
                report <- c(
                    sprintf("coloc summary report — %s", gwas_stratum),
                    sprintf("Generated: %s", format(Sys.time())),
                    sprintf("%s: %s", self$gwas_strata_key, gwas_stratum),
                    "",
                    sprintf("Signal pairs written (susie): %d", self$bucket["susie_pairs"]),
                    sprintf("Signal pairs written (coloc.abf): %d", self$bucket["coloc_abf_pairs"]),
                    sprintf("Loci w/ overlap below %g: %d", self$min_overlap, self$bucket["low_overlap"]),
                    "",
                    sprintf("Priors: p1 = %g, p2 = %g, p12 = %g", self$coloc_priors$p1, self$coloc_priors$p2, self$coloc_priors$p12),
                    sprintf("SuSiE: coverage = %.2f, max_iter = %g, L = %g, max s_rss = %.2f", self$cred_coverage, self$susie_max_iter, self$susie_l, self$max_s_rss),
                    sprintf("MIN OVERLAP = %g, MIN SNPS SUSIE = %g, DROP AMBIGUOUS = %s", self$min_overlap, self$susie_min_snps, self$drop_ambiguous),
                    sprintf("GWAS sample size = %g, case fraction = %.3f", self$gwas_sample_size, self$gwas_s),
                    sprintf("QTL sample size = %g, sdY = %g", self$qtl_sample_size, self$qtl_sdY),
                    sprintf("LD panel note: %s", unique(loci_data[[self$manifest_note_key]])[1])
                )
                writeLines(report, report_path)
                .log("%s: Completed. susie_pairs = %d, abf = %d, low overlap = %d",
                    gwas_stratum, self$bucket["susie_pairs"], self$bucket["coloc_abf_pairs"], self$bucket["low_overlap"])
            }

            # Save the SuSiE variant information and result details tables, separately for GWAS and QTL
            for (source in c("gwas", "qtl")) {
                results_dir <- self$susie_dir_for(source)
                label <- if (source == "gwas") self$gwas_label else self$qtl_label
                variant_info_fp <- file.path(results_dir, sprintf("%s_susie_variant_info.tsv", label))
                result_details_fp <- file.path(results_dir, sprintf("%s_susie_result_details.tsv", label))
                .log("Saving %s variant information to %s", toupper(source), variant_info_fp)
                fwrite(self[[paste0(source, "_variant_info_dt")]], variant_info_fp, sep = "\t", na = "NA", quote = FALSE)
                .log("Saving %s result details to %s", toupper(source), result_details_fp)
                fwrite(self[[paste0(source, "_result_details_dt")]], result_details_fp, sep = "\t", na = "NA", quote = FALSE)
            }
        }
    )
)

# ------------------------------------- Command Line Interface -------------------------------------

# Only run the CLI when this file is the main script (not when source()d by another script)
if (!interactive() && sys.nframe() == 0) {

    parser <- ArgumentParser(description = "Colocalization analysis using SuSiE and coloc.abf with inherited cascading-rescue fine-mapping")
    parser$add_argument("--output_dir", required = TRUE, help = "Path to the directory to save output files")
    parser$add_argument("--ld_dir", required = TRUE, help = "Path to the directory containing LD matrices and BIM files.")
    parser$add_argument("--qc_dir", required = TRUE, help = "Path to the directory containing QC files.")
    parser$add_argument("--gwas_fp", required = TRUE, help = "Path to the standardized GWAS summary statistics file")
    parser$add_argument("--gwas_target_allele_key", default = "EFFECT", help = "Column name for the target allele in the GWAS summary statistics file (optional). Default is 'EFFECT'.")
    parser$add_argument("--gwas_strata_key", default = NULL, help = "Column name for the GWAS strata in the GWAS summary statistics file (optional).")
    parser$add_argument("--gwas_sample_size", type = "numeric", required = TRUE, help = "Sample size for the GWAS dataset")
    parser$add_argument("--gwas_case_fraction", type = "numeric", required = TRUE, help = "Case fraction for the GWAS dataset.")
    parser$add_argument("--qtl_fp", required = TRUE, help = "Path to the standardized QTL summary statistics file")
    parser$add_argument("--qtl_target_allele_key", default = "EFFECT", help = "Column name for the target allele in the QTL summary statistics file (optional). Default is 'EFFECT'.")
    parser$add_argument("--qtl_strata_key", default = NULL, help = "Column name for the QTL strata in the QTL summary statistics file (optional).")
    parser$add_argument("--qtl_sample_size", type = "numeric", required = TRUE, help = "Sample size for the QTL dataset")
    parser$add_argument("--qtl_sdY", type = "numeric", default = 1, help = "Standard deviation of the trait for the QTL dataset (optional). Default is 1.")
    parser$add_argument("--ld_manifest", default = "ld_manifest.tsv", help = "File name of the ld manifest file containing locus information (optional). Default is 'ld_manifest.tsv'.")
    parser$add_argument("--drop_ambiguous", default = TRUE, type = "logical", help = "Drop ambiguous SNPs (A/T or C/G) from the analysis (optional). Default is TRUE.")
    parser$add_argument("--min_overlap", type = "numeric", default = 50, help = "Minimum number of overlapping SNPs required for colocalization analysis (optional). Default is 50.")
    parser$add_argument("--cred_coverage", type = "numeric", default = 0.95, help = "Credible set coverage for SuSiE analysis (optional). Default is 0.95.")
    parser$add_argument("--pp_h4_threshold", type = "numeric", default = 0.70, help = "Posterior probability threshold for strong colocalization (optional). Default is 0.70.")
    parser$add_argument("--max_s_rss", type = "numeric", default = 0.7, help = "Maximum estimated s_rss (after any rescue) for which SuSiE results are used (optional). Default is 0.3.")
    parser$add_argument("--susie_min_snps", type = "numeric", default = 5, help = "Minimum number of SNPs required for SuSiE analysis (optional). Default is 5.")
    parser$add_argument("--susie_max_iter", type = "numeric", default = 1000, help = "Maximum number of iterations for SuSiE analysis (optional). Default is 1000.")
    parser$add_argument("--susie_l", type = "numeric", default = 10, help = "Maximum number of causal variants for SuSiE analysis (optional). Default is 10.")
    parser$add_argument("--susie_repeat_until_converged", default = FALSE, type = "logical", help = "Repeat SuSiE until convergence (optional). Default is FALSE.")
    parser$add_argument("--coloc_prior_p1", type = "numeric", default = 1e-4, help = "Prior probability for hypothesis H1 (optional). Default is 1e-4.")
    parser$add_argument("--coloc_prior_p2", type = "numeric", default = 1e-4, help = "Prior probability for hypothesis H2 (optional). Default is 1e-4.")
    parser$add_argument("--coloc_prior_p12", type = "numeric", default = 1e-5, help = "Prior probability for hypothesis H4 (optional). Default is 1e-5.")
    parser$add_argument("--manifest_locus_id_key", default = "LOCUS_ID", help = "Column name for the locus ID in the LD manifest file (optional). Default is 'LOCUS_ID'.")
    parser$add_argument("--manifest_left_bound_key", default = "LEFT_500KB", help = "Column name for the left bound of the locus in the LD manifest file (optional). Default is 'LEFT_500KB'.")
    parser$add_argument("--manifest_right_bound_key", default = "RIGHT_500KB", help = "Column name for the right bound of the locus in the LD manifest file (optional). Default is 'RIGHT_500KB'.")
    parser$add_argument("--manifest_ld_key", default = "LD", help = "Column name for the LD matrix file path in the LD manifest file (optional). Default is 'LD'.")
    parser$add_argument("--manifest_bim_key", default = "BIM", help = "Column name for the BIM file path in the LD manifest file (optional). Default is 'BIM'.")
    parser$add_argument("--manifest_note_key", default = "note", help = "Column name for the LD panel note in the LD manifest file (optional). Default is 'note'.")
    parser$add_argument("--standardized_variant_id_key", default = "VAR", help = "Column name for the standardized variant ID in the GWAS and QTL summary statistics files (optional). Default is 'VAR'.")
    parser$add_argument("--standardized_chr_key", default = "CHR", help = "Column name for the standardized chromosome in the GWAS and QTL summary statistics files (optional). Default is 'CHR'.")
    parser$add_argument("--standardized_pos_key", default = "BP", help = "Column name for the standardized position in the GWAS and QTL summary statistics files (optional). Default is 'BP'.")
    parser$add_argument("--standardized_gene_id_key", default = "GENE", help = "Column name for the standardized gene ID in the QTL summary statistics file (optional). Default is 'GENE'.")
    parser$add_argument("--standardized_non_effect_allele_key", default = "NON_EFFECT", help = "Column name for the standardized reference allele in the GWAS and QTL summary statistics files (optional). Default is 'NON_EFFECT'.")
    parser$add_argument("--standardized_effect_allele_key", default = "EFFECT", help = "Column name for the standardized effect allele in the GWAS and QTL summary statistics files (optional). Default is 'EFFECT'.")
    parser$add_argument("--standardized_p_key", default = "P", help = "Column name for the standardized p-value in the GWAS and QTL summary statistics files (optional). Default is 'P'.")
    parser$add_argument("--standardized_beta_key", default = "BETA", help = "Column name for the standardized beta in the GWAS and QTL summary statistics files (optional). Default is 'BETA'.")
    parser$add_argument("--standardized_var_beta_key", default = "VARBETA", help = "Column name for the standardized variance of beta in the GWAS and QTL summary statistics files (optional). Default is 'VARBETA'.")
    parser$add_argument("--standardized_se_key", default = "SE", help = "Column name for the standardized standard error in the GWAS and QTL summary statistics files (optional). Default is 'SE'.")
    parser$add_argument("--standardized_eaf_key", default = "EAF", help = "Column name for the standardized effect allele frequency in the GWAS and QTL summary statistics files (optional). Default is 'EAF'.")
    parser$add_argument("--skip_MHC", default = TRUE, type = "logical", help = "Skip loci in the MHC region (optional). Default is TRUE.")
    parser$add_argument("--gwas_susie_results_dir", default = NULL, help = "Directory to save GWAS SuSiE results (optional). Default is <output_dir>/gwas_susie_results.")
    parser$add_argument("--qtl_susie_results_dir", default = NULL, help = "Directory to save QTL SuSiE results (optional). Default is <output_dir>/qtl_susie_results.")
    parser$add_argument("--gwas_label", default = "gwas", help = "Label for the GWAS dataset (optional). Default is 'gwas'.")
    parser$add_argument("--qtl_label", default = "qtl", help = "Label for the QTL dataset (optional). Default is 'qtl'.")

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

    # Create an instance of the ColocalizationAnalyzer class
    analyzer <- ColocalizationAnalyzer$new(args)

    # Run the colocalization analysis
    analyzer$run_coloc_analysis()
}
