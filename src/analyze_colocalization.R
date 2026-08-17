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

# Logging function to print messages with timestamps and write to a log file
.log <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), " | ", sprintf(...))
  cat(msg, "\n"); write(msg, file = log_file, append = TRUE)
}

# Standardize chromosome labels
.norm_chr <- function(x) paste0("chr", str_replace(toupper(as.character(x)), "^CHR", ""))

# Build an orientation-independent variant key
.canon_key <- function(chr, pos, a1, a2) {
  a1 <- toupper(a1); a2 <- toupper(a2)
  lo <- pmin(a1, a2); hi <- pmax(a1, a2)
  paste0(.norm_chr(chr), ":", as.integer(pos), ":", lo, "_", hi)
}

# Write a data frame to a file in a streaming manner (append if file exists)
.stream_write <- function(df, path) {
  if (!file.exists(path)) data.table::fwrite(df, path, sep = "\t", quote = FALSE)
  else data.table::fwrite(df, path, sep = "\t", quote = FALSE, append = TRUE, col.names = FALSE)
}

.credset_for_row <- function(coloc_result, row_k) {
    out <- list(size = NA_integer_, top = NA_character_)
    rr <- tryCatch(as.data.frame(coloc_result$results), error = function(e) NULL)
    if (is.null(rr) || !nrow(rr) || !"snp" %in% names(rr)) return(out)
    h4cols <- grep("^SNP\\.PP\\.H4", names(rr), value = TRUE)
    if (!length(h4cols)) return(out)
    col <- if (length(h4cols) >= row_k) {
        byname <- paste0("SNP.PP.H4.row", row_k)
        if (byname %in% h4cols) byname else h4cols[row_k]
    } else h4cols[1]
    v <- suppressWarnings(as.numeric(rr[[col]]))
    if (!any(is.finite(v))) return(out)
    o <- order(v, decreasing = TRUE)
    credible_set <- cumsum(v[o]); w <- which(credible_set >= CRED_COVERAGE)[1]
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
            self$qtl <- args$qtl
            self$qtl_path <- paste0(self$output_dir, args$qtl)
            self$ld_manifest <- args$ld_manifest
            self$ld_manifest_path <- paste0(self$ld_dir, args$ld_manifest)

            # Parameters
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
        },

        # Read and align LD
        # Returns list(M = panel r matrix [canon-key dimnames], A1 = named vector of the
        # allele the r is counted on). Sign alignment to a target allele happens later.
        read_ld <- function(ld_path, bim_path) {
            if (is.na(ld_path) || is.na(bim_path) ||
                !file.exists(ld_path) || !file.exists(bim_path)) {
                .log("  LD missing: ld=%s bim=%s", as.character(ld_path), as.character(bim_path))
                return(NULL)
            }
            bim <- tryCatch(fread(bim_path, header = FALSE, data.table = FALSE),
                            error = function(e) NULL)
            if (is.null(bim) || !nrow(bim)) {
                .log("  BIM missing/empty: %s", as.character(bim_path))
                return(NULL)
            }
            # PLINK .bim: V1 chr, V2 id, V3 cM, V4 bp, V5 A1, V6 A2 ; --r counts A1
            chr <- .norm_chr(bim$V1); pos <- as.integer(bim$V4)
            a1  <- toupper(bim$V5);   a2  <- toupper(bim$V6)
            keys <- .canon_key(chr, pos, a1, a2)

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
            dimnames(M) <- list(keys, keys)
            storage.mode(M) <- "numeric"
            M[!is.finite(M)] <- 0
            list(M = M, A1 = setNames(a1, keys), A2 = setNames(a2, keys))
        }

        # Align panel LD to a TARGET allele per variant, subset to `snps`, return matrix
        # in the order of `snps`. target_allele: named (by canon key) vector of the
        # allele the dataset's beta is on (GWAS effect allele / GTEx ALT).
        align_ld <- function(ld, snps, target_allele) {
            if (is.null(ld)) return(NULL)
            present <- snps[snps %in% rownames(ld$M)]
            if (length(present) < 2) return(NULL)
            M  <- ld$M[present, present, drop = FALSE]
            a1 <- ld$A1[present]; a2 <- ld$A2[present]
            tgt <- target_allele[present]
            # Align the LD sign to the effect allele
            s <- rep(NA_real_, length(present))
            s[a1 == tgt] <- 1
            s[a2 == tgt] <- -1
            ok <- is.finite(s)                          # drop variants whose alleles don't match panel
            if (sum(ok) < 2) return(NULL)
            present <- present[ok]; M <- M[ok, ok, drop = FALSE]; s <- s[ok]
            M <- (s %o% s) * M                          # flip rows+cols for negated variants
            diag(M) <- 1
            dimnames(M) <- list(present, present)
            M
        }

        # Build coloc datasets
        build_dataset <- function(tbl, ld, type, N, s = NULL, sdY = NULL) {
            # tbl: snp(KEY), beta, varbeta, MAF, EA(target allele), OA, POS
            D0 <- tbl %>% transmute(snp = KEY, beta = BETA_USE, varbeta = VARBETA_USE,
                                    MAF = MAF, EA = EA, OA = OA, POS = POS)
            D0 <- .collapse_keys(D0)
            if (DROP_AMBIGUOUS) D0 <- D0 %>% filter(!.is_ambiguous(EA, OA))
            if (nrow(D0) < 2) return(NULL)

            LD <- .align_ld(ld, D0$snp, setNames(D0$EA, D0$snp))
            if (is.null(LD)) return(NULL)

            keep <- D0$snp %in% rownames(LD)
            D0 <- D0[keep, , drop = FALSE]
            D0 <- D0[match(rownames(LD), D0$snp), , drop = FALSE]   # exact same order as LD
            if (nrow(D0) < MIN_SNPS_SUSIE) {
                return(list(too_few = TRUE, n = nrow(D0)))
            }

            D <- list(
                beta     = D0$beta,
                varbeta  = D0$varbeta,
                snp      = D0$snp,
                position = D0$POS,
                type     = type,
                N        = N,
                LD       = LD
            )
            if (!is.null(s))   D$s   <- s
            if (!is.null(sdY)) D$sdY <- sdY
            if (any(is.finite(D0$MAF))) D$MAF <- ifelse(is.finite(D0$MAF), D0$MAF, 0.5)
            list(D = D, n = nrow(D0))
        }

        # Function to safely run SuSiE on a dataset, handling errors and logging
        safe_runsusie <- function(D, label) {
            if (!is.null(check_dataset(D, req = "LD"))) {     # NULL means OK
                .log("  [%s] check_dataset failed (LD) -> skip SuSiE", label); return(NULL)
            }
            S <- tryCatch(
                suppressWarnings(runsusie(
                    D, coverage = self$susie_cred_coverage, max_iter = self$susie_max_iter, L = self$susie_l,
                    repeat_until_convergence = self$susie_repeat_until_converged
                )),
                error = function(e) { .log("  [%s] runsusie error: %s", label, conditionMessage(e)); NULL }
            )
            if (is.null(S)) return(NULL)
            ncs <- tryCatch(length(S$sets$cs), error = function(e) 0L)
            .log("  [%s] SuSiE credible sets: %d", label, ncs %||% 0L)
            if (is.null(ncs) || ncs < 1) return(NULL)   # nothing to colocalise
            S
        },

        # Function to run SuSiE Colocalization analysis (Multiple Causal Variants)
        run_susie_coloc <- function(self, gwas_susie_fit, qtl_susie_fit, gwas_stratum, qtl_stratum, locus, ld, gene_id, gene_name, n_overlap, n_gwas_ld, n_eqtl_ld, n_gwas_signals, n_eqtl_signals, top_gwas_variant, top_gwas_pval, lead_eqtl_id, lead_eqtl_p, full_path) {
            .log("Running SuSiE colocalization...")
            
            if (!is.null(gwas_susie_fit) && !is.null(qtl_susie_fit)) {
                coloc_result <- tryCatch(coloc.susie(gwas_susie_fit, qtl_susie_fit, p1 = self$coloc_priors$p1, p2 = self$coloc_priors$p2, p12 = self$coloc_priors$p12),
                                error = function(e) { .log("  coloc.susie error: %s", conditionMessage(e)); NULL })
                if (!is.null(coloc_result) && !is.null(coloc_result$summary) && nrow(coloc_result$summary)) {
                    coloc_summary <- as.data.frame(coloc_result$summary)
                    for (result_index in seq_len(nrow(coloc_summary))) {
                        credible_set <- .credset_for_row(coloc_result, result_index)
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
                        bucket["susie_pairs"] <- bucket["susie_pairs"] + 1L
                        wrote_any <- TRUE
                    }
                }
            }
            
        }
    )
)