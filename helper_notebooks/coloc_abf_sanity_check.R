#!/usr/bin/env Rscript

# Global Variable Definitions
standardized_keys = list(
    chr = "CHR",
    pos = "BP",
    rsid = "SNP",
    var_id = "VAR",
    gene_id = "GENE",
    non_effect_allele = "NON_EFFECT",
    effect_allele = "EFFECT",
    p_value = "P",
    beta = "BETA",
    se = "SE",
    statistic = "STAT",
    maf = "MAF",
    var_beta = "VARBETA"
)

min_overlap <- 10

coloc_priors <- list(
    p1 = 1e-4,
    p2 = 1e-4,
    p12 = 1e-5
)


# Imports

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(purrr)
  library(coloc)
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

run_coloc_abf = function(gwas_table, qtl_table, gwas_sample_size, gwas_s, qtl_sample_size, qtl_sdY, standardized_keys, coloc_priors, min_overlap) {
    
    # Runs coloc.abf analysis as a fallback method for colocalization between GWAS and QTL datasets, filtering for shared SNPs and handling errors.

    # Args:
    #     gwas_table: A data frame containing the GWAS dataset with standardized column names for SNPs, beta, variance of beta, and MAF.
    #     qtl_table: A data frame containing the QTL dataset with standardized column names for SNPs, beta, variance of beta, and MAF.
    #     gwas_sample_size: An integer representing the sample size for the GWAS dataset.
    #     gwas_s: An optional numeric vector representing the standard errors of the beta estimates for the GWAS dataset.
    #     qtl_sample_size: An integer representing the sample size for the QTL dataset.
    #     qtl_sdY: An optional numeric value representing the standard deviation of the trait for the QTL dataset.
    #     standardized_keys: A list of standardized column names for the GWAS and QTL datasets.
    #     coloc_priors: A list of prior probabilities for the coloc.abf analysis.
    #     min_overlap: An integer representing the minimum number of shared SNPs required for the analysis.

    # Returns:
    #     A list containing:
    #         - n: The number of shared SNPs between the GWAS and QTL datasets.
    #         - PP: A numeric vector of posterior probabilities for hypotheses H0 to H4.
    #         - top: The SNP with the highest posterior probability for hypothesis H4.
    #     Returns NULL if there are fewer than the minimum overlap required for ABF analysis or if there are issues with the coloc.abf analysis.
    
    .log("Running coloc.abf")

    # Prepare the GWAS and QTL datasets by selecting relevant columns and renaming them for compatibility with coloc.abf analysis
    t1 <- gwas_table %>% transmute(
        snp = .data[[standardized_keys$var_id]],
        gwas_beta = .data[[standardized_keys$beta]],
        gwas_variance = .data[[standardized_keys$var_beta]],
        gwas_maf = .data[[standardized_keys$maf]]
    )

    t2 <- qtl_table %>% transmute(
        snp = .data[[standardized_keys$var_id]],
        qtl_beta = .data[[standardized_keys$beta]],
        qtl_variance = .data[[standardized_keys$var_beta]],
        qtl_maf = .data[[standardized_keys$maf]]
    )

    # Identify shared SNPs between the GWAS and QTL datasets, filtering for finite beta and variance values
    raw_shared <- inner_join(t1 %>% select(snp), t2 %>% select(snp), by = "snp")
    
    shared <- inner_join(t1, t2, by = "snp") %>%
        filter(is.finite(gwas_beta), is.finite(gwas_variance), gwas_variance > 0,
                is.finite(qtl_beta), is.finite(qtl_variance), qtl_variance > 0)

    .log("  coloc.abf: shared raw=%d filtered=%d (min_overlap=%d)", nrow(raw_shared), nrow(shared), min_overlap)

    # Check if the number of shared SNPs is less than the minimum overlap required for ABF analysis, and return NULL if so
    if (nrow(shared) < min_overlap) {
        .log("  coloc.abf: too few shared SNPs (%d < %d) -> skip", nrow(shared), min_overlap)
        return(NULL)
    }

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
    coloc_result <- tryCatch(coloc.abf(d1, d2, p1 = coloc_priors$p1, p2 = coloc_priors$p2, p12 = coloc_priors$p12),
                    error = function(e) NULL)

    # Check if the coloc result is NULL and return NULL if it is
    if (is.null(coloc_result)) {
        .log("  coloc.abf failed -> skip")
        return(NULL)
    }

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

    # Return the wrote_any flag to indicate whether any results were written to the specified file
    return(coloc_abf_result) 
}


process_locus_gene_pair = function(locus_id, locus_chr, locus_left_bound, locus_right_bound, gene_id, gwas, gwas_sample_size, gwas_s, qtl, qtl_sample_size, qtl_sdY, standardized_keys, coloc_priors, min_overlap, result_fp) {
    
    # Processes a single locus for colocalization analysis between GWAS and QTL datasets, filtering variants, running coloc.abf, and recording results.

    # Args:
    #     locus_id: A string representing the unique identifier for the locus.
    #     locus_chr: A string representing the chromosome of the locus.
    #     locus_left_bound: An integer representing the left boundary position of the locus.
    #     locus_right_bound: An integer representing the right boundary position of the locus.
    #     gene_id: A string representing the unique identifier for the gene associated with the locus
    #     gwas: A data frame containing the GWAS dataset.
    #     gwas_sample_size: An integer representing the sample size for the GWAS dataset.
    #     gwas_s: An optional numeric vector representing the standard errors of the beta estimates for the GWAS dataset.
    #     qtl: A data frame containing the QTL dataset.
    #     qtl_sample_size: An integer representing the sample size for the QTL dataset.
    #     qtl_sdY: An optional numeric value representing the standard deviation of the trait for the QTL dataset.
    #     standardized_keys: A list containing standardized column names for the analysis.
    #     coloc_priors: A list containing prior distributions for the colocalization analysis.
    #     min_overlap: An integer representing the minimum overlap required between GWAS and QTL variants.
    #     result_fp: The file path to write the colocalization results.

    # Returns:
    #     coloc_abf_result: A list containing the results of the coloc.abf analysis, including the number of shared SNPs, posterior probabilities for hypotheses H0 to H4, and the top SNP. Returns NULL if there are insufficient variants for analysis.

    # Filter the GWAS and QTL datasets to include only variants within the specified locus boundaries (chromosome and position range)
    locus_gwas <- gwas %>% filter(
        .data[[standardized_keys$chr]] == locus_chr,
        .data[[standardized_keys$pos]] >= locus_left_bound,
        .data[[standardized_keys$pos]] <= locus_right_bound
    )
    locus_qtl <- qtl %>% filter(
        .data[[standardized_keys$chr]] == locus_chr,
        .data[[standardized_keys$pos]] >= locus_left_bound,
        .data[[standardized_keys$pos]] <= locus_right_bound
    )

    # Count the number of variants in the filtered GWAS and QTL datasets, skipping the locus if either dataset has no variants
    n_gwas = nrow(locus_gwas)
    n_qtl = nrow(locus_qtl)
    .log("%s: GWAS=%s QTL=%s", locus_id, format(n_gwas, big.mark = ","), format(n_qtl, big.mark = ","))
    if (!n_gwas || !n_qtl) {
        .log("Skipping locus %s: insufficient variants (GWAS=%s, QTL=%s)", 
            locus_id, n_gwas, n_qtl)
        return(NULL)
    }

    # Record the lead GWAS variant
    top_gwas <- locus_gwas %>% arrange(.data[[standardized_keys$p_value]]) %>% slice(1)
    top_gwas_variant <- top_gwas[[standardized_keys$var_id]] %||% NA_character_
    top_gwas_pval <- top_gwas[[standardized_keys$p_value]] %||% NA_real_
        
    # Filter the QTL dataset for the current gene id
    gene_qtl  <- locus_qtl %>% filter(.data[[standardized_keys$gene_id]] == gene_id)

    # Record the lead QTL variant
    lead_qtl <- gene_qtl %>% filter(is.finite(.data[[standardized_keys$p_value]])) %>% arrange(.data[[standardized_keys$p_value]]) %>% slice(1)
    lead_qtl_variant <- lead_qtl[[standardized_keys$var_id]][1] %||% NA_character_
    lead_qtl_pval <- lead_qtl[[standardized_keys$p_value]][1] %||% NA_real_

    # Check GWAS and QTL variant overlap
    n_overlap <- length(intersect(locus_gwas[[standardized_keys$var_id]], gene_qtl[[standardized_keys$var_id]]))

    # Run Coloc ABF
    coloc_abf_result <- run_coloc_abf(
        gwas_table = locus_gwas,
        qtl_table = gene_qtl,
        gwas_sample_size = gwas_sample_size,
        gwas_s = gwas_s,
        qtl_sample_size = qtl_sample_size,
        qtl_sdY = qtl_sdY,
        standardized_keys = standardized_keys,
        coloc_priors = coloc_priors,
        min_overlap = min_overlap
    )

    if (!is.null(coloc_abf_result)) {
        # Write the coloc.abf results to the specified file path in a streaming manner, appending to the file if it already exists
        .stream_write(tibble(
                locus_id = locus_id,
                gene_id = gene_id,
                method = "coloc.abf",
                nsnps = coloc_abf_result$n,
                PP.H0 = coloc_abf_result$PP[1],
                PP.H1 = coloc_abf_result$PP[2],
                PP.H2 = coloc_abf_result$PP[3],
                PP.H3 = coloc_abf_result$PP[4],
                PP.H4 = coloc_abf_result$PP[5],
                n_snps_overlap = n_overlap,
                n_gwas = n_gwas,
                n_qtl = n_qtl,
                cred_set_size_95 = NA_integer_,
                top_snp_h4 = coloc_abf_result$top,
                top_gwas_variant = top_gwas_variant,
                top_gwas_pval = top_gwas_pval,
                lead_qtl_variant = lead_qtl_variant,
                lead_qtl_pval = lead_qtl_pval
            ), result_fp)

            # Log the successful writing of coloc.abf results for the current gene and locus
            .log("Coloc.abf results written for gene %s in locus %s", gene_id, locus_id)
    } else {
        # Log the failure of coloc.abf analysis for the current gene and locus
        .log("Coloc.abf analysis failed for gene %s in locus %s", gene_id, locus_id)
    }

    # Return the coloc.abf result
    .log("Finished processing gene %s in locus %s", gene_id, locus_id)
    return(coloc_abf_result)
}

run_coloc_analysis = function(gwas_path, qtl_path, locus_gene_pairs, gwas_sample_size, gwas_s, qtl_sample_size, qtl_sdY, output_dir, standardized_keys, coloc_priors, min_overlap) {

    # Runs colocalization analysis between GWAS and QTL datasets for specified loci, processing each locus and recording results.

    # Args:
    #     gwas_path: Path to the GWAS dataset with standardized column names for SNPs, beta, variance of beta, and MAF.
    #     qtl_path: Path to the QTL dataset with standardized column names for SNPs, beta, variance of beta, and MAF.
    #     locus_gene_pairs: A tibble containing locus/gene pair information, including locus ID, chromosome, left bound, right bound, and gene id.
    #     gwas_sample_size: An integer representing the sample size for the GWAS dataset.
    #     gwas_s: A numeric vector representing the summary statistics for the GWAS dataset.
    #     qtl_sample_size: An integer representing the sample size for the QTL dataset.
    #     qtl_sdY: A numeric vector representing the standard deviation of the QTL dataset.
    #     output_dir: The directory path to write the colocalization summary results.
    #     standardized_keys: A list containing standardized column names for GWAS and QTL datasets.
    #     coloc_priors: A list containing prior probabilities for colocalization analysis.
    #     min_overlap: An integer representing the minimum overlap required for colocalization.

    # Returns:
    #     A tibble summarizing the results of the colocalization analysis across all loci and strata.

    # Load GWAS, QTL, and loci data from the specified file paths
    .log("Loading GWAS Summary Statistics...")
    gwas_data <- .load_table(gwas_path)
    .log("GWAS Summary Statistics loaded: %d rows", nrow(gwas_data))
    .log("Loading QTL Summary Statistics...")
    qtl_data <- .load_table(qtl_path)
    .log("QTL Summary Statistics loaded: %d rows", nrow(qtl_data))

    # Define results file path
    result_fp <- file.path(output_dir, "coloc_results.tsv")

    # Ensure fresh outputs for this run: remove prior files so .stream_write will write headers
    if (file.exists(result_fp)) file.remove(result_fp)

    # Loop through each combination of GWAS and QTL strata
    for (i in seq_len(nrow(locus_gene_pairs))) {
        
        # Extract current locus-gene pair information
        locus_id <- locus_gene_pairs$id[i]
        locus_chr <- locus_gene_pairs$chr[i]
        locus_left <- locus_gene_pairs$left[i]
        locus_right <- locus_gene_pairs$right[i]
        gene_id <- locus_gene_pairs$gene[i]

        .log("Running colocalization analysis [%d/%d] for %s / %s", i, nrow(locus_gene_pairs), locus_id, gene_id)

        locus_gene_result <- process_locus_gene_pair(locus_id = locus_id,
                                                    locus_chr = locus_chr,
                                                    locus_left_bound = locus_left,
                                                    locus_right_bound = locus_right,
                                                    gene_id = gene_id,
                                                    gwas = gwas_data,
                                                    gwas_sample_size = gwas_sample_size,
                                                    gwas_s = gwas_s,
                                                    qtl = qtl_data,
                                                    qtl_sample_size = qtl_sample_size,
                                                    qtl_sdY = qtl_sdY,
                                                    standardized_keys = standardized_keys,
                                                    coloc_priors = coloc_priors,
                                                    min_overlap = min_overlap,
                                                    result_fp = result_fp)
        if (!is.null(locus_gene_result)) {
            saveRDS(locus_gene_result, file.path(output_dir, sprintf("%s_%s_locus_gene_result.rds", tolower(gene_id), tolower(locus_id))))
        }
    }

    # Log the completion of the colocalization analysis
    .log("Colocalization analysis completed for all loci and genes.")
}

# ------------------------------------- Command Line Interface -------------------------------------

if (!interactive()) {

    # Variable Definitions
    gwas_path <- "/mnt/output/output/coloc/not_validated/GP2_et_al_2025_PD_case_control_EUR_ALL_hg38_rsID.standardized.tsv"
    qtl_path <- "/mnt/output/output/coloc/amp_ad/cortex_meta/Cortex_MetaAnalysis_ROSMAP_CMC_HBCC_Mayo_cis_eQTL_release.standardized.updated.tsv"
    gwas_sample_size <- 226196
    gwas_s <- 0.207306937
    qtl_sample_size <- 1694
    qtl_sdY <- 1.0
    output_dir <- "/mnt/output/output/coloc/sanity_checks"
    locus_gene_pairs <- tibble(
        id = c('7:22570400-23872051'),
        chr = c('7'),
        left = c(22570400),
        right = c(23872051),
        gene = c('GPNMB')
    )

    # Run the colocalization analysis
    run_coloc_analysis(
        gwas_path = gwas_path,
        qtl_path = qtl_path,
        locus_gene_pairs = locus_gene_pairs,
        gwas_sample_size = gwas_sample_size,
        gwas_s = gwas_s,
        qtl_sample_size = qtl_sample_size,
        qtl_sdY = qtl_sdY,
        output_dir = output_dir,
        standardized_keys = standardized_keys,
        coloc_priors = coloc_priors,
        min_overlap = min_overlap
    )
}