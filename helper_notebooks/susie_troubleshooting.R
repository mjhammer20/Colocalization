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
    msg <- sprintf(paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), " | ", "%s"), sprintf(...))
    cat(msg, "\n", file = stderr())
}


.norm_chr <- function(x) {

    # Standardizes chromosome labels by converting them to uppercase and ensuring they manifest_left_bound_key with 'chr'.

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

    # Log SNP match percentage
    pct_match <- 100 * length(present) / length(snps)
    if (pct_match < 80) {
        .log("    WARNING align_ld: only %.1f%% SNP match (present=%d, input=%d)", 
            pct_match, length(present), length(snps))
    }

    # If there are fewer than 2 SNPs present, return NULL
    if (length(present) < 2) return(NULL)

    # Extract the LD matrix and allele information for the present SNPs
    M  <- ld$M[present, present, drop = FALSE]
    a1 <- toupper(trimws(ld$A1[present]))  # ✓ Standardize
    a2 <- toupper(trimws(ld$A2[present]))  # ✓ Standardize
    tgt <- toupper(trimws(target_allele[present]))  # ✓ Standardize


    # Align the LD sign to the effect allele
    # If target allele is the same as the LD panels counted allele (a2), no flip is needed
    # If target allele is the same as the LD panels non-counted allele (a1), flip is needed
    s <- rep(NA_real_, length(present))
    s[a1 == tgt] <- -1
    s[a2 == tgt] <- 1

    # Log allele mismatches. If more than 10% of SNPs have mismatched alleles, return NULL to force fallback to coloc.abf
    unmatched <- sum(is.na(s))
    if (unmatched > 0) {
        pct_unmatched <- 100 * unmatched / length(present)
        .log("    WARNING align_ld: %.1f%% allele mismatch (%d/%d SNPs have target not in A1/A2)", 
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


build_dataset = function(tbl, ld, type, N, s = NULL, sdY = NULL, target_allele_key, standardized_keys, susie_min_snps = 10) {
    
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
        snp = .data[[standardized_keys$variant_id]],
        beta = .data[[standardized_keys$beta]],
        varbeta = .data[[standardized_keys$varbeta]],
        MAF = .data[[standardized_keys$maf]],
        EA = toupper(trimws(.data[[standardized_keys$effect_allele]])),
        OA = toupper(trimws(.data[[standardized_keys$non_effect_allele]])),
        POS = .data[[standardized_keys$pos]]
    )

    # Collapse duplicate SNPs by keeping the one with the smallest p-value (largest absolute z-score)
    D0 <- .collapse_keys(D0)
    D0 <- D0 %>% filter(!.is_ambiguous(EA, OA))

    # Remove SNPs with zero or unreliable variance/SE
    D0 <- D0 %>% 
        filter(is.finite(beta), is.finite(varbeta), varbeta > 0, !is.na(beta), beta != 0)
    
    # Check if there are fewer than 2 SNPs in the dataset and return NULL if so
    if (nrow(D0) < 2) return(NULL)

    # Align the LD matrix to the effect allele and subset it to the SNPs in D0
    if (target_allele_key == standardized_keys$effect_allele) {
        LD <- align_ld(ld, D0$snp, setNames(D0$EA, D0$snp))
    }
    else if (target_allele_key == standardized_keys$non_effect_allele) {
        LD <- align_ld(ld, D0$snp, setNames(D0$OA, D0$snp))
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
    if (nrow(D0) < susie_min_snps) {
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
        LD = LD, 
        MAF = D0$MAF,
        z = D0$beta / sqrt(D0$varbeta)
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
}

#------------- Usage -----------------

# Define standardized keys for column names
standardized_keys <- list(
    chr = "CHR",
    pos = "BP",
    left_bound = "LEFT_500KB",
    right_bound = "RIGHT_500KB",
    gene_id = "GENE",
    variant_id = "VAR",
    beta = "BETA",
    varbeta = "VARBETA",
    maf = "MAF",
    effect_allele = "EFFECT",
    non_effect_allele = "NON_EFFECT"
)

# Define locus and gene of interest
locus <- list(
    ID = "1:52194738-53458606",
    CHR = "1",
    LEFT_500KB = 52144738,
    RIGHT_500KB = 53458606
)
gene_id <- "RNF11"

# Define file paths for QTL, LD, and BIM data
qtl_path <- "/mnt/output/output/coloc/amp_ad/cortex_meta/1kg_v3_hg38_EUR/Cortex_MetaAnalysis_ROSMAP_CMC_HBCC_Mayo_cis_eQTL_release.standardized.tsv"
ld_path <- "/mnt/output/output/coloc/ld_1kg_v3_hg38_EUR/1_52194738_53458606.gwas.ld"
bim_path <- "/mnt/output/output/coloc/ld_1kg_v3_hg38_EUR/1_52194738_53458606.gwas.bim"

# Load QTL data and filter for the locus and gene of interest
qtl_data <- .load_table(qtl_path)
locus_qtl <- qtl_data %>% filter(
    .data[[standardized_keys$chr]] == locus[[standardized_keys$chr]],
    .data[[standardized_keys$pos]] >= locus[[standardized_keys$left_bound]],
    .data[[standardized_keys$pos]] <= locus[[standardized_keys$right_bound]]
)
gene_qtl  <- locus_qtl %>% filter(.data[[standardized_keys$gene_id]] == gene_id)

# Load LD matrix and BIM file
ld <- read_ld(
    ld_path  = file.path(ld_path),
    bim_path = file.path(bim_path),
    panel_note = "1kg_v3_hg38_EUR"
)

# Build the dataset for the gene QTL analysis
gene_qtl_dataset <- build_dataset(tbl = gene_qtl, ld = ld, type = "quant", N = 1694, sdY = 1.0, standardized_keys = standardized_keys, target_allele_key = standardized_keys$effect_allele)

# Calculate z-scores, estimate standard errors, and perform kriging for the gene QTL dataset
z_gene_qtl = gene_qtl_dataset$D$beta / sqrt(gene_qtl_dataset$D$varbeta)
gene_qtl_s <- estimate_s_rss(z_gene_qtl, gene_qtl_dataset$D$LD, gene_qtl_dataset$n)
gene_qtl_kriging <- kriging_rss(z_gene_qtl, gene_qtl_dataset$D$LD, gene_qtl_dataset$n, 1e-08, gene_qtl_s)

# Identify flagged SNPs based on kriging results for the gene QTL dataset
flag <- abs(gene_qtl_kriging$conditional_dist$z_std_diff) > 3 & gene_qtl_kriging$conditional_dist$logLR > 2
flagged_snps_gene_qtl <- gene_qtl_dataset$D$snp[flag]

# Run SuSiE with different configurations for the gene QTL dataset
susie_rss(
    z = z_gene_qtl, 
    R = gene_qtl_dataset$D$LD, 
    n = gene_qtl_dataset$n, 
    coverage = 0.95, 
    max_iter = 1000, 
    L = 5,
    repeat_until_convergence = FALSE
)

susie_rss(
    z = z_gene_qtl, 
    R = gene_qtl_dataset$D$LD, 
    n = gene_qtl_dataset$n, 
    coverage = 0.95, 
    max_iter = 1000, 
    L = 10,
    estimate_prior_variance = TRUE,  # Enable prior variance estimation for regularization
    estimate_prior_method = "EM",     # Use EM algorithm for stable estimation
    scaled_prior_variance = 0.01,     # Regularization strength (smaller = more regularization)
    repeat_until_convergence = FALSE
)

runsusie(
    gene_qtl_dataset$D, coverage = 0.95, max_iter = 1000, L = 5,
    repeat_until_convergence = FALSE
)