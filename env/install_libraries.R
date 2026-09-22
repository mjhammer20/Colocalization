.required_cran <- c(
    "data.table",
    "dplyr",
    "readr",
    "stringr",
    "tidyr",
    "purrr",
    "R6",
    "argparse",
    "tools"
)

.required_bioc <- c(
    "susieR",
    "coloc"
)

# Install missing CRAN packages
.missing_cran <- .required_cran[!sapply(.required_cran, requireNamespace, quietly = TRUE)]
if (length(.missing_cran)) {
    message(sprintf("Installing missing CRAN packages: %s", paste(.missing_cran, collapse = ", ")))
    install.packages(.missing_cran, repos = "https://cloud.r-project.org", quiet = TRUE)
}

# Install missing Bioconductor/GitHub packages
.missing_bioc <- .required_bioc[!sapply(.required_bioc, requireNamespace, quietly = TRUE)]
if (length(.missing_bioc)) {
    message(sprintf("Installing missing packages via BiocManager: %s", paste(.missing_bioc, collapse = ", ")))
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", repos = "https://cloud.r-project.org", quiet = TRUE)
    }
    BiocManager::install(.missing_bioc, ask = FALSE, update = FALSE)
}

# Final check — stop if any required package still cannot be loaded
.all_required <- c(.required_cran, .required_bioc)
.still_missing <- .all_required[!sapply(.all_required, requireNamespace, quietly = TRUE)]
if (length(.still_missing)) {
    stop(sprintf(
        "The following required packages could not be installed and are still missing:\n  %s",
        paste(.still_missing, collapse = "\n  ")
    ))
}
