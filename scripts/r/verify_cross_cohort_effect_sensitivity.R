#!/usr/bin/env Rscript

# Independent base-R verification of Hedges' g and leave-one-out correlations.
# All inputs are read-only. The output records numerical agreement with the
# Python implementation used for the primary analysis.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[1]))
root <- normalizePath(file.path(dirname(script_path), "..", ".."))

read_tsv <- function(path, row_names = FALSE) {
  if (row_names) {
    read.delim(path, check.names = FALSE, row.names = 1, stringsAsFactors = FALSE)
  } else {
    read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  }
}

clean_meta <- function(path) {
  d <- read_tsv(path)
  d <- d[, !grepl("^Unnamed", names(d)), drop = FALSE]
  rownames(d) <- d$sample
  d
}

hedges <- function(case, control) {
  n1 <- length(case); n0 <- length(control); df <- n1 + n0 - 2
  sp2 <- ((n1 - 1) * var(case) + (n0 - 1) * var(control)) / df
  d <- (mean(case) - mean(control)) / sqrt(sp2)
  j <- 1 - 3 / (4 * df - 1)
  g <- j * d
  v <- (n1 + n0) / (n1 * n0) + g^2 / (2 * df)
  c(g = g, variance = v, ci_low = g - qnorm(.975) * sqrt(v),
    ci_high = g + qnorm(.975) * sqrt(v))
}

effects <- list()

pab_m <- read_tsv(file.path(root, "data/processed/GSE210180_log2fpkm_matrix.tsv"), TRUE)
pab_meta <- clean_meta(file.path(root, "data/processed/GSE210180_metadata.tsv"))
pab_x <- as.numeric(pab_m["C1qtnf9", rownames(pab_meta)])
effects[["GSE210180"]] <- hedges(pab_x[pab_meta$group == "pab"], pab_x[pab_meta$group == "sham"])

mct_m <- read_tsv(file.path(root, "data/processed/GSE119754_counts_symbol.tsv"), TRUE)
mct_meta <- clean_meta(file.path(root, "data/processed/GSE119754_metadata.tsv"))
mct_log <- log2(sweep(mct_m, 2, colSums(mct_m), "/") * 1e6 + 1)
mct_x <- as.numeric(mct_log["C1qtnf9", rownames(mct_meta)])
effects[["GSE119754"]] <- hedges(mct_x[mct_meta$group == "MCT"], mct_x[mct_meta$group == "control"])

tof_m <- read_tsv(file.path(root, "data/processed/GSE36761_log2_matrix.tsv"), TRUE)
prov <- read_tsv(file.path(root, "analysis/GSE36761_sample_provenance_corrected.tsv"))
prov <- prov[prov$chamber == "RV", ]
tof_x <- as.numeric(tof_m["C1QTNF9", prov$sample])
effects[["GSE36761"]] <- hedges(tof_x[prov$group == "TOF"], tof_x[prov$group == "NH"])

hf_m <- read_tsv(file.path(root, "analysis/GSE57338_gene_matrix.tsv"), TRUE)
con <- gzfile(file.path(root, "data/raw/GSE57338_series_matrix.txt.gz"), "rt")
titles <- NULL; accessions <- NULL
repeat {
  line <- readLines(con, n = 1)
  if (length(line) == 0) break
  if (startsWith(line, "!Sample_title")) {
    titles <- strsplit(gsub('"', '', line), "\t", fixed = FALSE)[[1]][-1]
  }
  if (startsWith(line, "!Sample_geo_accession")) {
    accessions <- strsplit(gsub('"', '', line), "\t", fixed = FALSE)[[1]][-1]
  }
  if (!is.null(titles) && !is.null(accessions)) break
}
close(con)
groups <- ifelse(grepl("Non-failing", titles), "NF",
                 ifelse(grepl("Ischemic", titles), "ICM", "DCM"))
stopifnot(identical(colnames(hf_m), accessions))
hf_x <- as.numeric(hf_m["C1QTNF9", ])
effects[["GSE57338"]] <- hedges(hf_x[groups != "NF"], hf_x[groups == "NF"])

python_effect <- read_tsv(file.path(root, "analysis/cross_cohort_CTRP9_effect_sizes.tsv"))
effect_check <- do.call(rbind, lapply(names(effects), function(dataset) {
  r <- effects[[dataset]]
  p <- python_effect[python_effect$dataset == dataset, ]
  data.frame(
    check = "Hedges g",
    dataset = dataset,
    marker = "C1QTNF9/C1qtnf9",
    max_abs_difference = max(abs(c(
      r["g"] - p$hedges_g,
      r["variance"] - p$sampling_variance,
      r["ci_low"] - p$ci_low,
      r["ci_high"] - p$ci_high
    )))
  )
}))

full <- read_tsv(file.path(root, "analysis/correlation_multiplicity.tsv"))
python_loo <- read_tsv(file.path(root, "analysis/rat_marker_leave_one_out_correlations.tsv"))
matrices <- list(GSE210180 = pab_m[, rownames(pab_meta), drop = FALSE],
                 GSE119754 = mct_log[, rownames(mct_meta), drop = FALSE])
loo_checks <- list()
for (dataset in unique(full$dataset)) {
  mat <- matrices[[dataset]]
  samples <- colnames(mat)
  x <- as.numeric(mat["C1qtnf9", ])
  for (marker in full$marker[full$dataset == dataset]) {
    y <- as.numeric(mat[marker, ])
    r_values <- sapply(seq_along(samples), function(i) cor(x[-i], y[-i], method = "pearson"))
    rho_values <- sapply(seq_along(samples), function(i) cor(x[-i], y[-i], method = "spearman"))
    py <- python_loo[python_loo$dataset == dataset & python_loo$marker == marker, ]
    py <- py[match(samples, py$omitted_sample), ]
    loo_checks[[length(loo_checks) + 1]] <- data.frame(
      check = "leave-one-out correlations",
      dataset = dataset,
      marker = marker,
      max_abs_difference = max(abs(c(r_values - py$pearson_r, rho_values - py$spearman_rho)))
    )
  }
}

result <- rbind(effect_check, do.call(rbind, loo_checks))
result$status <- ifelse(result$max_abs_difference < 1e-10, "PASS", "FAIL")
write.table(result,
            file.path(root, "analysis/cross_cohort_effect_sensitivity_R_check.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("R cross-check max absolute difference:", max(result$max_abs_difference), "\n")
cat("Checks:", nrow(result), "PASS:", sum(result$status == "PASS"),
    "FAIL:", sum(result$status == "FAIL"), "\n")
if (any(result$status == "FAIL")) quit(status = 1)
