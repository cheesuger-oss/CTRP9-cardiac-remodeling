#!/usr/bin/env Rscript

# Primary count-data analysis for GSE119754.
#
# Frozen inputs:
#   data/processed/GSE119754_counts_symbol.tsv
#   data/processed/GSE119754_metadata.tsv
#
# Outputs (outside data/):
#   analysis/GSE119754_TMM_voom_primary.tsv
#   analysis/GSE119754_TMM_voom_primary_summary.tsv
#   analysis/GSE119754_TMM_voom_library_QC.tsv
#   analysis/GSE119754_TMM_voom_method.txt
#
# Method: retain genes with CPM >= 1 in at least five samples (the smaller
# biological group), calculate Robinson-Oshlack TMM factors, use the effective
# library sizes in limma::voom, fit the MCT-control contrast, apply robust
# empirical Bayes moderation, and control genome-wide FDR by Benjamini-Hochberg.

suppressPackageStartupMessages(library(limma))

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) normalizePath(args[[1]]) else normalizePath(".")
out_dir <- if (length(args) >= 2) args[[2]] else file.path(root, "analysis")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

counts_path <- file.path(root, "data/processed/GSE119754_counts_symbol.tsv")
metadata_path <- file.path(root, "data/processed/GSE119754_metadata.tsv")

counts <- read.delim(counts_path, row.names = 1, check.names = FALSE)
metadata <- read.delim(metadata_path, check.names = FALSE)
stopifnot(all(metadata$sample %in% colnames(counts)))
counts <- as.matrix(counts[, metadata$sample, drop = FALSE])
storage.mode(counts) <- "double"

group <- factor(metadata$group, levels = c("control", "MCT"))
stopifnot(all(table(group) == c(control = 6, MCT = 5)))
design <- model.matrix(~ group)
colnames(design) <- c("Intercept", "MCT_vs_control")

# Filtering was fixed before inspecting the target-gene result.
min_cpm <- 1
min_samples <- min(table(group))
raw_library <- colSums(counts)
raw_cpm <- t(t(counts) / raw_library) * 1e6
keep <- rowSums(raw_cpm >= min_cpm) >= min_samples
filtered <- counts[keep, , drop = FALSE]
if (!"C1qtnf9" %in% rownames(filtered)) stop("C1qtnf9 failed the prespecified filter")

# Robinson-Oshlack TMM, implemented locally because edgeR is unavailable in
# the frozen project runtime. Trimming is rank based, matching the published
# TMM definition; factors are centered to geometric mean 1.
tmm_factor_one <- function(obs, ref, lib_obs, lib_ref,
                           logratio_trim = 0.30, sum_trim = 0.05,
                           acutoff = -1e10) {
  keep_pos <- (obs > 0) & (ref > 0)
  obs <- obs[keep_pos]
  ref <- ref[keep_pos]
  if (length(obs) == 0) return(1)
  log_r <- log2((obs / lib_obs) / (ref / lib_ref))
  abs_e <- 0.5 * (log2(obs / lib_obs) + log2(ref / lib_ref))
  variance <- (lib_obs - obs) / lib_obs / obs +
    (lib_ref - ref) / lib_ref / ref
  finite <- is.finite(log_r) & is.finite(abs_e) & is.finite(variance) &
    variance > 0 & abs_e > acutoff
  log_r <- log_r[finite]
  abs_e <- abs_e[finite]
  variance <- variance[finite]
  n <- length(log_r)
  if (n == 0) return(1)
  rank_m <- rank(log_r, ties.method = "average")
  rank_a <- rank(abs_e, ties.method = "average")
  lo_m <- floor(n * logratio_trim) + 1
  hi_m <- n + 1 - lo_m
  lo_a <- floor(n * sum_trim) + 1
  hi_a <- n + 1 - lo_a
  use <- rank_m >= lo_m & rank_m <= hi_m & rank_a >= lo_a & rank_a <= hi_a
  if (!any(use)) return(1)
  2^(sum(log_r[use] / variance[use]) / sum(1 / variance[use]))
}

upper_quartile_cpm <- apply(raw_cpm[keep, , drop = FALSE], 2, function(x) {
  quantile(x[x > 0], probs = 0.75, names = FALSE)
})
ref_index <- which.min(abs(upper_quartile_cpm - mean(upper_quartile_cpm)))
ref_sample <- colnames(filtered)[ref_index]
norm_factors <- vapply(seq_len(ncol(filtered)), function(i) {
  tmm_factor_one(filtered[, i], filtered[, ref_index],
                 sum(filtered[, i]), sum(filtered[, ref_index]))
}, numeric(1))
norm_factors <- norm_factors / exp(mean(log(norm_factors)))
names(norm_factors) <- colnames(filtered)

filtered_library <- colSums(filtered)
effective_library <- filtered_library * norm_factors
voom_object <- voom(filtered, design = design, lib.size = effective_library,
                    normalize.method = "none", plot = FALSE)
fit <- lmFit(voom_object, design)
fit <- eBayes(fit, robust = TRUE)
tab <- topTable(fit, coef = "MCT_vs_control", number = Inf,
                sort.by = "none", confint = TRUE, adjust.method = "BH")
tab$gene <- rownames(tab)
tab <- tab[, c("gene", "logFC", "CI.L", "CI.R", "AveExpr", "t",
               "P.Value", "adj.P.Val", "B")]

full_out <- file.path(out_dir, "GSE119754_TMM_voom_primary.tsv")
write.table(tab, full_out, sep = "\t", quote = FALSE, row.names = FALSE)

target <- tab[tab$gene == "C1qtnf9", , drop = FALSE]
summary <- data.frame(
  dataset = "GSE119754",
  contrast = "MCT-control",
  n_control = sum(group == "control"),
  n_case = sum(group == "MCT"),
  genes_input = nrow(counts),
  genes_retained = nrow(filtered),
  min_cpm = min_cpm,
  min_samples = min_samples,
  gene = target$gene,
  log2FC = target$logFC,
  CI_low = target$CI.L,
  CI_high = target$CI.R,
  P_value = target$P.Value,
  BH_FDR = target$adj.P.Val,
  method = "CPM filter + TMM effective libraries + voom-limma robust empirical Bayes"
)
summary_out <- file.path(out_dir, "GSE119754_TMM_voom_primary_summary.tsv")
write.table(summary, summary_out, sep = "\t", quote = FALSE, row.names = FALSE)

qc <- data.frame(
  sample = colnames(filtered),
  group = as.character(group),
  raw_library_size = raw_library[colnames(filtered)],
  filtered_library_size = filtered_library,
  TMM_factor = norm_factors,
  effective_library_size = effective_library,
  reference_sample = colnames(filtered) == ref_sample
)
qc_out <- file.path(out_dir, "GSE119754_TMM_voom_library_QC.tsv")
write.table(qc, qc_out, sep = "\t", quote = FALSE, row.names = FALSE)

method_lines <- c(
  paste0("Frozen counts: ", counts_path),
  paste0("Frozen metadata: ", metadata_path),
  paste0("Filter: CPM >= ", min_cpm, " in >= ", min_samples,
         " samples (smaller biological group)"),
  paste0("Genes retained: ", nrow(filtered), "/", nrow(counts)),
  paste0("TMM reference sample: ", ref_sample),
  "Normalization: Robinson-Oshlack TMM factors centered to geometric mean 1",
  "Model: limma::voom with TMM effective library sizes; robust empirical Bayes",
  "Contrast: MCT minus control",
  "Multiplicity: Benjamini-Hochberg across all retained genes",
  paste0("R: ", R.version.string),
  paste0("limma: ", as.character(packageVersion("limma")))
)
writeLines(method_lines, file.path(out_dir, "GSE119754_TMM_voom_method.txt"))

print(summary)
cat("Outputs:\n", full_out, "\n", summary_out, "\n", qc_out, "\n", sep = "")
