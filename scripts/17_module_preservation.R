#!/usr/bin/env Rscript

suppressMessages(library(WGCNA))
options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) normalizePath(args[[1]], mustWork = TRUE) else getwd()
setwd(root)

seed <- 20261115L
n_perm <- 200L

read_expr <- function(path) {
  x <- read.delim(path, row.names = 1, check.names = FALSE)
  x <- t(as.matrix(x))
  storage.mode(x) <- "double"
  as.data.frame(x, check.names = FALSE)
}

read_labels <- function(path) {
  x <- read.delim(path, check.names = FALSE)
  stopifnot(all(c("gene", "module") %in% colnames(x)))
  if (anyDuplicated(x$gene)) stop(paste("Duplicated genes in", path))
  setNames(as.integer(x$module), x$gene)
}

lv <- read_expr("analysis/GSE57338_gene_matrix.tsv")
rv_all <- read_expr("data/processed/GSE36761_log2_matrix.tsv")
provenance <- read.delim("analysis/GSE36761_sample_provenance_corrected.tsv")
rv_samples <- provenance$sample[provenance$chamber == "RV"]
if (!all(rv_samples %in% rownames(rv_all))) stop("RV sample missing from expression matrix")
rv <- rv_all[rv_samples, , drop = FALSE]

lv_labels <- read_labels("analysis/wgcna_GSE57338_gene_modules.tsv")
rv_labels <- read_labels("results/wgcna_GSE36761_human_gene_modules.tsv")
common <- Reduce(intersect, list(colnames(lv), colnames(rv), names(lv_labels), names(rv_labels)))
common <- sort(common)
if (length(common) < 1000L) stop("Unexpectedly small common gene universe")

lv <- lv[, common, drop = FALSE]
rv <- rv[, common, drop = FALSE]
lv_labels <- lv_labels[common]
rv_labels <- rv_labels[common]

gsg_lv <- goodSamplesGenes(lv, verbose = 0)
gsg_rv <- goodSamplesGenes(rv, verbose = 0)
good_genes <- gsg_lv$goodGenes & gsg_rv$goodGenes
if (!all(good_genes)) {
  lv <- lv[, good_genes, drop = FALSE]
  rv <- rv[, good_genes, drop = FALSE]
  lv_labels <- lv_labels[good_genes]
  rv_labels <- rv_labels[good_genes]
  common <- common[good_genes]
}

multi_data <- list(GSE57338_LV = list(data = lv), GSE36761_RV = list(data = rv))

run_direction <- function(reference_index, colors, reference_name, test_name, target_module) {
  set.seed(seed)
  color_list <- list(colors)
  names(color_list) <- reference_name
  result <- modulePreservation(
    multiData = multi_data,
    multiColor = color_list,
    referenceNetworks = reference_index,
    nPermutations = n_perm,
    randomSeed = seed,
    networkType = "unsigned",
    quickCor = 0,
    verbose = 2
  )
  z_container <- result$preservation$Z[[paste0("ref.", reference_name)]][[paste0("inColumnsAlsoPresentIn.", test_name)]]
  rank_container <- result$preservation$observed[[paste0("ref.", reference_name)]][[paste0("inColumnsAlsoPresentIn.", test_name)]]
  z_table <- as.data.frame(z_container)
  rank_table <- as.data.frame(rank_container)
  z_table$module <- rownames(z_table)
  rank_table$module <- rownames(rank_table)
  merged <- merge(z_table, rank_table, by = "module", suffixes = c("_Z", "_observed"), all = TRUE)
  merged$reference_dataset <- reference_name
  merged$test_dataset <- test_name
  merged$prespecified_target_module <- as.integer(merged$module == as.character(target_module))
  merged
}

lv_to_rv <- run_direction(1L, lv_labels, "GSE57338_LV", "GSE36761_RV", 3L)
rv_to_lv <- run_direction(2L, rv_labels, "GSE36761_RV", "GSE57338_LV", 4L)
full <- rbind(lv_to_rv, rv_to_lv)

write.table(
  full,
  "analysis/wgcna_cross_dataset_module_preservation_full.tsv",
  sep = "\t", row.names = FALSE, quote = FALSE, na = "NA"
)

targets <- full[full$prespecified_target_module == 1L, , drop = FALSE]
write.table(
  targets,
  "analysis/wgcna_cross_dataset_module_preservation_targets.tsv",
  sep = "\t", row.names = FALSE, quote = FALSE, na = "NA"
)

manifest <- data.frame(
  field = c("random_seed", "permutations", "common_gene_n", "GSE57338_sample_n",
            "GSE36761_RV_sample_n", "GSE57338_module3_common_n", "GSE36761_module4_common_n",
            "network_type", "WGCNA_version"),
  value = c(seed, n_perm, length(common), nrow(lv), nrow(rv), sum(lv_labels == 3L),
            sum(rv_labels == 4L), "unsigned", as.character(packageVersion("WGCNA")))
)
write.table(
  manifest,
  "analysis/wgcna_cross_dataset_module_preservation_manifest.tsv",
  sep = "\t", row.names = FALSE, quote = FALSE
)

print(targets)
cat("MODULE_PRESERVATION_DONE\n")
