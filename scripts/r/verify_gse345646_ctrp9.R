#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) normalizePath(args[[1]], mustWork = TRUE) else getwd()
setwd(root)

input <- "analysis/external_validation/GSE345646/results/C1QTNF9_pseudobulk_by_donor_celltype.tsv"
output <- "analysis/external_validation/GSE345646/results/C1QTNF9_stage_tests_R_crosscheck.tsv"
d <- read.delim(input, check.names = FALSE)
order <- c("NF", "pRV", "RVF")

labelings <- list()
idx <- seq_len(11)
counter <- 1L
for (nf in combn(idx, 4, simplify = FALSE)) {
  remain <- setdiff(idx, nf)
  for (prv in combn(remain, 4, simplify = FALSE)) {
    labels <- rep(2L, 11)
    labels[nf] <- 0L
    labels[prv] <- 1L
    labelings[[counter]] <- labels
    counter <- counter + 1L
  }
}

rows <- list()
for (cell in sort(unique(d$cell_type))) {
  x <- d[d$cell_type == cell & d$informative_nuclei_ge_20 == "True", ]
  group_n <- table(factor(x$group, levels = order))
  eligible <- sum(x$C1QTNF9_raw_count) >= 20 && all(group_n >= 3)
  row <- data.frame(cell_type = cell, eligible = eligible,
                    kruskal_H = NA_real_, kruskal_exact_p = NA_real_,
                    spearman_rho = NA_real_, spearman_exact_p = NA_real_)
  if (eligible) {
    x$group <- factor(x$group, levels = order, ordered = TRUE)
    x <- x[order(x$group, x$donor), ]
    values <- x$C1QTNF9_log2_CPM_plus_0.5
    labels <- as.integer(x$group) - 1L
    h_obs <- as.numeric(kruskal.test(values ~ factor(labels))$statistic)
    rho_obs <- cor(values, labels, method = "spearman")
    h_null <- vapply(labelings, function(z) as.numeric(kruskal.test(values ~ factor(z))$statistic), numeric(1))
    rho_null <- vapply(labelings, function(z) cor(values, z, method = "spearman"), numeric(1))
    row$kruskal_H <- h_obs
    row$kruskal_exact_p <- mean(h_null >= h_obs - 1e-12)
    row$spearman_rho <- rho_obs
    row$spearman_exact_p <- mean(abs(rho_null) >= abs(rho_obs) - 1e-12)
  }
  rows[[length(rows) + 1L]] <- row
}
out <- do.call(rbind, rows)
out$kruskal_BH_FDR <- NA_real_
out$spearman_BH_FDR <- NA_real_
keep <- which(out$eligible)
out$kruskal_BH_FDR[keep] <- p.adjust(out$kruskal_exact_p[keep], method = "BH")
out$spearman_BH_FDR[keep] <- p.adjust(out$spearman_exact_p[keep], method = "BH")
write.table(out, output, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
print(out)
cat("GSE345646_R_CROSSCHECK_DONE\n")
