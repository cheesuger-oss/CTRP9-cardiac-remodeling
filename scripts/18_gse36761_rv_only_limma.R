suppressMessages(library(limma))

script_arg <- commandArgs(trailingOnly=FALSE)[grep("--file=", commandArgs(trailingOnly=FALSE))]
root <- if (length(script_arg)) normalizePath(file.path(dirname(sub("--file=", "", script_arg)), "..")) else normalizePath(".")
setwd(root)

rv_controls <- c("NH-02", "NH-04", "NH-06", "NH-08")
mat <- read.delim("data/processed/GSE36761_log2_matrix.tsv", row.names=1, check.names=FALSE)
ctrp9 <- read.delim("data/processed/GSE36761_C1QTNF9_transcript_sum.tsv", row.names=1, check.names=FALSE)
mat["C1QTNF9", ] <- log2(as.numeric(ctrp9[1, ]) + 1)
tof <- grep("^TOF-", colnames(mat), value=TRUE)
keep <- c(rv_controls, tof)
mat <- mat[, keep, drop=FALSE]
group <- factor(c(rep("NH_RV", length(rv_controls)), rep("TOF_RV", length(tof))), levels=c("NH_RV", "TOF_RV"))
design <- model.matrix(~group)
fit <- eBayes(lmFit(mat, design))
res <- topTable(fit, coef=2, number=Inf, sort.by="none")
res$gene <- rownames(res)
res$n_NH_RV <- length(rv_controls)
res$n_TOF_RV <- length(tof)
write.table(res, "analysis/GSE36761_RV_only_limma.tsv", sep="\t", row.names=FALSE, quote=FALSE)

cat("RV-only limma:", sum(res$adj.P.Val < 0.05), "DE genes;",
    sum(res$adj.P.Val < 0.05 & res$logFC > 0), "up;",
    sum(res$adj.P.Val < 0.05 & res$logFC < 0), "down\n")
print(res[res$gene %in% c("C1QTNF9", "ACTA1", "POSTN", "FN1", "COL1A1", "COL3A1", "TGFB1"),
          c("gene", "logFC", "P.Value", "adj.P.Val")])
