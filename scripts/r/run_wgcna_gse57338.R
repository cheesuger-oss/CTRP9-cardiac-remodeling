suppressMessages({library(WGCNA)})
options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly=TRUE)
root <- if (length(args)) normalizePath(args[[1]], mustWork=TRUE) else normalizePath(".")
setwd(root)
allowWGCNAThreads()

datExpr <- read.table("analysis/GSE57338_gene_matrix.tsv", sep="\t", header=TRUE, row.names=1, check.names=FALSE)
datExpr <- t(datExpr)  # samples x genes
gsg <- goodSamplesGenes(datExpr, verbose=0)
if (!gsg$allOK) datExpr <- datExpr[, gsg$goodGenes]
v <- apply(datExpr, 2, var, na.rm=TRUE)
keep <- names(sort(v, decreasing=TRUE))[1:min(8000, ncol(datExpr))]
datExpr0 <- datExpr[, keep]

powers <- c(seq(1,10), seq(12,30,2))
sft <- pickSoftThreshold(datExpr0, powerVector=powers, verbose=0)
fit <- -sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq
power <- sft$fitIndices$power[which.max(fit >= 0.8)]
if (length(power)==0 || is.na(power)) power <- 6
cat("chosen power =", power, "scaleFree R2 =", max(fit), "\n")

net <- blockwiseModules(datExpr0, power=power, TOMType="unsigned",
                        minModuleSize=40, reassignThreshold=0, mergeCutHeight=0.25,
                        numericLabels=TRUE, pamRespectsDendro=FALSE, verbose=0)
mods <- net$colors
names(mods) <- colnames(datExpr0)
write.table(data.frame(gene=names(mods), module=mods),
            "analysis/wgcna_GSE57338_gene_modules.tsv", sep="\t", row.names=FALSE, quote=FALSE)
mod_sizes <- table(mods)
cat("module sizes:", paste(names(mod_sizes), mod_sizes, sep="=", collapse=", "), "\n")
if ("C1QTNF9" %in% names(mods)) {
  m <- mods["C1QTNF9"]
  cat("C1QTNF9 in module", m, "\n")
  mod_genes <- names(mods)[mods == m]
  write.table(data.frame(gene=mod_genes), "analysis/wgcna_GSE57338_C1QTNF9_module_genes.txt",
              sep="\t", row.names=FALSE, quote=FALSE, col.names=FALSE)
  cat("C1QTNF9 module size =", length(mod_genes), "\n")
} else {
  cat("C1QTNF9 not in kept genes\n")
}
cat("WGCNA_GSE57338_DONE\n")
