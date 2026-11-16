suppressMessages({library(limma)})
args <- commandArgs(trailingOnly=TRUE)
root <- if (length(args)) normalizePath(args[[1]], mustWork=TRUE) else normalizePath(".")
setwd(root)

# ---------- GSE210180 rat PAB (log2 FPKM) ----------
mat <- read.table("data/processed/GSE210180_log2fpkm_matrix.tsv", sep="\t", header=TRUE, row.names=1, check.names=FALSE)
meta <- read.delim("data/processed/GSE210180_metadata.tsv", sep="\t")
grp <- factor(meta$group, levels=c("sham","pab"))
design <- model.matrix(~grp)
fit <- lmFit(mat, design)
fit <- eBayes(fit)
res <- topTable(fit, coef=2, number=Inf, sort.by="none")
res$gene <- rownames(res)
write.table(res, "results/limma_GSE210180_rat.tsv", sep="\t", row.names=FALSE, quote=FALSE)
show_gene <- function(res, gene) {
  sub <- res[res$gene == gene, c("logFC","P.Value","adj.P.Val")]
  if (nrow(sub) > 0) { cat(gene, ":"); print(round(sub, 4)) } else cat(gene, ": absent\n")
}
cat("GSE210180:", sum(res$adj.P.Val < 0.05), "DE (padj<0.05)\n")
show_gene(res, "C1qtnf9"); show_gene(res, "Nppb"); show_gene(res, "Myh7"); show_gene(res, "Col1a1")

# ---------- GSE36761 human TOF (log2 normalized) ----------
mat2 <- read.table("data/processed/GSE36761_log2_matrix.tsv", sep="\t", header=TRUE, row.names=1, check.names=FALSE)
ctrp9 <- read.table("data/processed/GSE36761_C1QTNF9_transcript_sum.tsv", sep="\t", header=TRUE, row.names=1, check.names=FALSE)
mat2["C1QTNF9", ] <- log2(as.numeric(ctrp9[1,]) + 1)
meta2 <- read.delim("data/processed/GSE36761_metadata.tsv", sep="\t")
grp2 <- factor(ifelse(grepl("^NH", meta2$sample), "NH", "TOF"), levels=c("NH","TOF"))
design2 <- model.matrix(~grp2)
fit2 <- lmFit(mat2, design2)
fit2 <- eBayes(fit2)
res2 <- topTable(fit2, coef=2, number=Inf, sort.by="none")
res2$gene <- rownames(res2)
write.table(res2, "results/limma_GSE36761_human.tsv", sep="\t", row.names=FALSE, quote=FALSE)
cat("GSE36761:", sum(res2$adj.P.Val < 0.05), "DE (padj<0.05)\n")
show_gene(res2, "C1QTNF9"); show_gene(res2, "NPPB"); show_gene(res2, "COL1A1"); show_gene(res2, "POSTN"); show_gene(res2, "FN1"); show_gene(res2, "TGFB1")

# ---------- GSE119754 rat MCT (logCPM) ----------
cnt <- read.table("data/processed/GSE119754_counts_symbol.tsv", sep="\t", header=TRUE, row.names=1, check.names=FALSE)
libsize <- colSums(cnt)
cpm <- t(t(cnt)/libsize)*1e6
logcpm <- log2(cpm+1)
meta3 <- read.delim("data/processed/GSE119754_metadata.tsv", sep="\t")
grp3 <- factor(meta3$group, levels=c("control","MCT"))
design3 <- model.matrix(~grp3)
fit3 <- lmFit(logcpm, design3)
fit3 <- eBayes(fit3)
res3 <- topTable(fit3, coef=2, number=Inf, sort.by="none")
res3$gene <- rownames(res3)
write.table(res3, "results/limma_GSE119754_rat.tsv", sep="\t", row.names=FALSE, quote=FALSE)
cat("GSE119754:", sum(res3$adj.P.Val < 0.05), "DE (padj<0.05)\n")
show_gene(res3, "C1qtnf9"); show_gene(res3, "Postn"); show_gene(res3, "Col1a1"); show_gene(res3, "Nppb")

# ---------- GSE224428 pig PAB (logCPM) ----------
pig <- read.table("data/processed/GSE224428_counts_symbol.tsv", sep="\t", header=TRUE, row.names=1, check.names=FALSE)
meta4 <- read.delim("data/processed/GSE224428_metadata.tsv", sep="\t")
meta4 <- meta4[meta4$sample %in% colnames(pig),]
pig <- pig[, meta4$sample]
libsize4 <- colSums(pig)
cpm4 <- t(t(pig)/libsize4)*1e6
logcpm4 <- log2(cpm4+1)
grp4 <- factor(meta4$group, levels=c("control","PAB"))
design4 <- model.matrix(~grp4)
fit4 <- lmFit(logcpm4, design4)
fit4 <- eBayes(fit4)
res4 <- topTable(fit4, coef=2, number=Inf, sort.by="none")
res4$gene <- rownames(res4)
write.table(res4, "results/limma_GSE224428_pig.tsv", sep="\t", row.names=FALSE, quote=FALSE)
cat("GSE224428:", sum(res4$adj.P.Val < 0.05), "DE (padj<0.05)\n")
for (g in c("NPPB","POSTN","COL1A1","FN1")) show_gene(res4, g)
cat("LIMMA_DONE\n")
