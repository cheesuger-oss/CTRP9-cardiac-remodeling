suppressMessages({library(WGCNA)})
options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly=TRUE)
root <- if (length(args)) normalizePath(args[[1]], mustWork=TRUE) else normalizePath(".")
setwd(root)
allowWGCNAThreads()

run_wgcna <- function(mat_file, meta_file, out_prefix, n_genes = 8000) {
  datExpr <- read.table(mat_file, sep="\t", header=TRUE, row.names=1, check.names=FALSE)
  datExpr <- t(datExpr)  # samples x genes
  datExpr0 <- datExpr
  gsg <- goodSamplesGenes(datExpr0, verbose=0)
  if (!gsg$allOK) datExpr0 <- datExpr0[, gsg$goodGenes]
  # filter to most variable genes
  v <- apply(datExpr0, 2, var, na.rm=TRUE)
  keep <- names(sort(v, decreasing=TRUE))[1:min(n_genes, ncol(datExpr0))]
  datExpr0 <- datExpr0[, keep]
  meta <- read.delim(meta_file, sep="\t")
  meta <- meta[meta$sample %in% rownames(datExpr0),]
  meta <- meta[match(rownames(datExpr0), meta$sample),]
  
  # soft threshold
  powers <- c(seq(1,10), seq(12,30,2))
  sft <- pickSoftThreshold(datExpr0, powerVector=powers, verbose=0)
  fit <- -sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq
  power <- sft$fitIndices$power[which.max(fit >= 0.8)]
  if (length(power)==0 || is.na(power)) power <- 6
  cat(out_prefix, ": chosen power =", power, " scaleFree R2=", max(fit), "\n")
  
  net <- blockwiseModules(datExpr0, power=power, TOMType="unsigned",
                          minModuleSize=30, reassignThreshold=0, mergeCutHeight=0.25,
                          numericLabels=TRUE, pamRespectsDendro=FALSE, verbose=0)
  mods <- net$colors
  names(mods) <- colnames(datExpr0)
  mod_labels <- setNames(paste0("ME", sort(unique(mods))), sort(unique(mods)))
  
  # traits: group indicator + CTRP9 + key markers
  trait_cols <- c("group01")
  datTraits <- data.frame(group01 = ifelse(meta$group %in% c("pab","MCT","TOF","PAB","treat"), 1, 0))
  rownames(datTraits) <- rownames(datExpr0)
  
  MEs <- net$MEs
  colnames(MEs) <- paste0("ME", sub("ME", "", colnames(MEs)))
  nSamples <- nrow(datExpr0)
  moduleTraitCor <- cor(MEs, datTraits, use="p")
  moduleTraitP <- corPvalueStudent(moduleTraitCor, nSamples)
  
  write.table(mod_labels, paste0("results/", out_prefix, "_module_labels.tsv"), sep="\t", quote=FALSE, col.names=FALSE)
  write.table(moduleTraitCor, paste0("results/", out_prefix, "_module_trait_cor.tsv"), sep="\t", quote=FALSE)
  write.table(moduleTraitP, paste0("results/", out_prefix, "_module_trait_p.tsv"), sep="\t", quote=FALSE)
  write.table(data.frame(gene=names(mods), module=mods), paste0("results/", out_prefix, "_gene_modules.tsv"),
              sep="\t", row.names=FALSE, quote=FALSE)
  
  # CTRP9 module
  ctrp9_gene <- if (grepl("rat", out_prefix)) "C1qtnf9" else if (grepl("human", out_prefix)) "C1QTNF9" else "C1QTNF9"
  if (ctrp9_gene %in% names(mods)) {
    cat(out_prefix, ": CTRP9 in module", mods[ctrp9_gene], "\n")
    m <- mods[ctrp9_gene]
    mod_genes <- names(mods)[mods == m]
    write.table(data.frame(gene=mod_genes), paste0("results/", out_prefix, "_ctrp9_module_genes.txt"),
                sep="\t", row.names=FALSE, quote=FALSE, col.names=FALSE)
    cat(out_prefix, ": module size =", length(mod_genes), "\n")
  }
  cat(out_prefix, ": modules =", paste(sort(unique(mods)), collapse=","), "\n")
  list(power=power, mods=mods, MEs=MEs, moduleTraitCor=moduleTraitCor, moduleTraitP=moduleTraitP)
}

run_wgcna("data/processed/GSE36761_log2_matrix.tsv", "data/processed/GSE36761_metadata.tsv", "wgcna_GSE36761_human")
cat("WGCNA_GSE36761_DONE\n")
