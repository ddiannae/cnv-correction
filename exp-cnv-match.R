setwd("/mnt/ddisk/transpipeline-data/breast-data-cnvs")

cnvs <- read.table("cnv-matrix.tsv", header = T, check.names = F, row.names = 1)
exp  <- read.table("exp-matrix.tsv", header = T, check.names = F, row.names = 1)

cases.intersect <- intersect(colnames(cnvs), colnames(exp))
genes.intersect <- intersect(rownames(cnvs), rownames(exp))

cnvs.filtered <- cnvs[genes.intersect, cases.intersect]
exp.filtered <- exp[genes.intersect, cases.intersect]

cnvs.filtered["Gene.Symbol"] <- rownames(cnvs.filtered)
exp.filtered["Gene.Symbol"] <- rownames(exp.filtered)
cnvs.filtered <- cnvs.filtered[, c(ncol(cnvs.filtered), 1:ncol(cnvs.filtered)-1)]
exp.filtered <- exp.filtered[, c(ncol(exp.filtered), 1:ncol(exp.filtered)-1)]
write.table(cnvs.filtered, file = "cnv-filtered-matrix.tsv", sep="\t", quote=FALSE, row.names=FALSE)
write.table(exp.filtered, file = "exp-filtered-matrix.tsv", sep="\t", quote=FALSE, row.names=FALSE)
