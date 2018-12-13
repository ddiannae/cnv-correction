setwd("/mnt/ddisk/transpipeline-data/breast-data-cnvs")

cnvs.caseids <- read.csv("cancer-caseid-aliquotes.tsv", sep="\t", stringsAsFactors = F, header = F, col.names = c("CaseID", "SubmitterID", "AliquotID"))
cnvs <- read.csv("cancer-gistic-results/all_data_by_genes.txt", sep = "\t",
                 header = T, stringsAsFactors = F, check.names =F)
rownames(cnvs.caseids) <- cnvs.caseids$AliquotID
colnames(cnvs) <- cnvs.caseids[colnames(cnvs), "SubmitterID"]
cnvs <- cnvs[, -c(2,3)]
colnames(cnvs)[1] <- "Gene.Symbol"
which(duplicated(cnvs$Gene.Symbol))

write.table(cnvs, file = "cnv-matrix.tsv", sep="\t", quote=FALSE, row.names=FALSE)
