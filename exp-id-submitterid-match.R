setwd("/mnt/ddisk/transpipeline-data/breast-data-cnvs")

load("/mnt/ddisk/transpipeline-data/breast-data/rdata/Length.full_GC.full_Between.tmm_Norm_cpm10_arsyn.RData")

norm.data_cpm10_arsyn$Targets$FileName = paste(substring(norm.data_cpm10_arsyn$Targets$File,20, length(norm.data_cpm10_arsyn$Targets$File)), "gz", sep=".")
exp.caseids = read.csv("exp-caseids-cancer.csv", header = T, col.names = c("FileName", "SubmitterID", "CaseID", "FileID"), stringsAsFactors = F)
targets.merged <- merge(exp.caseids, norm.data_cpm10_arsyn$Targets, by = "FileName" )
rownames(targets.merged) <- targets.merged$ID

tumor <- as.data.frame(norm.data_cpm10_arsyn$M[,norm.data_cpm10_arsyn$Targets$Group == "T"])
colnames(tumor) <- targets.merged[colnames(tumor), "SubmitterID"]

genes = read.csv("/mnt/ddisk/transpipeline-data/biomarts/Biomart_EnsemblG94_GRCh38_p12.txt", header = T, sep="\t",  stringsAsFactors = F)
duplicated.gene.ids <- which(duplicated(genes$Gene.stable.ID))
genes[genes$Gene.stable.ID %in% genes[duplicated.gene.ids, "Gene.stable.ID"], ]
genes <- genes[!duplicated(genes$Gene.stable.ID), ]
rownames(genes) <- genes$Gene.stable.ID
genes <- genes[order(genes$Gene.stable.ID), c("Gene.stable.ID","HGNC.symbol")]

symbols <- genes[, "HGNC.symbol"]
names(symbols) <- genes[, "Gene.stable.ID"]
symbols[532]
symbols[which(symbols == "")] = names(which(symbols == ""))
symbols[532]
tumor.symbols <- cbind(symbols[rownames(tumor)], tumor, stringsAsFactors = F)
colnames(tumor.symbols)[1] <- "Gene.Symbol"
which(duplicated(tumor.symbols$Gene.Symbol))
tumor.symbols[13812, "Gene.Symbol"] <- rownames(tumor.symbols)[13812]
tumor.symbols[14391, "Gene.Symbol"] <- rownames(tumor.symbols)[14391]
which(is.na(tumor.symbols$Gene.Symbol))
tumor.symbols[12290, "Gene.Symbol"] <- rownames(tumor.symbols)[12290]
rownames(tumor.symbols) <- tumor.symbols$Gene.Symbol

write.table(tumor.symbols, file = "exp-matrix.tsv", sep="\t", quote=FALSE, row.names=FALSE)
