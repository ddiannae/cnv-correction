library(ComplexHeatmap)
library(circlize)
setwd("/mnt/ddisk/transpipeline-data/breast-data-cnvs")

exp <- read.delim("exp-filtered-matrix.tsv", header = T, stringsAsFactors = F, row.names = 1, check.names = F)
cnv <- read.delim("cnv-filtered-matrix.tsv", header = T, stringsAsFactors = F, row.names = 1, check.names = F)
genes.annot <- read.delim("/mnt/ddisk/transpipeline-data/biomarts/Biomart_EnsemblG94_GRCh38_p12.txt", header = T, stringsAsFactors = F)
colnames(genes.annot) <- c("ID", "Chr", "start", "end", "GC", "type", "symbol")
chrs <- as.character(seq(1, 22))
  
lapply(chrs, function(chr) {
      
    genes <- genes.annot[genes.annot$Chr == chr, c("symbol", "start")]
    genes <- genes[!duplicated(genes$symbol), ]
    exp_df <- merge(genes, exp, by.x = "symbol", by.y = "row.names")
    exp_df <-exp_df[order(exp_df$start), ]
    rownames(exp_df) <- exp_df$symbol
    exp_matrix <- as.matrix(exp_df[, -c(1,2)])
    write.table(exp_matrix, file = paste("chr-files/chr", chr, "exp.tsv", sep="_"), quote = F, row.names = T, col.names = T)
    
    cnvs_df <- merge(genes, cnv, by.x = "symbol", by.y = "row.names")
    cnvs_df <-cnvs_df[order(cnvs_df$start), ]
    rownames(cnvs_df) <- cnvs_df$symbol
    cnvs_matrix <- as.matrix(cnvs_df[, -c(1,2)])
    write.table(cnvs_matrix, file = paste("chr-files/chr", chr, "cnv.tsv", sep="_"), quote = F, row.names = T, col.names = T)
    
    exp_pearson <- cor(t(exp_matrix), method = "pearson")
    exp_pearson_pos <- apply(exp_pearson, c(1,2), function(x) {max(0, x)})
    exp_pearson_neg <- apply(exp_pearson, c(1,2), function(x) {min(0, x)})
    
    min_neg <- if(min(exp_pearson_neg) == 0) -1 else min(exp_pearson_neg)
    
    h1 <- Heatmap(exp_pearson_pos, col = colorRamp2(c(min_neg, 0, 1), c("red", "white", "blue")),  column_title = "Positive",  column_title_gp = gpar(fontsize = 24, fontface = "bold"),
              cluster_rows = FALSE, cluster_columns = FALSE, show_heatmap_legend = F, show_column_names = F, show_row_names = F)
    
    h2 <- Heatmap(exp_pearson_neg, col = colorRamp2(c(min_neg, 0, 1), c("red", "white", "blue")),  column_title = paste("Negative - min val:", min(exp_pearson_neg)),  column_title_gp = gpar(fontsize = 24, fontface = "bold"),
            cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = F, show_row_names = F)

    ht_list = h1 + h2
    png(paste("chr-files/chr", chr, "exp.png", sep="_"), width = 2400, height = 1200)
    draw(ht_list, column_title = paste( "Chromosome", chr, "exp", sep="_"),  column_title_gp = gpar(fontsize = 24, fontface = "bold"))
    dev.off()
    
    cnvs_pearson <- cor(t(cnvs_matrix), method = "pearson")
    cnvs_pearson_pos <- apply(cnvs_pearson, c(1,2), function(x) {max(0, x)})
    cnvs_pearson_neg <- apply(cnvs_pearson, c(1,2), function(x) {min(0, x)})
    
    min_neg <- if(min(cnvs_pearson_neg) == 0) -1 else min(cnvs_pearson_neg)
    
    h1 <- Heatmap(cnvs_pearson_pos, col = colorRamp2(c(min_neg, 0, 1), c("red", "white", "blue")), column_title = "Positive",  column_title_gp = gpar(fontsize = 24, fontface = "bold"),
                  cluster_rows = FALSE, cluster_columns = FALSE, show_heatmap_legend = F, show_column_names = F, show_row_names = F)
    
    h2 <- Heatmap(cnvs_pearson_neg, col = colorRamp2(c(min_neg, 0, 1), c("red", "white", "blue")), column_title = paste("Negative - min val:", min(cnvs_pearson_neg)),  column_title_gp = gpar(fontsize = 24, fontface = "bold"),
                  cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = F, show_row_names = F)
    
    ht_list = h1 + h2
    png(paste( "chr-files/chr", chr, "cnv.png", sep="_"), width = 2400, height = 1200)
    draw(ht_list, column_title = paste( "Chromosome", chr, "cnv", sep="_"),  column_title_gp = gpar(fontsize = 24, fontface = "bold"))
    dev.off()
    
    cols_sorted <- sort(colnames(cnvs_matrix))
    cnvs_matrix <- cnvs_matrix[, cols_sorted]
    exp_matrix <- exp_matrix[, cols_sorted]
    sum(rownames(exp_matrix) == rownames(cnvs_matrix)) == nrow(exp_matrix)
    
    gnames <- rownames(exp_matrix)
    residuals <-parallel::mclapply(X = gnames, mc.cores = 4,  mc.cleanup = FALSE, FUN= function(i){
      recta <- lm(exp_matrix[i, ] ~ cnvs_matrix[i, ])
      a = coefficients(recta)["(Intercept)"]
      b = coefficients(recta)["cnvs_matrix[i, ]"]
      c(i, a, b)
    })
    
    df <- plyr::ldply(residuals)
    colnames(df) <- c("symbol", "b0", "b1")
    df$b0 <- as.numeric(df$b0)
    df$b1 <- as.numeric(df$b1)
    rownames(df) <- df$symbol
    
    exp_corr <-parallel::mclapply(X = gnames, mc.cores = 4,  mc.cleanup = FALSE, FUN= function(i){
      c(i, exp_matrix[i, ]  - df[i, "b0"] - df[i, "b1"]* cnvs_matrix[i,])
    })
    exp_corr <- plyr::ldply(exp_corr)
    exp_matrix_corr <- data.matrix(exp_corr[, 2:ncol(exp_corr)])
    rownames(exp_matrix_corr) <- exp_corr[, 1]
    
    write.table(exp_matrix_corr, file = paste( "chr-files/chr", chr, "exp_corrected.tsv", sep="_"), quote = F, row.names = T, col.names = T)
    
    exp_pearson <- cor(t(exp_matrix_corr), method = "pearson")
    exp_pearson_pos <- apply(exp_pearson, c(1,2), function(x) {max(0, x)})
    exp_pearson_neg <- apply(exp_pearson, c(1,2), function(x) {min(0, x)})
    
    min_neg <- if(min(exp_pearson_neg) == 0) -1 else min(exp_pearson_neg)
    
    
    h1 <- Heatmap(exp_pearson_pos, col = colorRamp2(c(min_neg, 0, 1), c("red", "white", "blue")), column_title = "Positive",  column_title_gp = gpar(fontsize = 24, fontface = "bold"),
                  cluster_rows = FALSE, cluster_columns = FALSE, show_heatmap_legend = F, show_column_names = F, show_row_names = F)
    
    h2 <- Heatmap(exp_pearson_neg, col = colorRamp2(c(min_neg, 0, 1), c("red", "white", "blue")), column_title = paste("Negative - min val:", min(exp_pearson_neg)),  column_title_gp = gpar(fontsize = 24, fontface = "bold"),
                  cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = F, show_row_names = F)
    
    ht_list = h1 + h2
    png(paste( "chr-files/chr", chr, "exp_corrected.png", sep="_"), width = 2400, height = 1200)
    draw(ht_list, column_title = paste( "Chromosome", chr, "exp", "corrected", sep="_"),  column_title_gp = gpar(fontsize = 24, fontface = "bold"))
    dev.off()

  })


