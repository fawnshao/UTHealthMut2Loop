setwd("/Volumes/GoogleDrive/My Drive/hkg_tsg/v2.feature/test.ChIPs/")
library(data.table)
library(pheatmap)
loopcanfacs <- c("HDGF","CGGBP1","CTCF","YY1","ADNP","GABPA","MORC3",
                 "ZNF512","GATAD2A","ZBTB11","NRF1","ZNF143","MGA","MORC2",
                 "HMGB2","JARID2","PRDM10","ZNF281","ZNF296","POU5F1","ESRRB","NR0B1","DPPA2")
# range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}
range01 <- function(x, ...){(x - min(x, ...)) / (quantile(x, probs = 0.999, ...) - min(x, ...))}

#### if there is a peak
args <- c("highGM12878expr.1k.promoter.GM12878.ChIP.annos.tsv")
input <- fread(args[1], header = T, sep = "\t")
rowannos <- data.frame(Type = input$Type)
rownames(rowannos) <- input$ID
datamat <- data.matrix(input[,-c(1:2)])
rownames(datamat) <- input$ID
colors <- colorRampPalette(c("blue", "white", "red"))(100)

# png(filename = paste(args[1], "pheatmap.png", sep = "."), width = 2500, height = 1500)
datamat[datamat > 1] <- 1
pdf(file = paste(args[1], "pheatmap.pdf", sep = "."), width = 12, height = 10)
myplot <- pheatmap(datamat, scale = "none", annotation_row = rowannos,
                   show_rownames = F, show_colnames = T, color = colors, 
                   cluster_cols = T, cluster_rows = F, 
                   gaps_row = c(102, 827, 1558), fontsize_col = 5)
dev.off()

#### read depth
args <- c("hg19.GM12878.covereddepth.tsv")
input <- fread(args[1], header = T, sep = "\t")
rowannos <- data.frame(Type = input$Type)
rownames(rowannos) <- input$Gene
datamat <- data.matrix(input[,-c(1:2)])
rownames(datamat) <- input$Gene
colors <- colorRampPalette(c("blue", "white", "red"))(100)
rowannos$Type <- factor(rowannos$Type, levels = c("102_exprHKG.1", "729_exprHKG.2", 
                                                  "736_exprOther", "120_exprTSG", "762_exprEnh"))
mycols <- grep("HCFC1|TAF1|SP|ELK|ELF|ETS|KLF|SCL|MAX|E2F4|H3K|HDGF|CGGBP1|CTCF|YY1|ADNP|GABPA|MORC3|ZNF512|GATAD2A|ZBTB11|NRF1|ZNF143|MGA|MORC2|HMGB2|JARID2|PRDM10|ZNF281|ZNF296|POU5F1|ESRRB|NR0B1|DPPA2", colnames(datamat))
simdata <- datamat[,mycols]

# png(filename = paste(args[1], "pheatmap.png", sep = "."), width = 2500, height = 1500)
# datamat[datamat > 5] <- 5
datamat.scaled <- apply(datamat, 2, range01)
datamat.scaled[datamat.scaled > 1] <- 1
pdf(file = paste(args[1], "pheatmap.pdf", sep = "."), width = 16, height = 10)
myplot <- pheatmap(datamat.scaled, scale = "none", annotation_row = rowannos,
                   show_rownames = F, show_colnames = T, color = colors, 
                   cluster_cols = T, cluster_rows = F, 
                   gaps_row = c(102, 831, 1567, 1687), fontsize_col = 5)
dev.off()

# simdata[simdata > 30] <- 30
pdf(file = paste(args[1], "loopcanfacs.pheatmap.byrow.pdf", sep = "."), width = 6, height = 6)
myplot <- pheatmap(simdata, scale = "row", annotation_row = rowannos,
                   show_rownames = F, show_colnames = T, color = colors, 
                   cluster_cols = T, cluster_rows = F, 
                   gaps_row = c(102, 831, 1567, 1687), fontsize_col = 5)
dev.off()
pdf(file = paste(args[1], "loopcanfacs.pheatmap.bycol.pdf", sep = "."), width = 6, height = 6)
myplot <- pheatmap(simdata, scale = "column", annotation_row = rowannos,
                   show_rownames = F, show_colnames = T, color = colors, 
                   cluster_cols = T, cluster_rows = F, 
                   gaps_row = c(102, 831, 1567, 1687), fontsize_col = 5)
dev.off()

pdf(file = paste(args[1], "loopcanfacs.pheatmap.log2.pdf", sep = "."), width = 6, height = 6)
myplot <- pheatmap(log2(simdata + 1), scale = "column", annotation_row = rowannos,
                   show_rownames = F, show_colnames = T, color = colors, 
                   cluster_cols = T, cluster_rows = F, 
                   gaps_row = c(102, 831, 1567, 1687), fontsize_col = 5)
dev.off()

pdf(file = paste(args[1], "loopcanfacs.pheatmap.log2.none.pdf", sep = "."), width = 6, height = 6)
myplot <- pheatmap(log2(simdata + 1), scale = "none", annotation_row = rowannos,
                   show_rownames = F, show_colnames = T, color = colors, 
                   cluster_cols = T, cluster_rows = F, 
                   gaps_row = c(102, 831, 1567, 1687), fontsize_col = 5)
dev.off()

simdata.scaled <- apply(simdata, 2, range01)
simdata.scaled[simdata.scaled > 1] <- 1
pdf(file = paste(args[1], "loopcanfacs.pheatmap.range01.pdf", sep = "."), width = 6, height = 6)
myplot <- pheatmap(simdata.scaled, scale = "none", annotation_row = rowannos,
                   show_rownames = F, show_colnames = T, color = colors, 
                   cluster_cols = T, cluster_rows = F, 
                   clustering_distance_cols = "euclidean", clustering_method = "ward.D2", 
                   gaps_row = c(102, 831, 1567, 1687), fontsize_col = 5)
dev.off()

mytfs <- grep("YY1|E2F4|GABPA|NRF1|HCFC1|ELF1", colnames(simdata.scaled))
tfdata <- simdata.scaled[,mytfs]
pdf(file = paste(args[1], "tfs.pheatmap.range01.pdf", sep = "."), width = 6, height = 6)
myplot <- pheatmap(tfdata, scale = "none", annotation_row = rowannos,
                   show_rownames = F, show_colnames = T, color = colors, 
                   cluster_cols = T, cluster_rows = F, 
                   clustering_distance_cols = "euclidean", clustering_method = "ward.D2", 
                   gaps_row = c(102, 831, 1567, 1687), fontsize_col = 5)
dev.off()
pdf(file = paste(args[1], "tfs.pheatmap.range01.clst.pdf", sep = "."), width = 6, height = 6)
myplot <- pheatmap(tfdata, scale = "none", annotation_row = rowannos,
                   show_rownames = F, show_colnames = T, color = colors, 
                   cluster_cols = T, cluster_rows = T, 
                   clustering_distance_cols = "euclidean", 
                   clustering_distance_rows = "euclidean", clustering_method = "ward.D2", 
                   fontsize_col = 5)
dev.off()


#### K562
# E2F4 archived
#### read depth
args <- c("addE2F4.hg19.K562.covereddepth.tsv")
input <- fread(args[1], header = T, sep = "\t")
rowannos <- data.frame(Type = input$Type)
rownames(rowannos) <- input$Gene
datamat <- data.matrix(input[,-c(1:2)])
rownames(datamat) <- input$Gene
colors <- colorRampPalette(c("blue", "white", "red"))(100)
rowannos$Type <- factor(rowannos$Type, levels = c("84_exprHKG.1", "572_exprHKG.2", 
                                                  "439_exprOther", "57_exprTSG", "432_exprEnh"))
mycols <- grep("HCFC1|TAF1|SP|ELK|ELF|ETS|KLF|SCL|MAX|E2F4|H3K|HDGF|CGGBP1|CTCF|YY1|ADNP|GABPA|MORC3|ZNF512|GATAD2A|ZBTB11|NRF1|ZNF143|MGA|MORC2|HMGB2|JARID2|PRDM10|ZNF281|ZNF296|POU5F1|ESRRB|NR0B1|DPPA2", colnames(datamat))
simdata <- datamat[,mycols]

# png(filename = paste(args[1], "pheatmap.png", sep = "."), width = 2500, height = 1500)
# datamat[datamat > 5] <- 5
datamat.scaled <- apply(datamat, 2, range01)
datamat.scaled[datamat.scaled > 1] <- 1
pdf(file = paste(args[1], "pheatmap.pdf", sep = "."), width = 20, height = 10)
myplot <- pheatmap(datamat.scaled, scale = "none", annotation_row = rowannos,
                   show_rownames = F, show_colnames = T, color = colors, 
                   cluster_cols = T, cluster_rows = F, 
                   gaps_row = c(84, 656, 1095, 1152), fontsize_col = 4)
dev.off()

# simdata[simdata > 30] <- 30
pdf(file = paste(args[1], "loopcanfacs.pheatmap.byrow.pdf", sep = "."), width = 6, height = 6)
myplot <- pheatmap(simdata, scale = "row", annotation_row = rowannos,
                   show_rownames = F, show_colnames = T, color = colors, 
                   cluster_cols = T, cluster_rows = F, 
                   gaps_row = c(84, 656, 1095, 1152), fontsize_col = 5)
dev.off()
pdf(file = paste(args[1], "loopcanfacs.pheatmap.bycol.pdf", sep = "."), width = 6, height = 6)
myplot <- pheatmap(simdata, scale = "column", annotation_row = rowannos,
                   show_rownames = F, show_colnames = T, color = colors, 
                   cluster_cols = T, cluster_rows = F, 
                   gaps_row = c(84, 656, 1095, 1152), fontsize_col = 5)
dev.off()

pdf(file = paste(args[1], "loopcanfacs.pheatmap.log2.pdf", sep = "."), width = 6, height = 6)
myplot <- pheatmap(log2(simdata + 1), scale = "column", annotation_row = rowannos,
                   show_rownames = F, show_colnames = T, color = colors, 
                   cluster_cols = T, cluster_rows = F, 
                   gaps_row = c(84, 656, 1095, 1152), fontsize_col = 5)
dev.off()

pdf(file = paste(args[1], "loopcanfacs.pheatmap.log2.none.pdf", sep = "."), width = 6, height = 6)
myplot <- pheatmap(log2(simdata + 1), scale = "none", annotation_row = rowannos,
                   show_rownames = F, show_colnames = T, color = colors, 
                   cluster_cols = T, cluster_rows = F, 
                   gaps_row = c(84, 656, 1095, 1152), fontsize_col = 5)
dev.off()

simdata.scaled <- apply(simdata, 2, range01)
simdata.scaled[simdata.scaled > 1] <- 1
pdf(file = paste(args[1], "loopcanfacs.pheatmap.range01.pdf", sep = "."), width = 6, height = 6)
myplot <- pheatmap(simdata.scaled, scale = "none", annotation_row = rowannos,
                   show_rownames = F, show_colnames = T, color = colors, 
                   cluster_cols = T, cluster_rows = F, 
                   clustering_distance_cols = "euclidean", clustering_method = "ward.D2", 
                   gaps_row = c(84, 656, 1095, 1152), fontsize_col = 5)
dev.off()