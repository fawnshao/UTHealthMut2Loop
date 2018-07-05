args <- commandArgs(TRUE)
# args <- c("v2.hkg.tsg.vert.known.motifs.mat")
library(pheatmap)
library(data.table)
input <- fread(args[1], sep = "\t", header = T, na.strings = "/")
scores <- data.matrix(input[,-c(1:3)])
colors <- colorRampPalette(c("white", "red"))(10)
rownames(scores) <- as.matrix(input[,1])[,1]
annos_row <- as.data.frame(input[,2:3])
rownames(annos_row) <- rownames(scores)
scores[scores > 10] <- 10

png(filename = paste(args[1], "pheatmap.png", sep = "."), width = 2500, height = 2000)
myplot <- pheatmap(scores, scale = "none", annotation_row = annos_row, 
	show_rownames = F, show_colnames = T, color = colors,
	cluster_cols = F, cluster_rows = F)
dev.off()
