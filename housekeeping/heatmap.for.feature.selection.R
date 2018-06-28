args <- commandArgs(TRUE)
library(pheatmap)
library(data.table)
input <- fread(args[1], sep = "\t", header = T, na.strings = "/")
scores <- data.matrix(input[,-c(1:3)])
colors <- colorRampPalette(c("white", "red"))(10)
rownames(scores) <- as.matrix(input[,1])[,1]
annos_row <- as.data.frame(input[,2:3])
rownames(annos_row) <- rownames(scores)

png(filename = paste(args[1], "pheatmap.png", sep = "."), width = 1800, height = 1000)
myplot <- pheatmap(scores, scale = "column", annotation_row = annos_row, 
	show_rownames = F, show_colnames = T, color = colors,
	cluster_cols = F, cluster_rows = F)
dev.off()
