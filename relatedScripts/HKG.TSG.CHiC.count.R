library(pheatmap)
args <- commandArgs(TRUE)
data <- as.matrix(read.table(args[1], sep = "\t", row.names = 1, header = T, na.strings = "/"))
data[is.na(data)] <- 0
# breaklists <- c(seq(0, 2, by = 0.01),seq(2.1, 5.7, by = 0.1))
# colorn <- length(breaklists)
# colors <- colorRampPalette(c("blue", "yellow", "red"))(colorn)
colors <- colorRampPalette(c("blue", "yellow", "red"))(100)
x <- log2(data + 1)
if(args[2] == "y"){
	n <- ncol(data)
	x <- data.frame(log2(data[,-n] + 1), data[,n])
}
png(filename = paste(args[1], "nocluster.pheatmap.png", sep = "."), width = 800, height = 1000)
pheatmap(x, scale = "none", show_rownames = F, show_colnames = T, 
         color = colors, cluster_rows = F, cluster_cols = F
         )
dev.off()

png(filename = paste(args[1], "pheatmap.png", sep = "."), width = 800, height = 1000)
p1 <- pheatmap(x, scale = "none", show_rownames = F, show_colnames = T, 
         color = colors, 
         clustering_distance_cols = "euclidean", clustering_distance_rows = "euclidean", 
         clustering_method = "ward.D2"
         )
dev.off()
cluster <- cutree(p1$tree_row, k = 10)
write.table(data.frame(cluster[p1$tree_row$order], 
	rownames(data)[p1$tree_row$order],
	data[p1$tree_row$order, p1$tree_col$order]), 
	file = paste(args[1], "pheatmap.tsv", sep = "."), 
	sep = "\t", row.names = FALSE, quote = FALSE)
