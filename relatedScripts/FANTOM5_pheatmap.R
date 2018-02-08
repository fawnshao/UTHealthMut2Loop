args <- c("fantom5.v19.promoter.within1k.div.pairs.txt", "fantom5.v19.promoter.within1k.div.pairs.tsv")
anno <- as.matrix(read.table(args[1], sep = "\t"))
data <- as.matrix(read.table(args[2], sep = "\t", row.names = 1))
output <- "fantom5.div.tissues"
library(pheatmap)

breaklists <- c(seq(0, 10, by = 0.01), seq(10.1, 15, by = 0.1))
colorn <- length(breaklists)
colors <- colorRampPalette(c("blue", "yellow", "red"))(colorn)

anno[anno[,4] != "protein-coding",4] <- "ncRNA"
anno[anno[,5] != "protein-coding",5] <- "ncRNA"
anno[anno[,4] == "protein-coding",4] <- "proteincoding"
anno[anno[,5] == "protein-coding",5] <- "proteincoding"
row_annos <- data.frame(fwd = factor(anno[,5]), rev = factor(anno[,4]))
rownames(row_annos) <- rownames(data)
ann_colors = list(fwd = c(ncRNA = "#1B9E77", proteincoding = "#D95F02"),
	rev = c(ncRNA = "#1B9E77", proteincoding = "#D95F02"))
dist_methods <- c("correlation", "euclidean", "maximum", "manhattan", "canberra", "minkowski")
clust_methods <- c("ward.D", "ward.D2", "complete", "average", "mcquitty", "median", "centroid")
i <- 2
j <- 2
png(filename = paste(output, dist_methods[i], clust_methods[j], "png", sep = "."), width = 1800, height = 2000)
p1 <- pheatmap(log(data + 1, 2), scale = "none", show_rownames = F, show_colnames = F, 
         color = colors, cluster_cols = F, clustering_distance_rows = dist_methods[i], 
         clustering_method = clust_methods[j], breaks = breaklists,
         annotation_row = row_annos, annotation_colors = ann_colors,
         cutree_row = 10)
dev.off()
cluster <- cutree(p1$tree_row, k = 10)
rev.mean <- apply(data[,1:188], 1, mean)
fwd.mean <- apply(data[,189:376], 1, mean)
##
anno <- as.matrix(read.table(args[1], sep = "\t"))
##
write.table(data.frame(cluster[p1$tree_row$order], 
	rownames(data)[p1$tree_row$order],anno[p1$tree_row$order,],
	rev.mean[p1$tree_row$order], fwd.mean[p1$tree_row$order],
	data[p1$tree_row$order,]), 
	file = paste(output, "pheatmap.tsv", sep = "."), 
	sep = "\t", row.names = FALSE, quote = FALSE)
save.image("fantom5.RData")

res <- matrix(ncol = 4, nrow = nrow(data))
for(i in 1:nrow(data)){
	x1 <- data[i,1:188]
	x2 <- data[i,189:376]
	cor.x <- cor.test(x1, x2, method = "pearson")
	cor.y <- cor.test(x1, x2, method = "spearman")
	res[i,] <- c(cor.x$estimate[[1]], cor.x$p.value, cor.y$estimate[[1]], cor.y$p.value)
}
write.table(data.frame(cluster[p1$tree_row$order], 
	rownames(data)[p1$tree_row$order],anno[p1$tree_row$order,],
	rev.mean[p1$tree_row$order], fwd.mean[p1$tree_row$order],
	res[p1$tree_row$order,]), 
	file = paste(output, "pheatmap.cor.tsv", sep = "."), 
	sep = "\t", row.names = FALSE, quote = FALSE)

# coding <- data[anno[,1] == "protein_coding",]
# png(filename = paste(output, dist_methods[i], clust_methods[j], "coding.png", sep = "."), width = 1800, height = 2000)
# p1 <- pheatmap(log(coding + 1, 2), scale = "none", show_rownames = F, show_colnames = T, 
#          color = colors, cluster_cols = F, clustering_distance_rows = dist_methods[i], 
#          clustering_method = clust_methods[j], breaks = breaklists,
#          cutree_row = 10)
# # p1 <- pheatmap(log(data + 1, 2), scale = "none", show_rownames = F, show_colnames = T, 
# #          color = colors, cluster_cols = T, clustering_distance_rows = dist_methods[i], 
# #          clustering_distance_cols = dist_methods[i], 
# #          clustering_method = clust_methods[j], breaks = breaklists,
# #          cutree_row = 8)
# dev.off()
# codingcluster <- cutree(p1$tree_row, k = 10)
# write.table(data.frame(codingcluster, rownames(coding), coding), file = paste(output, "coding.pheatmap.tsv", sep = "."), sep = "\t", row.names = FALSE, quote = FALSE)

# noncoding <- data[anno[,1] != "protein_coding",]
# png(filename = paste(output, dist_methods[i], clust_methods[j], "noncoding.png", sep = "."), width = 1800, height = 2000)
# p2 <- pheatmap(log(noncoding + 1, 2), scale = "none", show_rownames = F, show_colnames = T, 
#          color = colors, cluster_cols = F, clustering_distance_rows = dist_methods[i], 
#          clustering_method = clust_methods[j], breaks = breaklists,
#          cutree_row = 10)
# dev.off()
# noncodingcluster <- cutree(p2$tree_row, k = 10)
# write.table(data.frame(noncodingcluster, rownames(noncoding), noncoding), file = paste(output, "noncoding.pheatmap.tsv", sep = "."), sep = "\t", row.names = FALSE, quote = FALSE)

