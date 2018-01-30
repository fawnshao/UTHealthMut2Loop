args <- c("GTEx.gene.all.annotation.txt", "GTEx_sample.tissue.txt")
anno <- as.matrix(read.table(args[1], sep = "\t", header = T, row.names = 1))
info <- as.matrix(read.table(args[2], sep = "\t"))
tissues <- unique(info[,2])
data <- read.csv("GTEx.mean.tpm.tsv", header = T, row.names = 1)
output <- "GTEx.mean.10"
library(pheatmap)

breaklists <- c(seq(0, 10, by = 0.01), seq(10.1, 18, by = 0.1))
colorn <- length(breaklists)
colors <- colorRampPalette(c("blue", "yellow", "red"))(colorn)

anno[anno[,1] != "protein_coding",1] <- "ncRNA"
# row_annos <- data.frame(Class = factor(anno[,1]))
# rownames(row_annos) <- rownames(data)
# ann_colors = list(
#     Class = c(ncRNA = "#1B9E77", protein_coding = "#D95F02")
# )
dist_methods <- c("correlation", "euclidean", "maximum", "manhattan", "canberra", "minkowski")
clust_methods <- c("ward.D", "ward.D2", "complete", "average", "mcquitty", "median", "centroid")
i <- 2
j <- 2

coding <- data[anno[,1] == "protein_coding",]
png(filename = paste(output, dist_methods[i], clust_methods[j], "coding.png", sep = "."), width = 1800, height = 2000)
p1 <- pheatmap(log(coding + 1, 2), scale = "none", show_rownames = F, show_colnames = T, 
         color = colors, cluster_cols = F, clustering_distance_rows = dist_methods[i], 
         clustering_method = clust_methods[j], breaks = breaklists,
         cutree_row = 10)
# p1 <- pheatmap(log(data + 1, 2), scale = "none", show_rownames = F, show_colnames = T, 
#          color = colors, cluster_cols = T, clustering_distance_rows = dist_methods[i], 
#          clustering_distance_cols = dist_methods[i], 
#          clustering_method = clust_methods[j], breaks = breaklists,
#          cutree_row = 8)
dev.off()
codingcluster <- cutree(p1$tree_row, k = 10)
write.table(data.frame(codingcluster, rownames(coding), coding), file = paste(output, "coding.pheatmap.tsv", sep = "."), sep = "\t", row.names = FALSE, quote = FALSE)

noncoding <- data[anno[,1] != "protein_coding",]
png(filename = paste(output, dist_methods[i], clust_methods[j], "noncoding.png", sep = "."), width = 1800, height = 2000)
p2 <- pheatmap(log(noncoding + 1, 2), scale = "none", show_rownames = F, show_colnames = T, 
         color = colors, cluster_cols = F, clustering_distance_rows = dist_methods[i], 
         clustering_method = clust_methods[j], breaks = breaklists,
         cutree_row = 10)
dev.off()
noncodingcluster <- cutree(p2$tree_row, k = 10)
write.table(data.frame(noncodingcluster, rownames(noncoding), noncoding), file = paste(output, "noncoding.pheatmap.tsv", sep = "."), sep = "\t", row.names = FALSE, quote = FALSE)

save.image("GTEx.2.RData")

# randrows <- sample(1:nrow(data),1000)
# png(filename = paste(output, dist_methods[i], clust_methods[j], "r1000.png", sep = "."), width = 1800, height = 2000)
# pheatmap(log(data[randrows,] + 1, 2), scale = "none", show_rownames = F, show_colnames = T, 
#          color = colors, cluster_cols = F, clustering_distance_rows = dist_methods[i], 
#          clustering_method = clust_methods[j], breaks = breaklists,
#          annotation_row = row_annos[randrows,], annotation_colors = ann_colors)
# dev.off()
