library(pheatmap)
args <- commandArgs(TRUE)
data <- as.matrix(read.table(args[1], sep = "\t", row.names = 1, header = T, na.strings = "/"))
x <- data
x[is.na(x)] <- 0
x[x>0] <- 1
# breaklists <- c(seq(0, 2, by = 0.01),seq(2.1, 5.7, by = 0.1))
# colorn <- length(breaklists)
# colors <- colorRampPalette(c("blue", "yellow", "red"))(colorn)
colors <- colorRampPalette(c("blue", "yellow", "red"))(100)

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

### the following is same
# png(filename = paste(args[1], "a.png", sep = "."), width = 800, height = 1000)
# p1 <- pheatmap(x, scale = "none", show_rownames = F, show_colnames = T, 
#          color = colors, 
#          clustering_distance_cols = "correlation", clustering_distance_rows = "correlation", 
#          clustering_method = "average"
#          )
# dev.off()

# library(gplots)
# hc <- hclust(as.dist(1 - cor(t(x))), method="average")
# hc2 <- hclust(as.dist(1 - cor(x)), method="average")
# png(filename = paste(args[1], "b.png", sep = "."), width = 800, height = 1000)
# heatmap(as.matrix(x), Rowv = as.dendrogram(hc), Colv = as.dendrogram(hc2), 
# 	col=colors, labRow="", scale = "none")
# dev.off()


### use loop to get the best pattern
# dist_methods <- c("correlation", "euclidean", "maximum", "manhattan", "canberra", "minkowski")
# clust_methods <- c("ward.D", "ward.D2", "complete", "average", "mcquitty", "median", "centroid")
# for(i in 1:length(dist_methods)){
# 	for(j in 1:length(clust_methods)){
# 		png(filename = paste(args[1], dist_methods[i], clust_methods[j], "png", sep = "."), width = 1200, height = 1000)
# 		p1 <- pheatmap(x, scale = "none", 
# 			show_rownames = F, show_colnames = T, color = colors, 
# 		    clustering_distance_cols = dist_methods[i], 
# 		    clustering_distance_rows = dist_methods[i], 
# 		    clustering_method = clust_methods[j]
# 		    )
# 		dev.off()
# 	}
# }

### cases
# library(pheatmap)
# args <- c("housekeepinggene.Cell_Javierre_17cells.bait.oe.mean.HPA.tsv","34")
# args <- c("tissuespecificgene.Cell_Javierre_17cells.bait.oe.mean.HPA.tsv", "34")
# data <- as.matrix(read.table(args[1], sep = "\t", row.names = 1, header = T, na.strings = "/"))
# data[is.na(data)] <- 0
# colors <- colorRampPalette(c("blue", "yellow", "red"))(100)
# count <- as.numeric(args[2])
# n <- ncol(data)
# x <- data.frame(log2(data[,1:count] + 1), data[,c((count + 1):n)])

# dist_methods <- c("correlation", "euclidean", "maximum", "manhattan", "canberra", "minkowski")
# clust_methods <- c("ward.D", "ward.D2", "complete", "average", "mcquitty", "median", "centroid")
# for(i in 1:length(dist_methods)){
# 	for(j in 1:length(clust_methods)){
# 		png(filename = paste(args[1], dist_methods[i], clust_methods[j], "png", sep = "."), width = 1200, height = 1000)
# 		p1 <- pheatmap(x, scale = "none", 
# 			show_rownames = F, show_colnames = T, color = colors, 
# 		    clustering_distance_cols = dist_methods[i], 
# 		    clustering_distance_rows = dist_methods[i], 
# 		    clustering_method = clust_methods[j]
# 		    )
# 		dev.off()
# 	}
# }

# a <- data[,1:34]
# data[rowSums(a)==0,35]
# sort(data[rowSums(a)==0,35])