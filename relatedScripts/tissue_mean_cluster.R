args <- commandArgs(TRUE)
library(pheatmap)
# args <- c("GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct", "GTEx_sample.tissue.txt")
exprs <- as.matrix(read.table(args[1], sep = "\t", row.names = 1, header = T))
info <- as.matrix(read.table(args[2], sep = "\t"))
tissues <- unique(info[,2])
ave.tissues <- matrix(ncol = length(tissues), nrow = nrow(exprs))
output <- "GTEx.mean"

for(i in 1:nrow(exprs)){
	for(j in 1:length(tissues)){
		ave.tissues[i,j] <- mean(exprs[i,info[,2]==tissues[j]])
	}
}
colnames(ave.tissues) <- tissues
rownames(ave.tissues) <- rownames(exprs)
write.csv(ave.tissues, file = paste(output, "tpm.tsv", sep = "."))
# save.image("GTEx.RData")

args <- c("GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct", "GTEx_sample.tissue.txt")
info <- as.matrix(read.table(args[2], sep = "\t"))
tissues <- unique(info[,2])
ave.tissues <- read.csv("GTEx.mean.tpm.tsv", header = T, row.names = 1)
output <- "GTEx.mean"
library(pheatmap)

nacount <- apply(ave.tissues, 1, function(x){sum(is.na(x))})
data <- ave.tissues[nacount < length(tissues) / 2,]
breaklists <- c(seq(0, 8, by = 0.1), seq(9, 18, by = 1))
colorn <- length(breaklists)
colors <- colorRampPalette(c("blue", "white", "red"))(colorn)

dist_methods <- c("correlation", "euclidean", "maximum", "manhattan", "canberra", "minkowski")
clust_methods <- c("ward.D", "ward.D2", "complete", "average", "mcquitty", "median", "centroid")
i <- 2
j <- 2
png(filename = paste(output, dist_methods[i], clust_methods[j], "png", sep = "."), width = 1500, height = 2000)
p1 <- pheatmap(log(data + 1, 2), scale = "none", show_rownames = F, show_colnames = T, 
         color = colors, cluster_cols = F, clustering_distance_rows = dist_methods[i], 
         clustering_method = clust_methods[j], breaks = breaklists,
         cutree_row = 8)
# p1 <- pheatmap(log(data + 1, 2), scale = "none", show_rownames = F, show_colnames = T, 
#          color = colors, cluster_cols = T, clustering_distance_rows = dist_methods[i], 
#          clustering_distance_cols = dist_methods[i], 
#          clustering_method = clust_methods[j], breaks = breaklists,
#          cutree_row = 8)
dev.off()

pairscluster <- cutree(p1$tree_row, k = 8)
write.table(data.frame(pairscluster, rownames(data), data), file = paste(output, "mean.pheatmap.tsv", sep = "."), sep = "\t", row.names = FALSE, quote = FALSE)

save.image("GTEx.1.RData")

