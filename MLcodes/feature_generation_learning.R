library(data.table)
library(ggplot2)
library(pheatmap)
library(cluster)
args <- commandArgs(TRUE)
# args <- c("GTEx.GTRD.mat","GTRD.maxcells.GTEx.hkg.tsg.txt")
# args <- c("housekeepinggene.GTRD.mat","GTRD.maxcells.GTEx.hkg.tsg.txt")
scores <- as.matrix(read.table(args[1], sep = "\t", header = T, row.names = 1))
class <- as.matrix(read.table(args[2], sep = "\t", header = T, row.names = 1))
types <- factor(class)
#### PCA
myPCA <- prcomp(scores, scale. = T, center = T)
pca.scores <- as.data.frame(myPCA$x)
png(filename = paste(args[1], "pca.png", sep = "."), width = 1500, height = 600)
ggplot(data = pca.scores[class[,1]!="other",], aes(x = PC1, y = PC2)) + #, label = rownames(pca.scores))) +
	geom_point(aes(colour = types[class[,1]!="other"], alpha = 1/10)) +
	ggtitle("PCA plot of GTRD maxcells")
dev.off()





# kmeans(x,1)$withinss # trivial one-cluster, (its W.SS == ss(x))

## random starts do help here with too many clusters
## (and are often recommended anyway!):
# subs <- sample(1:nrow(scores), 1000)
cl <- kmeans(scores, 3, nstart = 25)
png(filename = paste(args[1], "kmeans.3.png", sep = "."), width = 1500, height = 600)
# plot(scores, col = cl$cluster)
# points(cl$centers, col = 1:5, pch = 8)
# labels = 2, 
clusplot(scores, cl$cluster, color = TRUE, shade = TRUE, lines = 0)
dev.off()
