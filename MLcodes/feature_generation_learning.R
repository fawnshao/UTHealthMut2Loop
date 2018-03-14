library(data.table)
library(ggplot2)
library(pheatmap)
library(cluster)
library(randomForest)
args <- commandArgs(TRUE)
# args <- c("GTEx.GTRD.mat","GTRD.maxcells.GTEx.hkg.tsg.txt")
# args <- c("housekeepinggene.GTRD.mat","GTRD.maxcells.GTEx.hkg.tsg.txt")
# scores <- as.matrix(read.table(args[1], sep = "\t", header = T, row.names = 1))
# class <- as.matrix(read.table(args[2], sep = "\t", header = T, row.names = 1))
# args <- c("pc.hkg.subsettsg.genes.HiC.HiChIP.GM.Liver.CAGE.Homer")
# args <- c("total.type.HiC.HiChIP.GTRD.homer.cage.sequenceFeatures")
# input <- read.table(args[1], sep = "\t", header = T, row.names = 1)
# args <- c("pc.hkg.subsettsg.genes.HiC.HiChIP.roadmap")
input <- fread(args[1], sep = "\t", header = T)
# write.csv(data.frame(input[,c(1,3:ncol(input)),with=FALSE],input[,2,with=FALSE]),"pc.hkg.subsettsg.genes.HiC.HiChIP.GM.Liver.CAGE.Homer.csv", row.names = F)
### notice: HiC.CH12-LX.is.mouse.B-lymphoblasts
### just remove it
scores <- data.matrix(input[,-c(1:3)])
genes <- as.matrix(input[,1])
class <- as.matrix(input[,2])
types <- factor(class)
rownames(scores) <- genes
#### PCA
# [class[,1]!="other",]

####
# colors <- colorRampPalette(c("blue", "white", "red"))(10)
# annotation_row <- data.frame(type = types)
# rownames(annotation_row) <- genes
# data <- data.matrix(input[,-c(1:2)])
# rownames(data) <- genes
# # ann_colors <- list(type = c("white", "firebrick"))
# png(filename = paste(args[1], "heatmap.png", sep = "."), width = 1500, height = 1200)
# myplot <- pheatmap(data, scale = "column", 
# 	show_rownames = F, show_colnames = T, color = colors, 
# 	annotation_row = annotation_row, 
# 	cluster_cols = T, cluster_rows = T)
# dev.off()
# png(filename = paste(args[1], "heatmap.hkgVStsg.png", sep = "."), width = 1500, height = 1200)
# myplot <- pheatmap(rbind.data.frame(data[types=="hkg",],data[types=="tsg",]), scale = "column", 
# 	show_rownames = F, show_colnames = T, color = colors, 
# 	cluster_cols = T, cluster_rows = F)
# dev.off()
# png(filename = paste(args[1], "heatmap.mansrt.png", sep = "."), width = 1500, height = 1200)
# myplot <- pheatmap(rbind.data.frame(data[types=="hkg",],data[types=="other",],data[types=="tsg",]), scale = "column", 
# 	show_rownames = F, show_colnames = T, color = colors, 
# 	cluster_cols = T, cluster_rows = F)
# dev.off()
####

myPCA <- prcomp(scores, scale. = T, center = T)
pca.scores <- as.data.frame(myPCA$x)
png(filename = paste(args[1], "pca.png", sep = "."), width = 1500, height = 1200)
ggplot(data = pca.scores, aes(x = PC1, y = PC2)) + #, label = rownames(pca.scores))) +
	geom_point(aes(colour = types, alpha = 1/10)) +
	ggtitle(args[1])
dev.off()
myPCA.t <- prcomp(t(scores), scale. = T, center = T)
pca.scores.t <- as.data.frame(myPCA.t$x)
png(filename = paste(args[1], "pca.t.png", sep = "."), width = 1500, height = 1500)
ggplot(data = pca.scores.t, aes(x = PC1, y = PC2, label = rownames(pca.scores.t))) +
	geom_point(aes(alpha = 1/10)) + geom_text(colour = "tomato", alpha = 0.8, size = 4) +
	ggtitle(args[1])
dev.off()
# feature.select <- order(abs(pca.scores.t[,1]),decreasing=T)[1:50]
# colnames(scores)[feature.select]

rawname <- colnames(scores)
colnames(scores) <- paste("var", 1:ncol(scores),sep = "")
scores <- scale(scores)
rf <- randomForest(types ~ ., data = scores)
a <- importance(rf)
rownames(a) <- rawname
b <- a[order(a, decreasing = T),]#[1:50]
mydata <- data.frame(names(b),b)
colnames(mydata) <- c("group", "importance")
mydata$group <- factor(mydata$group, levels = mydata$group)
png(filename = paste(args[1], "rfimportance.png", sep = "."), width = 1500, height = 600)
ggplot(data = mydata, aes(x = group, y = importance)) +
	geom_bar(stat="identity") + 
	ggtitle(args[1]) +
	theme(axis.text.x = element_text(angle = 60, hjust = 1))
dev.off()

c <- a[order(a, decreasing = T),]#[1:50]
colors <- colorRampPalette(c("blue", "white", "red"))(100)
# data <- scores[, order(a, decreasing = T)[1:50]]
colnames(scores) <- rawname
data <- scores[, order(a, decreasing = T)]
png(filename = paste(args[1], "topimportance.heatmap.png", sep = "."), width = 1500, height = 1200)
myplot <- pheatmap(data, scale = "none", 
	show_rownames = F, show_colnames = T, color = colors, 
	cluster_cols = F, cluster_rows = F)
dev.off()
# colors <- colorRampPalette(c("white", "red"))(2)
# data[data > 0] <- 1 
# png(filename = paste(args[1], "topimportance.heatmap.bin.png", sep = "."), width = 1500, height = 1200)
# myplot <- pheatmap(data, scale = "none", 
# 	show_rownames = F, show_colnames = T, color = colors, 
# 	cluster_cols = F, cluster_rows = F)
# dev.off()


#### boxplot for top 50 features
#### add 2 for the real column in raw input
# columns <- order(a, decreasing = T)[1:50] + 2
# columns <- order(a, decreasing = T)[1:50] + 3
# for(i in columns){
# 	datax <- input[,c(2,i),with=FALSE]
# 	colnames(datax) <- c("type", "value")
# 	png(filename = paste(args[1], i, "boxplot.png", sep = "."), width = 400, height = 400)
# 	ymax <- quantile(datax$value, probs = 0.95)
# 	myplot <- ggplot(data = datax, aes(x = type, y = value, fill = type)) + 
# 		geom_boxplot(position = position_dodge(1), outlier.alpha = 0.1, outlier.size = 0.1) + 
# 		scale_y_continuous(limits = c(0, ymax)) +
# 		ggtitle(colnames(input)[i]) + theme(legend.position = "none")
# 	print(myplot)
# 	dev.off()
# }
# save.image("a.RData")

# # seprate rf
# iris[1:4], iris$Species,



# rf.1 <- randomForest(types ~ ., data = scores)



# plot for hkg:tsg:other
###### for sequence feature
# for (i in 3:6){
# 	png(filename = paste(args[1], colnames(input)[i], "boxplot.png", sep = "."), width = 400, height = 400)
# 	datax <- input[,c(2,i),with = FALSE]
# 	colnames(datax) <- c("type", "value")
# 	ymax <- quantile(datax$value,probs = 0.95)
# 	myplot <- ggplot(data = datax, aes(x = type, y = value, fill = type)) + 
# 		geom_boxplot(position = position_dodge(1), outlier.alpha = 0.1, outlier.size = 0.1) + 
# 		scale_y_continuous(limits = c(0, ymax)) +
# 		ggtitle(colnames(input)[i]) + theme(legend.position = "none")
# 	print(myplot)
# 	dev.off()
# }











# i <- 206
# datax <- input[,c(2,i),with=FALSE]
# colnames(datax) <- c("type", "count")
# png(filename = paste(args[1], i, "boxplot.png", sep = "."), width = 400, height = 400)
# ggplot(data = datax, aes(x = type, y = count, fill = type)) + 
# 	geom_boxplot(position = position_dodge(1), outlier.alpha = 0.1, outlier.size = 0.1) + 
# 	ggtitle("count") + theme(legend.position = "none")
# dev.off()



# scores.bak <- scores
# scores[scores > 0] <- 1
# rf.b <- randomForest(types ~ ., data = scores)
# aa <- importance(rf)
# rownames(aa) <- rawname
# bb <- aa[order(aa, decreasing = T),][1:50]









# colors <- colorRampPalette(c("blue", "white", "red"))(100)
# data <- scores
# png(filename = "GTRD.maxexp.png", width = 1500, height = 1200)
# myplot <- pheatmap(data, scale = "column", 
# 	show_rownames = F, show_colnames = F, color = colors, 
# 	cluster_cols = F, cluster_rows = F)
# dev.off()

# colors <- colorRampPalette(c("white", "red"))(2)
# data[data > 0] <- 1
# torm <- sample(1:2506,2230)
# png(filename = "GTRD.maxexp.binary.png", width = 1500, height = 1200)
# myplot <- pheatmap(data[-torm,], scale = "none", 
# 	show_rownames = F, show_colnames = F, color = colors, 
# 	cluster_cols = F, cluster_rows = F)
# dev.off()

# subdata <- data[-torm,feature.select]
# annotation_row <- data.frame(type = types[-torm])
# rownames(annotation_row) <- genes[-torm]
# # ann_colors <- list(type = c("white", "firebrick"))
# png(filename = "subGTRD.maxexp.binary.png", width = 1500, height = 1200)
# myplot <- pheatmap(data[-torm,feature.select], scale = "none", 
# 	show_rownames = F, show_colnames = T, color = colors, 
# 	annotation_row = annotation_row, 
# 	cluster_cols = T, cluster_rows = T)
# dev.off()
# png(filename = "subGTRD.maxexp.binary.1.png", width = 1500, height = 1200)
# myplot <- pheatmap(data[-torm,feature.select], scale = "none", 
# 	show_rownames = F, show_colnames = T, color = colors, 
# 	annotation_row = annotation_row, 
# 	cluster_cols = T, cluster_rows = F)
# dev.off()


# myPCA.b <- prcomp(data, scale. = T, center = T)
# pca.scores.b <- as.data.frame(myPCA.b$x)
# png(filename = paste(args[1], "pca.binary.png", sep = "."), width = 1500, height = 1200)
# ggplot(data = pca.scores.b, aes(x = PC1, y = PC2)) + #, label = rownames(pca.scores))) +
# 	geom_point(aes(colour = types, alpha = 1/10)) +
# 	ggtitle("PCA plot of HiC.HiChIP.GM.Liver.CAGE.Homer")
# dev.off()
# # rownames(pca.scores.t[abs(pca.scores.t[,1])> 500 | abs(pca.scores.t[,2])> 500,])
# # rownames(pca.scores.t[order(abs(pca.scores.t[,1]),decreasing=T),1:2])[1:30]
# # order(abs(pca.scores.t[,1]),decreasing=T)[1:30]

# # kmeans(x,1)$withinss # trivial one-cluster, (its W.SS == ss(x))

# ## random starts do help here with too many clusters
# ## (and are often recommended anyway!):
# # subs <- sample(1:nrow(scores), 1000)
# cl <- kmeans(scores, 3, nstart = 25)
# png(filename = paste(args[1], "kmeans.3.png", sep = "."), width = 1500, height = 600)
# # plot(scores, col = cl$cluster)
# # points(cl$centers, col = 1:5, pch = 8)
# # labels = 2, 
# clusplot(scores, cl$cluster, color = TRUE, shade = TRUE, lines = 0)
# dev.off()
