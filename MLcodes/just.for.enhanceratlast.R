# .libPaths("/work/04935/shaojf/stampede2/myTools/anaconda2/lib/R/library")
library(data.table)
library(ggplot2)
library(pheatmap)
# args <- commandArgs(TRUE)
args <- c("hkg.tsg.srtbyPCA.class.enhanceratlas")
input <- fread(args[1], sep = "\t", header = T)
### for new data.table
scores <- data.matrix(input[,-c(1:2)])
genes <- as.matrix(input[,1])
class <- as.matrix(input[,2])
types <- factor(class)
rownames(scores) <- genes
types.sim <- as.vector(types)
types.sim[1:2551] <- "HKG"
types.sim[2552:5286] <- "singleTSG"
types.sim <- as.factor(types.sim)

colnames(scores) <- gsub("-","_",colnames(scores))
enhanceratlas <- scores
enhanceratlas01 <- scores
enhanceratlas01[enhanceratlas01 > 0] <- 1
enhanceratlas <- apply(enhanceratlas, 2, as.numeric)
rownames(enhanceratlas) <- genes

############## FeatureSelection
RunFS <- function(x, y = as.numeric(types), pre = args[1]){
	require(FeatureSelection)
	trainx <- x
	trainy <- y
	pre <- pre
	params_glmnet <- list(alpha = 1, family = 'gaussian', nfolds = 5, parallel = TRUE)
	params_xgboost <- list(params = list("objective" = "reg:linear", "bst:eta" = 0.001, "subsample" = 0.75, "max_depth" = 5, "colsample_bytree" = 0.75, "nthread" = 60), nrounds = 1000, print.every.n = 250, maximize = FALSE)
	params_ranger <- list(dependent.variable.name = 'y', probability = FALSE, num.trees = 1000, verbose = TRUE, mtry = 5, min.node.size = 10, num.threads = 60, classification = FALSE, importance = 'permutation')
	params_features <- list(keep_number_feat = NULL, union = TRUE)
	feat <- wrapper_feat_select(X = trainx, y = trainy, params_glmnet = params_glmnet, params_xgboost = params_xgboost, params_ranger = params_ranger, xgb_sort = 'Gain', CV_folds = 5, stratified_regr = FALSE, scale_coefs_glmnet = FALSE, cores_glmnet = 5, params_features = params_features, verbose = TRUE)
	write.table(feat$all_feat$"glmnet-lasso", paste(pre, "FeatureSelection.glmnetlasso.tsv", sep = "."), row.names = F, sep = "\t")
	write.table(feat$all_feat$xgboost, paste(pre, "FeatureSelection.xgboost.tsv", sep = "."), row.names = F, sep = "\t")
	write.table(feat$all_feat$ranger, paste(pre, "FeatureSelection.ranger.tsv", sep = "."), row.names = F, sep = "\t")
	write.table(feat$union_feat, paste(pre, "FeatureSelection.union_feat.tsv", sep = "."), row.names = F, sep = "\t")
	return(feat)
}
enhanceratlas.feat <- RunFS(x = enhanceratlas, y = as.numeric(types), pre = paste(args[1],"enhanceratlas", sep = "."))
enhanceratlas01.feat <- RunFS(x = enhanceratlas01, y = as.numeric(types), pre = paste(args[1],"enhanceratlas01", sep = "."))
enhanceratlas.sim.feat <- RunFS(x = enhanceratlas, y = as.numeric(types.sim), pre = paste(args[1],"enhanceratlas.sim", sep = "."))
enhanceratlas01.sim.feat <- RunFS(x = enhanceratlas01, y = as.numeric(types.sim), pre = paste(args[1],"enhanceratlas01.sim", sep = "."))
save.image(paste(args[1], "FS.RData", sep = "."))
############## PCA
RunPCA <- function(x, y = types, z = types.sim, pre = args[1]){
	require(ggplot2)
	data <- x
	mycolors <- y
	mycolors.sim <- z
	pre <- pre
	myPCA.a <- prcomp(data, scale. = F, center = F)
	myPCA.a.t <- prcomp(t(data), scale. = F, center = F)
	pca.scores <- as.data.frame(myPCA.a$x)
	pdf(file = paste(pre, "pca.pdf", sep = "."), width = 15, height = 12)
	p1 <- ggplot(data = pca.scores, aes(x = PC1, y = PC2)) +
		geom_point(aes(colour = mycolors, alpha = 1/10)) +
		ggtitle(pre)
	print(p1)
	dev.off()
	pdf(file = paste(pre, "pca.sim.pdf", sep = "."), width = 15, height = 12)
	p1 <- ggplot(data = pca.scores, aes(x = PC1, y = PC2)) +
		geom_point(aes(colour = mycolors.sim, alpha = 1/10)) +
		ggtitle(pre)
	print(p1)
	dev.off()

	pca.scores <- as.data.frame(myPCA.a.t$x)
	pdf(file = paste(pre, "t.pca.pdf", sep = "."), width = 15, height = 12)
	p1 <- ggplot(data = pca.scores, aes(x = PC1, y = PC2, label = rownames(pca.scores))) +
		geom_point(aes(alpha = 1/10)) +
		geom_text(colour = "tomato", alpha = 0.8, size = 2) +
		ggtitle(args[1])
	print(p1)
	dev.off()
}
enhanceratlas.pca <- RunPCA(x = enhanceratlas, pre = paste(args[1],"enhanceratlas", sep = "."))
enhanceratlas01.pca <- RunPCA(x = enhanceratlas01, pre = paste(args[1],"enhanceratlas01", sep = "."))
save.image(paste(args[1], "FS.PCA.RData", sep = "."))

############## t-SNE
# conda.R
args <- c("hkg.tsg.srtbyPCA.all")
load(paste(args[1], "detail.data.RData", sep = "."))
.libPaths("/work/04935/shaojf/stampede2/myTools/anaconda2/lib/R/library")
library(ggplot2)
library(Rtsne)
set.seed(123)
RuntSNE <- function(x, pre = args[1]){
	data <- x
	pre <- pre
	tsne <- Rtsne(data, check_duplicates = FALSE, dims = 2, perplexity = 30, verbose = TRUE, max_iter = 500)
	v.tsne <- as.data.frame(tsne$Y)
	colnames(v.tsne) <- c("tSNE.1", "tSNE.2")
	write.table(v.tsne, paste(pre, "tSNE.tsv", sep = "."), row.names = T, sep = "\t")
	return(tsne)
}
enhanceratlas.tsne <- RuntSNE(x = enhanceratlas, pre = paste(args[1],"enhanceratlas", sep = "."))
enhanceratlas01.tsne <- RuntSNE(x = enhanceratlas01, pre = paste(args[1],"enhanceratlas01", sep = "."))

PlottSNE <- function(x, y = types, z = types.sim, pre = args[1]){
	require(ggplot2)
	data <- x
	mycolors <- y
	mycolors.sim <- z
	pre <- pre
	v.tsne <- as.data.frame(data$Y)
	colnames(v.tsne) <- c("tSNE.1", "tSNE.2")
	pdf(file = paste(pre, "tsne.pdf", sep = "."), width = 15, height = 12)
	p1 <- ggplot(data = v.tsne, aes(x = tSNE.1, y = tSNE.2)) + 
		geom_point(aes(colour = mycolors, alpha = 1/10)) +
		ggtitle(pre)
	print(p1)
	dev.off()
	pdf(file = paste(pre, "tsne.sim.pdf", sep = "."), width = 15, height = 12)
	p1 <- ggplot(data = v.tsne, aes(x = tSNE.1, y = tSNE.2)) + 
		geom_point(aes(colour = mycolors.sim, alpha = 1/10)) +
		ggtitle(pre)
	print(p1)
	dev.off()
}
enhanceratlas.tsne.plot <- PlottSNE(x = enhanceratlas.tsne, pre = paste(args[1],"enhanceratlas", sep = "."))
enhanceratlas01.tsne.plot <- PlottSNE(x = enhanceratlas01.tsne, pre = paste(args[1],"enhanceratlas01", sep = "."))
save.image(paste(args[1], "tSNE.RData", sep = "."))

# please also ref to: 
# http://r-statistics.co/Variable-Selection-and-Importance-With-R.html
###################################################
#####           look at the above             #####
###################################################

############## visualization of the important features
# cat onlyhkg.tsg.further.srt.all.dnase.FeatureSelection.union_feat.tsv onlyhkg.tsg.further.srt.all.pchic.FeatureSelection.union_feat.tsv onlyhkg.tsg.further.srt.all.hic.FeatureSelection.union_feat.tsv onlyhkg.tsg.further.srt.all.seqfeature.FeatureSelection.union_feat.tsv onlyhkg.tsg.further.srt.all.homermotif.FeatureSelection.union_feat.tsv onlyhkg.tsg.further.srt.all.meth.FeatureSelection.union_feat.tsv onlyhkg.tsg.further.srt.all.gtrdpeaks.FeatureSelection.union_feat.tsv onlyhkg.tsg.further.srt.all.histone.FeatureSelection.union_feat.tsv | awk '$2>0.8' > important.features.tsv
# cat onlyhkg.tsg.further.srt.all.dnase.sim.FeatureSelection.union_feat.tsv onlyhkg.tsg.further.srt.all.pchic.sim.FeatureSelection.union_feat.tsv onlyhkg.tsg.further.srt.all.hic.sim.FeatureSelection.union_feat.tsv onlyhkg.tsg.further.srt.all.seqfeature01.sim.FeatureSelection.union_feat.tsv onlyhkg.tsg.further.srt.all.homermotif.sim.FeatureSelection.union_feat.tsv onlyhkg.tsg.further.srt.all.meth.sim.FeatureSelection.union_feat.tsv onlyhkg.tsg.further.srt.all.gtrdpeaks.sim.FeatureSelection.union_feat.tsv onlyhkg.tsg.further.srt.all.histone.sim.FeatureSelection.union_feat.tsv | awk '$2>0.8' > sim.important.features.tsv
args <- c("hkg.tsg.srtbyPCA.all")
load(paste(args[1], "detail.data.RData", sep = "."))
feats <- read.table("important.features.tsv", sep = "\t", row.names = 1)
library(pheatmap)
data <- scores[,colnames(scores) %in% rownames(feats)]
# data.1 <- scores[,match(rownames(feats),colnames(scores))]
range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}
data.x <- apply(data, 2, range01)
annos <- data.frame(class = types)
rownames(annos) <- rownames(data.x)
colors <- colorRampPalette(c("blue", "white", "red"))(100)
png(filename = paste(args[1], "important.pheatmap.png", sep = "."), width = 2500, height = 1200)
myplot <- pheatmap(data.x, scale = "none", annotation_row = annos,
	show_rownames = F, show_colnames = T, color = colors, 
	cluster_cols = F, cluster_rows = F)
dev.off()
pdf(file = paste(args[1], "important.pheatmap.pdf", sep = "."), width = 25, height = 12)
myplot <- pheatmap(data.x, scale = "none", annotation_row = annos,
	show_rownames = F, show_colnames = T, color = colors, 
	cluster_cols = F, cluster_rows = F)
dev.off()
#########
feats <- read.table("sim.important.features.tsv", sep = "\t", row.names = 1)
library(pheatmap)
data <- scores[,colnames(scores) %in% rownames(feats)]
# data.1 <- scores[,match(rownames(feats),colnames(scores))]
range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}
data.x <- apply(data, 2, range01)
annos <- data.frame(class = types)
rownames(annos) <- rownames(data.x)
colors <- colorRampPalette(c("blue", "white", "red"))(100)
png(filename = paste(args[1], "sim.important.pheatmap.png", sep = "."), width = 2500, height = 1200)
myplot <- pheatmap(data.x, scale = "none", annotation_row = annos,
	show_rownames = F, show_colnames = T, color = colors, 
	cluster_cols = F, cluster_rows = F)
dev.off()
pdf(file = paste(args[1], "sim.important.pheatmap.pdf", sep = "."), width = 25, height = 12)
myplot <- pheatmap(data.x, scale = "none", annotation_row = annos,
	show_rownames = F, show_colnames = T, color = colors, 
	cluster_cols = F, cluster_rows = F)
dev.off()
#########
# all
data <- scores
range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}
data.x <- apply(data, 2, range01)
annos <- data.frame(class = types)
rownames(annos) <- rownames(data.x)
colors <- colorRampPalette(c("blue", "white", "red"))(100)
png(filename = paste(args[1], "pheatmap.png", sep = "."), width = 2500, height = 1200)
myplot <- pheatmap(data.x, scale = "none", annotation_row = annos,
	show_rownames = F, show_colnames = T, color = colors, 
	cluster_cols = F, cluster_rows = F)
dev.off()
pdf(file = paste(args[1], "pheatmap.pdf", sep = "."), width = 25, height = 12)
myplot <- pheatmap(data.x, scale = "none", annotation_row = annos,
	show_rownames = F, show_colnames = T, color = colors, 
	cluster_cols = F, cluster_rows = F)
dev.off()
#########




library(data.table)
library(pheatmap)
args <- c("hkg.tsg.srtbyPCA.class.part2")
input <- fread(args[1], sep = "\t", header = T)
### for new data.table
# 1 Gene
# 2 Type
# 3 Length
# 4 ExonCount
# 5 ExonLength
# 6 IntronLenth
# 7 FirstExonLength
# 8 FirstIntronLength
# 9 FirstExonCons
# 10 FirstIntronCons
# 11 PromoterCons
# 12 cpgIslandLength
# 13 simpleRepeatLength
# 14 CAGE.count

# scores <- data.matrix(input[,-c(1:2)])
scores <- data.matrix(input[,c(4,5,12,14,15:29,81:153)])
genes <- as.matrix(input[,1])
class <- as.matrix(input[,2])
types <- factor(class)
rownames(scores) <- genes
subs <- rbind.data.frame(scores[scores[,91]==0 & scores[,92]==0 & class=="hkg1",],
	scores[scores[,91]==0 & scores[,92]==0 & class=="hkg2",],
	scores[scores[,91]==0 & scores[,92]==0 & class=="hkg3",],
	scores[scores[,91]==0 & scores[,92]==0 & class=="hkg4",]
	)



scores[scores[,1]>100,1] <- 100
scores[scores[,2]>5000,2] <- 5000
scores[scores[,3]>1000,3] <- 1000
scores[scores[,4]>10,4] <- 10
a <- scores[,5:19]
a[a > 5] <- 5
scores[,5:19] <- a
a <- scores[,20:90]
a[a > 10] <- 10
scores[,20:90] <- a
a <- scores[,91:92]
a[a > 0] <- 1
scores[,91:92] <- a

range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}
data.x <- apply(scores, 2, range01)
annos <- data.frame(class = types)
rownames(annos) <- rownames(data.x)
colors <- colorRampPalette(c("blue", "white", "red"))(100)
png(filename = paste(args[1], "test.pheatmap.png", sep = "."), width = 2500, height = 1200)
myplot <- pheatmap(data.x, scale = "none", annotation_row = annos,
	show_rownames = F, show_colnames = T, color = colors, 
	cluster_cols = F, cluster_rows = F)
dev.off()
pdf(file = paste(args[1], "test.pheatmap.pdf", sep = "."), width = 25, height = 12)
myplot <- pheatmap(data.x, scale = "none", annotation_row = annos,
	show_rownames = F, show_colnames = T, color = colors, 
	cluster_cols = F, cluster_rows = F)
dev.off()




data.y <- data.x[1:5286,c(3,6,17,44,52:53,63,85,91:92)]
png(filename = paste(args[1], "sub.pheatmap.png", sep = "."), width = 1000, height = 1200)
myplot <- pheatmap(data.y, scale = "none", annotation_row = annos,
	show_rownames = F, show_colnames = T, color = colors, 
	cluster_cols = F, cluster_rows = F)
dev.off()
pdf(file = paste(args[1], "sub.pheatmap.pdf", sep = "."), width = 10, height = 12)
myplot <- pheatmap(data.y, scale = "none", annotation_row = annos,
	show_rownames = F, show_colnames = T, color = colors, 
	cluster_cols = F, cluster_rows = F)
dev.off()

data.z <- data.x[1:2551,c(3,6,17,44,52:53,63,85,91:92)]
png(filename = paste(args[1], "hkg.sub.pheatmap.png", sep = "."), width = 1000, height = 1200)
myplot <- pheatmap(data.z, scale = "none", annotation_row = annos,
	show_rownames = F, show_colnames = T, color = colors, 
	cluster_cols = F, cluster_rows = F)
dev.off()




subs <- rbind.data.frame(data.x[scores[,91]==0 & scores[,92]==0 & class=="hkg1",],
	data.x[scores[,91]==0 & scores[,92]==0 & class=="hkg2",],
	data.x[scores[,91]==0 & scores[,92]==0 & class=="hkg3",],
	data.x[scores[,91]==0 & scores[,92]==0 & class=="hkg4",]
	)
png(filename = paste(args[1], "hkgnotinTAD.sub.pheatmap.png", sep = "."), width = 1000, height = 1200)
myplot <- pheatmap(subs, scale = "none", annotation_row = annos,
	show_rownames = T, show_colnames = T, color = colors, 
	cluster_cols = F, cluster_rows = F)
dev.off()
subs2 <- rbind.data.frame(data.x[scores[,6]==0 & scores[,91]==0 & scores[,92]==0 & class=="hkg1",],
	data.x[scores[,6]==0 & scores[,6]==0 & scores[,91]==0 & scores[,92]==0 & class=="hkg2",],
	data.x[scores[,6]==0 & scores[,91]==0 & scores[,92]==0 & class=="hkg3",],
	data.x[scores[,6]==0 & scores[,91]==0 & scores[,92]==0 & class=="hkg4",]
	)
png(filename = paste(args[1], "hkgnotinTADnotinGMloop.sub.pheatmap.png", sep = "."), width = 1000, height = 1200)
myplot <- pheatmap(subs2, scale = "none", annotation_row = annos,
	show_rownames = T, show_colnames = T, color = colors, 
	cluster_cols = F, cluster_rows = F)
dev.off()


scores[scores[,91]==0 & scores[,92]==0 & class=="hkg1",c(1,3,6,91:92)]
scores[scores[,91]==0 & scores[,92]==0 & class=="hkg2",c(1,3,6,91:92)]
scores[scores[,91]==0 & scores[,92]==0 & class=="hkg3",c(1,3,6,91:92)]
scores[scores[,91]==0 & scores[,92]==0 & class=="hkg4",c(1,3,6,91:92)]






library(data.table)
library(pheatmap)
args <- c("hkg.tsg.srtbyPCA.class.enhanceratlas")
input <- fread(args[1], sep = "\t", header = T)
### for new data.table
scores <- data.matrix(input[,-c(1:2)])
genes <- as.matrix(input[,1])
class <- as.matrix(input[,2])
types <- factor(class)

enhanceratlas01 <- scores
enhanceratlas01[enhanceratlas01 > 0] <- 1
annos <- data.frame(class = types)
rownames(annos) <- rownames(enhanceratlas01)
colors <- colorRampPalette(c("blue", "white", "red"))(100)
png(filename = paste(args[1], "test.pheatmap.png", sep = "."), width = 2500, height = 1200)
myplot <- pheatmap(enhanceratlas01, scale = "none", 
	show_rownames = F, show_colnames = T, color = colors, 
	cluster_cols = F, cluster_rows = F)
dev.off()
pdf(file = paste(args[1], "test.pheatmap.pdf", sep = "."), width = 25, height = 12)
myplot <- pheatmap(enhanceratlas01, scale = "none", 
	show_rownames = F, show_colnames = T, color = colors, 
	cluster_cols = F, cluster_rows = F)
dev.off()
