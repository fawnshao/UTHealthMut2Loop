# .libPaths("/work/04935/shaojf/stampede2/myTools/anaconda2/lib/R/library")
library(data.table)
library(ggplot2)
library(pheatmap)
# args <- commandArgs(TRUE)
# args <- c("hkg.tsg.srtbyPCA.transcript.class.vert.known.motifs.mat")
args <- c("hkg.tsg.srtbyPCA.transcript.class.top30.vert.known.motifs.mat")
input <- fread(args[1], sep = "\t", header = T)
### for old data.table
# scores <- data.matrix(input[,.SD, .SDcols = c(3:170)])
# genes <- as.matrix(input[, .SD, .SDcols = 1])
# class <- as.matrix(input[, .SD, .SDcols = 2])
### for new data.table
scores <- data.matrix(input[,-c(1:3)])
genes <- as.matrix(input[,1])
class <- as.matrix(input[,2])
trans <- as.matrix(input[,3])
types <- factor(class)
rownames(scores) <- trans
types.sim <- as.vector(types)
types.sim[c(3663:5094,6050:6084)] <- "othersingleTSG"

annos <- data.frame(class = types)
rownames(annos) <- trans
colors <- colorRampPalette(c("blue", "white", "red"))(100)

Plotheatmap <- function(x, y = annos, z = colors, pre = args[1]){
	# png(filename = paste(pre, "pheatmap.png", sep = "."), width = 1000, height = 1600)
	pdf(file = paste(pre, "pheatmap.pdf", sep = "."), width = 12, height = 12)
	myplot <- pheatmap(x, scale = "none", annotation_row = y,
		show_rownames = F, show_colnames = T, color = z, gaps_row = c(3662, 6084),
		cluster_cols = F, cluster_rows = F)
	dev.off()
}
range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}

withmotif <- apply(scores, 2, function(x){length(x[x>0])})/nrow(scores)
withmotif.1 <- apply(scores[1:3662,], 2, function(x){length(x[x>0])})/nrow(scores[1:3662,])
withmotif.2 <- apply(scores[3663:6084,], 2, function(x){length(x[x>0])})/nrow(scores[3663:6084,])
withmotif.3 <- apply(scores[6085:7796,], 2, function(x){length(x[x>0])})/nrow(scores[6085:7796,])
withmotif[order(withmotif, decreasing=T)][1:10]
withmotif.1 <- withmotif.1[order(withmotif.1, decreasing=T)]
withmotif.2 <- withmotif.2[order(withmotif.2, decreasing=T)]
withmotif.3 <- withmotif.3[order(withmotif.3, decreasing=T)]

out <- data.frame(names(withmotif.1),withmotif.1,names(withmotif.2),withmotif.2,names(withmotif.3),withmotif.3)
write.table(out, file = "hkg.tsg.srtbyPCA.transcript.300up100down.homer.vert.known.tsv", sep ="\t", quote =T, row.names = F, col.names = T)

data.x <- apply(scores, 2, range01)
Plotheatmap(data.x, pre = args[1])







range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}
seqfeature <- scores[,1:6]
seqfeature01 <- apply(seqfeature, 2, range01)
homermotif <- scores[,7:369]
gtrdpeaks <- scores[,370:2539]
meth <- scores[,2540:2576]
histone <- scores[,2577:3462]
dnase <- scores[,3463:3501]
hic <- scores[,3502:3514]
pchic <- scores[,3515:3548]

save.image(paste(args[1], "detail.data.RData", sep = "."))

###################################################
#####         look at the following           #####
###################################################
## add hkg classification by PCA
############## FeatureSelection
args <- c("onlyhkg.tsg.further.srt.all")
load(paste(args[1], "detail.data.RData", sep = "."))
RunFS <- function(x, y = as.numeric(types), pre = args[1]){
	require(FeatureSelection)
	trainx <- x
	trainy <- y
	pre <- pre
	params_glmnet <- list(alpha = 1, family = 'gaussian', nfolds = 5, parallel = TRUE)
	params_xgboost <- list( params = list("objective" = "reg:linear", "bst:eta" = 0.001, "subsample" = 0.75, "max_depth" = 5, 
	                                     "colsample_bytree" = 0.75, "nthread" = 60),
	                                      nrounds = 1000, print.every.n = 250, maximize = FALSE)
	params_ranger <- list(dependent.variable.name = 'y', probability = FALSE, num.trees = 1000, verbose = TRUE, mtry = 5, 
	                     min.node.size = 10, num.threads = 60, classification = FALSE, importance = 'permutation')
	params_features <- list(keep_number_feat = NULL, union = TRUE)
	feat <- wrapper_feat_select(X = trainx, y = trainy, params_glmnet = params_glmnet, params_xgboost = params_xgboost, 
	                          params_ranger = params_ranger, xgb_sort = 'Gain', CV_folds = 5, stratified_regr = FALSE, 
	                          scale_coefs_glmnet = FALSE, cores_glmnet = 5, params_features = params_features, verbose = TRUE)
	# save.image(paste(pre, "RData", sep = ".")
	write.table(feat$all_feat$"glmnet-lasso", paste(pre, "FeatureSelection.glmnetlasso.tsv", sep = "."), row.names = F, sep = "\t")
	write.table(feat$all_feat$xgboost, paste(pre, "FeatureSelection.xgboost.tsv", sep = "."), row.names = F, sep = "\t")
	write.table(feat$all_feat$ranger, paste(pre, "FeatureSelection.ranger.tsv", sep = "."), row.names = F, sep = "\t")
	write.table(feat$union_feat, paste(pre, "FeatureSelection.union_feat.tsv", sep = "."), row.names = F, sep = "\t")
	return(feat)
}
seqfeature.feat <- RunFS(x = seqfeature, y = as.numeric(types), pre = paste(args[1],"seqfeature.raw", sep = "."))
seqfeature01.feat <- RunFS(x = seqfeature01, y = as.numeric(types), pre = paste(args[1],"seqfeature01", sep = "."))
homermotif.feat <- RunFS(x = homermotif, y = as.numeric(types), pre = paste(args[1],"homermotif", sep = "."))
gtrdpeaks.feat <- RunFS(x = gtrdpeaks, y = as.numeric(types), pre = paste(args[1],"gtrdpeaks", sep = "."))
meth.feat <- RunFS(x = meth, y = as.numeric(types), pre = paste(args[1],"meth", sep = "."))
histone.feat <- RunFS(x = histone, y = as.numeric(types), pre = paste(args[1],"histone", sep = "."))
dnase.feat <- RunFS(x = dnase, y = as.numeric(types), pre = paste(args[1],"dnase", sep = "."))
hic.feat <- RunFS(x = hic, y = as.numeric(types), pre = paste(args[1],"hic", sep = "."))
pchic.feat <- RunFS(x = pchic, y = as.numeric(types), pre = paste(args[1],"pchic", sep = "."))
save.image(paste(args[1], "FeatureSelection.RData", sep = ".")
# cat onlyhkg.tsg.further.srt.all.*.FeatureSelection.union_feat.tsv | awk '$2>0.5'
# cat onlyhkg.tsg.further.srt.all.dnase.FeatureSelection.union_feat.tsv onlyhkg.tsg.further.srt.all.pchic.FeatureSelection.union_feat.tsv onlyhkg.tsg.further.srt.all.hic.FeatureSelection.union_feat.tsv onlyhkg.tsg.further.srt.all.seqfeature.FeatureSelection.union_feat.tsv onlyhkg.tsg.further.srt.all.homermotif.FeatureSelection.union_feat.tsv onlyhkg.tsg.further.srt.all.meth.FeatureSelection.union_feat.tsv onlyhkg.tsg.further.srt.all.gtrdpeaks.FeatureSelection.union_feat.tsv onlyhkg.tsg.further.srt.all.histone.FeatureSelection.union_feat.tsv | awk '$2>0.5' > important.features.tsv
############## FeatureSelection use simplyfy classification
args <- c("onlyhkg.tsg.further.srt.all")
load(paste(args[1], "detail.data.RData", sep = "."))
types.sim <- as.vector(types)
types.sim[1:2506] <- "hkg"
types.sim[c(2507:4263)] <- "singleTSG"
RunFS <- function(x, y = as.numeric(factor(types.sim)), pre = args[1]){
	require(FeatureSelection)
	trainx <- x
	trainy <- y
	pre <- pre
	params_glmnet <- list(alpha = 1, family = 'gaussian', nfolds = 5, parallel = TRUE)
	params_xgboost <- list( params = list("objective" = "reg:linear", "bst:eta" = 0.001, "subsample" = 0.75, "max_depth" = 5, 
	                                     "colsample_bytree" = 0.75, "nthread" = 60),
	                                      nrounds = 1000, print.every.n = 250, maximize = FALSE)
	params_ranger <- list(dependent.variable.name = 'y', probability = FALSE, num.trees = 1000, verbose = TRUE, mtry = 5, 
	                     min.node.size = 10, num.threads = 60, classification = FALSE, importance = 'permutation')
	params_features <- list(keep_number_feat = NULL, union = TRUE)
	feat <- wrapper_feat_select(X = trainx, y = trainy, params_glmnet = params_glmnet, params_xgboost = params_xgboost, 
	                          params_ranger = params_ranger, xgb_sort = 'Gain', CV_folds = 5, stratified_regr = FALSE, 
	                          scale_coefs_glmnet = FALSE, cores_glmnet = 5, params_features = params_features, verbose = TRUE)
	# save.image(paste(pre, "RData", sep = ".")
	write.table(feat$all_feat$"glmnet-lasso", paste(pre, "FeatureSelection.glmnetlasso.tsv", sep = "."), row.names = F, sep = "\t")
	write.table(feat$all_feat$xgboost, paste(pre, "FeatureSelection.xgboost.tsv", sep = "."), row.names = F, sep = "\t")
	write.table(feat$all_feat$ranger, paste(pre, "FeatureSelection.ranger.tsv", sep = "."), row.names = F, sep = "\t")
	write.table(feat$union_feat, paste(pre, "FeatureSelection.union_feat.tsv", sep = "."), row.names = F, sep = "\t")
	return(feat)
}
seqfeature.feat <- RunFS(x = seqfeature, pre = paste(args[1],"seqfeature.raw.sim", sep = "."))
seqfeature01.feat <- RunFS(x = seqfeature01, pre = paste(args[1],"seqfeature01.sim", sep = "."))
homermotif.feat <- RunFS(x = homermotif, pre = paste(args[1],"homermotif.sim", sep = "."))
gtrdpeaks.feat <- RunFS(x = gtrdpeaks, pre = paste(args[1],"gtrdpeaks.sim", sep = "."))
meth.feat <- RunFS(x = meth, pre = paste(args[1],"meth.sim", sep = "."))
histone.feat <- RunFS(x = histone, pre = paste(args[1],"histone.sim", sep = "."))
dnase.feat <- RunFS(x = dnase, pre = paste(args[1],"dnase.sim", sep = "."))
hic.feat <- RunFS(x = hic, pre = paste(args[1],"hic.sim", sep = "."))
pchic.feat <- RunFS(x = pchic, pre = paste(args[1],"pchic.sim", sep = "."))
save.image(paste(args[1], "FeatureSelection.sim.RData", sep = ".")

############## PCA
args <- c("onlyhkg.tsg.further.srt.all")
load(paste(args[1], "detail.data.RData", sep = "."))
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
seqfeature.pca <- RunPCA(x = seqfeature, pre = paste(args[1],"seqfeature.raw", sep = "."))
seqfeature01.pca <- RunPCA(x = seqfeature01, pre = paste(args[1],"seqfeature01", sep = "."))
homermotif.pca <- RunPCA(x = homermotif, pre = paste(args[1],"homermotif", sep = "."))
gtrdpeaks.pca <- RunPCA(x = gtrdpeaks, pre = paste(args[1],"gtrdpeaks", sep = "."))
meth.pca <- RunPCA(x = meth, pre = paste(args[1],"meth", sep = "."))
histone.pca <- RunPCA(x = histone, pre = paste(args[1],"histone", sep = "."))
dnase.pca <- RunPCA(x = dnase, pre = paste(args[1],"dnase", sep = "."))
hic.pca <- RunPCA(x = hic, pre = paste(args[1],"hic", sep = "."))
pchic.pca <- RunPCA(x = pchic, pre = paste(args[1],"pchic", sep = "."))
save.image(paste(args[1], "PCA.RData", sep = "."))

############## t-SNE
# conda.R
args <- c("onlyhkg.tsg.further.srt.all")
load(paste(args[1], "detail.data.RData", sep = "."))
.libPaths("/work/04935/shaojf/stampede2/myTools/anaconda2/lib/R/library")
library(ggplot2)
library(Rtsne)
set.seed(123)
RuntSNE <- function(x, y = types, z = types.sim, pre = args[1]){
	# require(ggplot2)
	data <- x
	# mycolors <- y
	# mycolors.sim <- z
	pre <- pre
	# a <- floor(ncol(data) / 10)
	# tsne.t <- Rtsne(t(data), check_duplicates = FALSE, dims = 2, perplexity = a, verbose = TRUE, max_iter = 500)
	tsne <- Rtsne(data, check_duplicates = FALSE, dims = 2, perplexity = 30, verbose = TRUE, max_iter = 500)
	v.tsne <- as.data.frame(tsne$Y)
	colnames(v.tsne) <- c("tSNE.1", "tSNE.2")
	# pdf(file = paste(pre, "tsne.pdf", sep = "."), width = 15, height = 12)
	# p1 <- ggplot(data = v.tsne, aes(x = tSNE.1, y = tSNE.2)) + 
	# 	geom_point(aes(colour = y, alpha = 1/10)) +
	# 	ggtitle(pre)
	# print(p1)
	# dev.off()
	# pdf(file = paste(pre, "tsne.sim.pdf", sep = "."), width = 15, height = 12)
	# p1 <- ggplot(data = v.tsne, aes(x = tSNE.1, y = tSNE.2)) + 
	# 	geom_point(aes(colour = z, alpha = 1/10)) +
	# 	ggtitle(pre)
	# print(p1)
	# dev.off()
	write.table(v.tsne, paste(pre, "tSNE.tsv", sep = "."), row.names = T, sep = "\t")

	# v.tsne <- as.data.frame(tsne.t$Y)
	# colnames(v.tsne) <- c("tSNE.1", "tSNE.2")
	# pdf(file = paste(pre, "tsne.t.pdf", sep = "."), width = 15, height = 12)
	# p1 <- ggplot(data = v.tsne, aes(x = tSNE.1, y = tSNE.2)) + 
	# 	geom_point() +
	# 	ggtitle(pre)
	# print(p1)
	# dev.off()
	# write.table(v.tsne, paste(pre, "tSNE.t.tsv", sep = "."), row.names = T, sep = "\t")
	return(tsne)
}
seqfeature.tsne <- RuntSNE(x = seqfeature, pre = paste(args[1],"seqfeature.raw", sep = "."))
seqfeature01.tsne <- RuntSNE(x = seqfeature01, pre = paste(args[1],"seqfeature01", sep = "."))
homermotif.tsne <- RuntSNE(x = homermotif, pre = paste(args[1],"homermotif", sep = "."))
gtrdpeaks.tsne <- RuntSNE(x = gtrdpeaks, pre = paste(args[1],"gtrdpeaks", sep = "."))
meth.tsne <- RuntSNE(x = meth, pre = paste(args[1],"meth", sep = "."))
histone.tsne <- RuntSNE(x = histone, pre = paste(args[1],"histone", sep = "."))
dnase.tsne <- RuntSNE(x = dnase, pre = paste(args[1],"dnase", sep = "."))
hic.tsne <- RuntSNE(x = hic, pre = paste(args[1],"hic", sep = "."))
pchic.tsne <- RuntSNE(x = pchic, pre = paste(args[1],"pchic", sep = "."))
save.image(paste(args[1], "detail.tSNE.RData", sep = "."))
#### plot t-SNE
library(ggplot2)
args <- c("onlyhkg.tsg.further.srt.all")
load(paste(args[1], "detail.tSNE.RData", sep = "."))
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
seqfeature.tsne.plot <- PlottSNE(x = seqfeature.tsne, pre = paste(args[1],"seqfeature.raw", sep = "."))
seqfeature01.tsne.plot <- PlottSNE(x = seqfeature01.tsne, pre = paste(args[1],"seqfeature01", sep = "."))
homermotif.tsne.plot <- PlottSNE(x = homermotif.tsne, pre = paste(args[1],"homermotif", sep = "."))
gtrdpeaks.tsne.plot <- PlottSNE(x = gtrdpeaks.tsne.tsne, pre = paste(args[1],"gtrdpeaks", sep = "."))
meth.tsne.plot <- PlottSNE(x = meth.tsne, pre = paste(args[1],"meth", sep = "."))
histone.tsne.plot <- PlottSNE(x = histone.tsne, pre = paste(args[1],"histone", sep = "."))
dnase.tsne.plot <- PlottSNE(x = dnase.tsne, pre = paste(args[1],"dnase", sep = "."))
hic.tsne.plot <- PlottSNE(x = hic.tsne, pre = paste(args[1],"hic", sep = "."))
pchic.tsne.plot <- PlottSNE(x = pchic.tsne, pre = paste(args[1],"pchic", sep = "."))

# please also ref to: 
# http://r-statistics.co/Variable-Selection-and-Importance-With-R.html
###################################################
#####           look at the above             #####
###################################################

############## visualization of the important features
# cat onlyhkg.tsg.further.srt.all.dnase.FeatureSelection.union_feat.tsv onlyhkg.tsg.further.srt.all.pchic.FeatureSelection.union_feat.tsv onlyhkg.tsg.further.srt.all.hic.FeatureSelection.union_feat.tsv onlyhkg.tsg.further.srt.all.seqfeature.FeatureSelection.union_feat.tsv onlyhkg.tsg.further.srt.all.homermotif.FeatureSelection.union_feat.tsv onlyhkg.tsg.further.srt.all.meth.FeatureSelection.union_feat.tsv onlyhkg.tsg.further.srt.all.gtrdpeaks.FeatureSelection.union_feat.tsv onlyhkg.tsg.further.srt.all.histone.FeatureSelection.union_feat.tsv | awk '$2>0.8' > important.features.tsv
# cat onlyhkg.tsg.further.srt.all.dnase.sim.FeatureSelection.union_feat.tsv onlyhkg.tsg.further.srt.all.pchic.sim.FeatureSelection.union_feat.tsv onlyhkg.tsg.further.srt.all.hic.sim.FeatureSelection.union_feat.tsv onlyhkg.tsg.further.srt.all.seqfeature01.sim.FeatureSelection.union_feat.tsv onlyhkg.tsg.further.srt.all.homermotif.sim.FeatureSelection.union_feat.tsv onlyhkg.tsg.further.srt.all.meth.sim.FeatureSelection.union_feat.tsv onlyhkg.tsg.further.srt.all.gtrdpeaks.sim.FeatureSelection.union_feat.tsv onlyhkg.tsg.further.srt.all.histone.sim.FeatureSelection.union_feat.tsv | awk '$2>0.8' > sim.important.features.tsv
args <- c("onlyhkg.tsg.further.srt.all")
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


