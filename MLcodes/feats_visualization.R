# refine
.libPaths("/work/04935/shaojf/stampede2/myTools/anaconda2/lib/R/library")
library(data.table)
library(ggplot2)
library(pheatmap)
library(Rtsne)
set.seed(123)
args <- c("hkg.tsg.srtbyPCA.all2")
input <- fread(args[1], sep = "\t", header = T)
scores <- data.matrix(input[,-c(1:2)])
genes <- as.matrix(input[,1])
class <- as.matrix(input[,2])
types <- factor(class)
rownames(scores) <- genes
types.sim <- as.vector(types)
types.sim[1:2551] <- "HKG"
types.sim[2552:5286] <- "singleTSG"
colnames(scores)[3574:3644] <- gsub("-","_",colnames(scores)[3574:3644])
scores.x <- scores

myscale <- function(x){
	th <- as.vector(quantile(x, probs = 0.95))
	x[x > th] <- th
	return(x)
}
range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}

RuntSNE <- function(x, pre = args[1]){
	data <- x
	pre <- pre
	tsne <- Rtsne(data, check_duplicates = FALSE, dims = 2, perplexity = 30, verbose = TRUE, max_iter = 500)
	v.tsne <- as.data.frame(tsne$Y)
	colnames(v.tsne) <- c("tSNE.1", "tSNE.2")
	write.table(v.tsne, paste(pre, "tSNE.tsv", sep = "."), row.names = T, sep = "\t")
	return(tsne)
}
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

bigdiffs <- read.table("bigdiff.features.tsv", sep = "\t", row.names = 1)
feats1 <- read.table("important.features.tsv", sep = "\t", row.names = 1)
feats2 <- read.table("sim.important.features.tsv", sep = "\t", row.names = 1)
interfeats <- intersect(intersect(rownames(feats1),rownames(feats2)),rownames(bigdiffs))
interfeats.1 <- c(colnames(scores)[c(2,6,10,12,3538:3539)],interfeats)
data <- scores[,colnames(scores) %in% interfeats.1]
data <- apply(data, 2, myscale)
data <- apply(data, 2, range01)

IEratio <- scores.x[,4]/scores.x[,3]
IEratio.scale <- IEratio
IEratio.scale[IEratio.scale > 20] <- 20
IEratio.scale <- range01(IEratio)
data.x <- data.frame(IEratio.scale, data)

annos <- data.frame(class = types, class.sim = types.sim)
rownames(annos) <- rownames(data.x)
colors <- colorRampPalette(c("white", "blue"))(100)
png(filename = paste(args[1], "intersection.pheatmap.png", sep = "."), width = 2500, height = 1200)
myplot <- pheatmap(data.x, scale = "none", annotation_row = annos,
	show_rownames = F, show_colnames = T, color = colors, # gaps_row = c(2552, 5287),
	cluster_cols = F, cluster_rows = F)
dev.off()


mydata.tsne <- RuntSNE(x = data.x, pre = paste(args[1],"interfeats.1", sep = "."))
mydata.plot <- PlottSNE(x = mydata.tsne, pre = paste(args[1],"interfeats.1", sep = "."))

mydata.hkg.tsne <- RuntSNE(x = data.x[types.sim=="HKG",], pre = paste(args[1],"interfeats.1.hkg", sep = "."))
mydata.hkg.plot <- PlottSNE(x = mydata.hkg.tsne, y = types[types.sim=="HKG"], z = types.sim[types.sim=="HKG"], pre = paste(args[1],"interfeats.1.hkg", sep = "."))
png(filename = paste(args[1], "intersection.1.hkg.pheatmap.png", sep = "."), width = 2000, height = 2500)
myplot <- pheatmap(data.x[types.sim=="HKG",], scale = "none", annotation_row = annos[types.sim=="HKG",],
	show_rownames = T, show_colnames = T, color = colors, # gaps_row = c(2552, 5287),
	cluster_cols = T, cluster_rows = T)
dev.off()



data <- scores[,colnames(scores) %in% interfeats]
data <- apply(data, 2, myscale)
data <- apply(data, 2, range01)
data.tsne <- RuntSNE(x = data, pre = paste(args[1],"interfeats", sep = "."))
data.plot <- PlottSNE(x = mydata.tsne, pre = paste(args[1],"interfeats", sep = "."))
png(filename = paste(args[1], "intersection.pheatmap.png", sep = "."), width = 2500, height = 1200)
myplot <- pheatmap(data, scale = "none", annotation_row = annos,
	show_rownames = F, show_colnames = T, color = colors, # gaps_row = c(2552, 5287),
	cluster_cols = F, cluster_rows = F)
dev.off()
data.hkg.tsne <- RuntSNE(x = data[types.sim=="HKG",], pre = paste(args[1],"interfeats.hkg", sep = "."))
data.hkg.plot <- PlottSNE(x = data.hkg.tsne, y = types[types.sim=="HKG"], z = types.sim[types.sim=="HKG"], pre = paste(args[1],"interfeats.hkg", sep = "."))
png(filename = paste(args[1], "intersection.hkg.pheatmap.png", sep = "."), width = 2000, height = 2500)
myplot <- pheatmap(data[types.sim=="HKG",], scale = "none", annotation_row = annos[types.sim=="HKG",],
	show_rownames = T, show_colnames = T, color = colors, # gaps_row = c(2552, 5287),
	cluster_cols = T, cluster_rows = T)
dev.off()


ets.col <- colnames(scores.x)[grep("ETS",colnames(scores.x))]
a <- colnames(scores.x)[grep("GTRD.EHF",colnames(scores.x))]
b <- colnames(scores.x)[grep("GTRD.ERG",colnames(scores.x))]
c <- colnames(scores.x)[grep("GTRD.ETV",colnames(scores.x))]
d <- colnames(scores.x)[grep("GTRD.Elf",colnames(scores.x))]
e <- colnames(scores.x)[grep("GTRD.Elk",colnames(scores.x))]
f <- colnames(scores.x)[grep("GTRD.Fli",colnames(scores.x))]
g <- colnames(scores.x)[grep("GTRD.GABP",colnames(scores.x))]
h <- colnames(scores.x)[grep("GTRD.Spi",colnames(scores.x))]
i <- colnames(scores.x)[grep("GTRD.PU.1",colnames(scores.x))]
ets.col <- c(ets.col, a, b, c, d, e, f, g, h, i)
yy1.col <- colnames(scores.x)[grep("YY1",colnames(scores.x))]
sp1.col <- colnames(scores.x)[grep("Sp1",colnames(scores.x))]
loop.col <- colnames(scores.x)[grep("HiC",colnames(scores.x))]
selects <- c(colnames(scores.x)[c(6,10,12)],ets.col, yy1.col, sp1.col, loop.col)
data.s <- scores[,colnames(scores) %in% selects]
data.s <- apply(data.s, 2, myscale)
data.s <- apply(data.s, 2, range01)

png(filename = paste(args[1], "selects.pheatmap.png", sep = "."), width = 2500, height = 1500)
myplot <- pheatmap(data.s, scale = "none", annotation_row = annos,
	show_rownames = F, show_colnames = T, color = colors, # gaps_row = c(2552, 5287),
	cluster_cols = F, cluster_rows = F)
dev.off()
png(filename = paste(args[1], "selects.hkg.pheatmap.png", sep = "."), width = 2500, height = 3000)
myplot <- pheatmap(t(na.omit(t(data.s[types.sim=="HKG",]))), scale = "none", annotation_row = annos[types.sim=="HKG",],
	show_rownames = T, show_colnames = T, color = colors, fontsize_row = 3, # gaps_row = c(2552, 5287),
	cluster_cols = T, cluster_rows = T)
dev.off()
used <- t(na.omit(t(data.s[types.sim=="HKG",])))
used.ordered <- used[myplot$tree_row$order, myplot$tree_col$order]
used.tsne <- RuntSNE(x = used.ordered, pre = paste(args[1],"selects.hkg.1.hkg", sep = "."))
used.tsne.plot <- PlottSNE(x = used.tsne, y = types[types.sim=="HKG"], z = types.sim[types.sim=="HKG"], pre = paste(args[1],"selects.hkg.1.hkg", sep = "."))
out <- data.frame(rownames(used.ordered), used.ordered, used.tsne$Y)
write.table(out, paste(args[1], "selects.hkg.featureswithtSNE.tsv", sep = "."), sep = "\t", row.names = F)

tsne.scale <- apply(used.tsne$Y, 2, range01)
data.v <- data.frame(used.ordered,tsne.scale)
colnames(data.v)[206:207] <- c("tSNE.1", "tSNE.2")
png(filename = paste(args[1], "selects.hkg.tSNE.pheatmap.png", sep = "."), width = 2500, height = 1500)
pheatmap(data.v, scale = "none", annotation_row = annos[types.sim=="HKG",],
	show_rownames = F, show_colnames = T, color = colors, # gaps_row = c(2552, 5287),
	cluster_cols = F, cluster_rows = F)
dev.off()



#############################################
# FS
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
	# save.image(paste(pre, "RData", sep = ".")
	write.table(feat$all_feat$"glmnet-lasso", paste(pre, "FeatureSelection.glmnetlasso.tsv", sep = "."), row.names = F, sep = "\t")
	write.table(feat$all_feat$xgboost, paste(pre, "FeatureSelection.xgboost.tsv", sep = "."), row.names = F, sep = "\t")
	write.table(feat$all_feat$ranger, paste(pre, "FeatureSelection.ranger.tsv", sep = "."), row.names = F, sep = "\t")
	write.table(feat$union_feat, paste(pre, "FeatureSelection.union_feat.tsv", sep = "."), row.names = F, sep = "\t")
	return(feat)
}
data.s.na <- t(na.omit(t(data.s)))
data.s.feat <- RunFS(x = data.s.na, y = as.numeric(types), pre = paste(args[1],"selects", sep = "."))
#############################################



# data.s[is.na(data.s)] <- 0
data.sx <- data.s
data.s <- t(na.omit(t(data.s)))
data.s.tsne <- RuntSNE(x = data.s, pre = paste(args[1],"data.s", sep = "."))
data.s.plot <- PlottSNE(x = data.s.tsne, y = types, z = types.sim, pre = paste(args[1],"data.s", sep = "."))
data.s.hkg.tsne <- RuntSNE(x = data.s[types.sim=="HKG",], pre = paste(args[1],"data.s.hkg", sep = "."))
data.s.hkg.plot <- PlottSNE(x = data.s.hkg.tsne, y = types[types.sim=="HKG"], z = types.sim[types.sim=="HKG"], pre = paste(args[1],"data.s.hkg", sep = "."))

save.image(paste(args[1], "mytest.RData", sep = "."))
