# .libPaths("/work/04935/shaojf/stampede2/myTools/anaconda2/lib/R/library")
library(data.table)
library(ggplot2)
library(pheatmap)
library(cluster)
library(randomForest)
library(glmnet)
# install.packages("Rtsne")
# library(Rtsne)
# devtools::install_github('mlampros/FeatureSelection')
# library(FeatureSelection)
args <- commandArgs(TRUE)
# args <- c("onlyhkg.tsg.all")
input <- fread(args[1], sep = "\t", header = T)
### for old data.table
# scores <- data.matrix(input[,.SD, .SDcols = c(3:170)])
# genes <- as.matrix(input[, .SD, .SDcols = 1])
# class <- as.matrix(input[, .SD, .SDcols = 2])
### for new data.table
scores <- data.matrix(input[,-c(1:2)])
genes <- as.matrix(input[,1])
class <- as.matrix(input[,2])
types <- factor(class)
rownames(scores) <- genes
types.sim <- as.vector(types)
types.sim[types.sim!="hkg" & types.sim!="mixTSG" & types.sim!="Testis"] <- "othersingleTSG"

###################################################
#####         look at the following           #####
###################################################
############## prepare the data
range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}
scaled.scores <- apply(scores, 2, range01)
save.image("onlyhkg.tsg.data.RData")
############## PCA
load("onlyhkg.tsg.data.RData")
library(ggplot2)
scaled.scores[is.na(scaled.scores)] <- 0
myPCA.a <- prcomp(scaled.scores, scale. = F, center = F)
myPCA.a.t <- prcomp(t(scaled.scores), scale. = F, center = F)
pca.scores <- as.data.frame(myPCA.a$x)
pdf(file = paste(args[1], "range01.pca.pdf", sep = "."), width = 15, height = 12)
ggplot(data = pca.scores, aes(x = PC1, y = PC2)) +
	geom_point(aes(colour = types, alpha = 1/10)) +
	ggtitle(args[1])
dev.off()
pdf(file = paste(args[1], "range01.pca.sim.pdf", sep = "."), width = 15, height = 12)
ggplot(data = pca.scores, aes(x = PC1, y = PC2)) +
	geom_point(aes(colour = types.sim, alpha = 1/10)) +
	ggtitle(args[1])
dev.off()

pca.scores <- as.data.frame(myPCA.a.t$x)
pdf(file = paste(args[1], "range01.t.pca.pdf", sep = "."), width = 15, height = 12)
ggplot(data = pca.scores, aes(x = PC1, y = PC2, label = rownames(pca.scores))) +
	geom_point(aes(alpha = 1/10)) +
	geom_text(colour = "tomato", alpha = 0.8, size = 2) +
	ggtitle(args[1])
dev.off()
save.image("onlyhkg.tsg.PCA.RData")
## feature selection
write.csv(pca.scores, paste(args[1], "range01.t.pca.csv", sep = "."))
feature.select <- order(abs(pca.scores[,1]), decreasing = T)[1:500]
colnames(scores)[feature.select]
pca.scores[feature.select,1]
colors <- colorRampPalette(c("blue", "white", "red"))(10)
data <- scaled.scores[,feature.select]
library(pheatmap)
png(filename = paste(args[1], "range01.top500.pheatmap.png", sep = "."), width = 1500, height = 1200)
myplot <- pheatmap(data, scale = "none", 
	show_rownames = F, show_colnames = T, color = colors, 
	cluster_cols = F, cluster_rows = F)
dev.off()
############## LASSO
# load("all.PCA.RData")
# lmdata <- data.frame(subsets, as.numeric(subsets.types))
# my.lm <- lm(subsets.types ~ ., data = lmdata)
# anova(my.lm)
# summary(my.lm)
my.lasso <- glmnet(subsets, as.numeric(subsets.types), standardize = TRUE, alpha = 1)
# a <- coef(my.lasso, s = 0.1)
png(filename = paste(args[1], "lasso.lambda.png", sep = "."), width = 1500, height = 1200)
plot(my.lasso, xvar = "lambda", label = TRUE)
dev.off()
# my.lasso.cv <- cv.glmnet(subsets, as.numeric(subsets.types), standardize = TRUE, alpha = 1, type.measure = "mse", nfolds = 5)
save.image("all.LASSO.RData")
############## FeatureSelection
# install_github('mlampros/FeatureSelection') 
# load("all.PCA.RData")
library(FeatureSelection)
params_glmnet <- list(alpha = 1, family = 'gaussian', nfolds = 5, parallel = TRUE)
params_xgboost <- list( params = list("objective" = "reg:linear", "bst:eta" = 0.001, "subsample" = 0.75, "max_depth" = 5, 
                                     "colsample_bytree" = 0.75, "nthread" = 60),
                                      nrounds = 1000, print.every.n = 250, maximize = FALSE)
params_ranger <- list(dependent.variable.name = 'y', probability = FALSE, num.trees = 1000, verbose = TRUE, mtry = 5, 
                     min.node.size = 10, num.threads = 60, classification = FALSE, importance = 'permutation')
params_features <- list(keep_number_feat = NULL, union = TRUE)
feat <- wrapper_feat_select(X = scaled.scores, y = as.numeric(types), params_glmnet = params_glmnet, params_xgboost = params_xgboost, 
                          params_ranger = params_ranger, xgb_sort = 'Gain', CV_folds = 5, stratified_regr = FALSE, 
                          scale_coefs_glmnet = FALSE, cores_glmnet = 5, params_features = params_features, verbose = TRUE)
str(feat)
params_barplot < list(keep_features = 300, horiz = TRUE, cex.names = 1.0)
png(filename = paste(args[1], "FeatureSelection.bar.png", sep = "."), width = 1500, height = 1200)
barplot_feat_select(feat, params_barplot, xgb_sort = 'Cover')
dev.off()
dat <- data.frame(p = as.numeric(subsets.types), subsets)
cor_feat <- func_correlation(dat, target = 'p', correlation_thresh = 0.1, use_obs = 'complete.obs', correlation_method = 'pearson')
out_lst < lapply(feat$all_feat, function(x) which(rownames(cor_feat) %in% x[1:100, 1]))
cor_lasso <- func_correlation(subsets[, feat$all_feat$`glmnet-lasso`[, 1]], target = NULL, correlation_thresh = 0.9, 
                             use_obs = 'complete.obs', correlation_method = 'pearson')
cor_xgb = func_correlation(subsets[, feat$all_feat$xgboost[, 1][1:100]], target = NULL, correlation_thresh = 0.9, 
                           use_obs = 'complete.obs', correlation_method = 'pearson')
cor_rf = func_correlation(subsets[, feat$all_feat$ranger[, 1][1:100]], target = NULL, correlation_thresh = 0.9,
                          use_obs = 'complete.obs', correlation_method = 'pearson')
write.csv(cor_feat,paste(args[1], "FeatureSelection.cor_feat.csv", sep = "."))
write.csv(cor_lasso,paste(args[1], "FeatureSelection.cor_lasso.csv", sep = "."))
write.csv(cor_xgb,paste(args[1], "FeatureSelection.cor_xgb.csv", sep = "."))
write.csv(cor_rf,paste(args[1], "FeatureSelection.cor_rf.csv", sep = "."))
save.image("all.FeatureSelection.RData")
############## Caret
# load("all.PCA.RData")
# library(mlbench)
library(caret)
set.seed(123)
mydata <- data.frame(subsets, as.numeric(subsets.types))
# calculate correlation matrix
correlationMatrix <- cor(mydata)
# summarize the correlation matrix
print(correlationMatrix)
# find attributes that are highly corrected (ideally >0.75)
highlyCorrelated <- findCorrelation(correlationMatrix, cutoff = 0.5)
# print indexes of highly correlated attributes
print(highlyCorrelated)
# prepare training scheme
control <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
# train the model
# model <- train(subsets.types ~ ., data = mydata, method = "lvq", preProcess = "scale", trControl = control)
model <- train(subsets.types ~ ., data = mydata, method = "lvq", trControl = control)
# estimate variable importance
importance <- varImp(model, scale=FALSE)
# summarize importance
print(importance)
# plot importance
plot(importance)
############## MXM
# load("all.PCA.RData")
library(MXM)
set.seed(123)
mydata <- data.frame(subsets, as.numeric(subsets.types))
sesObject <- SES(target, dataset, max_k = 5, threshold = 0.2, test = "testIndFisher", hash = TRUE, hashObject = NULL))
# please also ref to: 
# http://r-statistics.co/Variable-Selection-and-Importance-With-R.html
###################################################
#####           look at the above             #####
###################################################














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
###################################################
subsets <- scores[types!="other",]
subsets.types <- factor(types[types!="other"], levels = unique(types[types!="other"]))
range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}
scaled.scores <- apply(scores, 2, range01)
### PCA
subsets.myPCA.a <- prcomp(subsets, scale. = T, center = T)
subsets.pca.scores <- as.data.frame(subsets.myPCA.a$x)
# png(filename = paste(args[1], "subsets.pca.png", sep = "."), width = 1500, height = 1200)
pdf(file = paste(args[1], "subsets.pca.pdf", sep = "."), width = 15, height = 12)
ggplot(data = subsets.pca.scores, aes(x = PC1, y = PC2)) + #, label = rownames(pca.scores))) +
	geom_point(aes(colour = types[types!="other"], alpha = 1/10)) +
	ggtitle(args[1])
dev.off()
subsets.myPCA.b <- prcomp(subsets, scale. = F, center = F)
subsets.pca.scores <- as.data.frame(subsets.myPCA.b$x)
# png(filename = paste(args[1], "subsets.noscale.pca.png", sep = "."), width = 1500, height = 1200)
pdf(file = paste(args[1], "subsets.noscale.pca.pdf", sep = "."), width = 15, height = 12)
ggplot(data = subsets.pca.scores, aes(x = PC1, y = PC2)) + #, label = rownames(pca.scores))) +
	geom_point(aes(colour = types[types!="other"], alpha = 1/10)) +
	ggtitle(args[1])
dev.off()
t.subsets.myPCA.a <- prcomp(t(subsets), scale. = T, center = T)
subsets.pca.scores <- as.data.frame(t.subsets.myPCA.a$x)
pdf(file = paste(args[1], "t.subsets.scale.pca.pdf", sep = "."), width = 15, height = 12)
ggplot(data = subsets.pca.scores, aes(x = PC1, y = PC2, label = rownames(subsets.pca.scores))) +
	geom_point(aes(alpha = 1/10)) +
	geom_text(colour = "tomato", alpha = 0.8, size = 2) +
	ggtitle(args[1])
dev.off()
t.subsets.myPCA.b <- prcomp(t(subsets), scale. = F, center = F)
subsets.pca.scores <- as.data.frame(t.subsets.myPCA.b$x)
pdf(file = paste(args[1], "t.subsets.noscale.pca.pdf", sep = "."), width = 15, height = 12)
ggplot(data = subsets.pca.scores, aes(x = PC1, y = PC2, label = rownames(subsets.pca.scores))) +
	geom_point(aes(alpha = 1/10)) +
	geom_text(colour = "tomato", alpha = 0.8, size = 2) +
	ggtitle(args[1])
dev.off()

# subsets.pca.scores <- as.data.frame(t.subsets.myPCA.a$x)
subsets <- scaled.scores[types!="other",]
t.scale.subsets.myPCA.b <- prcomp(t(subsets), scale. = F, center = F)
subsets.pca.scores <- as.data.frame(t.scale.subsets.myPCA.b$x)
feature.select <- order(abs(subsets.pca.scores[,1]), decreasing=T)[1:500]
colnames(subsets)[feature.select]
subsets.pca.scores[feature.select,1]
pdf(file = paste(args[1], "t.range01.subsets.noscale.pca.pdf", sep = "."), width = 15, height = 12)
ggplot(data = subsets.pca.scores, aes(x = PC1, y = PC2, label = rownames(subsets.pca.scores))) +
	geom_point(aes(alpha = 1/10)) +
	geom_text(colour = "tomato", alpha = 0.8, size = 2) +
	ggtitle(args[1])
dev.off()

### t-SNE
colors <- rainbow(length(unique(subsets.types)))
names(colors) <- unique(subsets.types)
## Executing the algorithm on curated data
set.seed(123)
tsne <- Rtsne(subsets, dims = 2, perplexity = 30, verbose = TRUE, max_iter = 500)
tsne.scale <- Rtsne(subsets, dims = 2, perplexity = 30, verbose = TRUE, max_iter = 500, pca_center = T, pca_scale = T)
t.tsne <- Rtsne(t(subsets), dims = 2, perplexity = 30, verbose = TRUE, max_iter = 500)
t.tsne.scale <- Rtsne(t(subsets), dims = 2, perplexity = 30, verbose = TRUE, max_iter = 500, pca_center = T, pca_scale = T)

# exeTimeTsne<- system.time(Rtsne(subsets, dims = 2, perplexity = 30, verbose = TRUE, max_iter = 500))

## Plotting
subsets.tsne <- as.data.frame(tsne$Y)
colnames(subsets.tsne) <- c("tSNE.1", "tSNE.2")
# png(filename = paste(args[1], "subsets.noscale.tsne.png", sep = "."), width = 1500, height = 1200)
pdf(file = paste(args[1], "subsets.noscale.tsne.pdf", sep = "."), width = 15, height = 12)
ggplot(data = subsets.tsne, aes(x = tSNE.1, y = tSNE.2)) + #, label = rownames(pca.scores))) +
	geom_point(aes(colour = subsets.types, alpha = 1/10)) +
	ggtitle(args[1])
dev.off()
subsets.tsne <- as.data.frame(tsne.scale$Y)
colnames(subsets.tsne) <- c("tSNE.1", "tSNE.2")
# png(filename = paste(args[1], "subsets.scale.tsne.png", sep = "."), width = 1500, height = 1200)
pdf(file = paste(args[1], "subsets.scale.tsne.pdf", sep = "."), width = 15, height = 12)
ggplot(data = subsets.tsne, aes(x = tSNE.1, y = tSNE.2)) + #, label = rownames(pca.scores))) +
	geom_point(aes(colour = subsets.types, alpha = 1/10)) +
	ggtitle(args[1])
dev.off()

subsets.tsne <- as.data.frame(t.tsne$Y)
# subsets.tsne <- data.frame(subsets.tsne, colnames(subsets))
colnames(subsets.tsne) <- c("tSNE.1", "tSNE.2") #, "labs")
rownames(subsets.tsne) <- colnames(subsets)
pdf(file = paste(args[1], "t.subsets.noscale.tsne.pdf", sep = "."), width = 15, height = 12)
ggplot(data = subsets.tsne, aes(x = tSNE.1, y = tSNE.2, label = rownames(subsets.tsne))) +
	geom_point(aes(alpha = 1/10)) +
	geom_text(colour = "tomato", alpha = 0.8, size = 2) +
	ggtitle(args[1])
dev.off()
subsets.tsne <- as.data.frame(t.tsne.scale$Y)
colnames(subsets.tsne) <- c("tSNE.1", "tSNE.2")
rownames(subsets.tsne) <- colnames(subsets)
pdf(file = paste(args[1], "t.subsets.scale.tsne.pdf", sep = "."), width = 15, height = 12)
ggplot(data = subsets.tsne, aes(x = tSNE.1, y = tSNE.2, label = rownames(subsets.tsne))) +
	geom_point(aes(alpha = 1/10)) +
	geom_text(colour = "tomato", alpha = 0.8, size = 2) +
	ggtitle(args[1])
dev.off()
# plot(tsne$Y, t='n', main="tsne")
# text(tsne$Y, labels = subsets.types, col = colors[subsets.types])

### randomForest
rawname <- colnames(subsets)
colnames(subsets) <- paste("var", 1:ncol(subsets),sep = "")
subsets.scale <- scale(subsets)
# subsets.scale <- data.frame(subsets.scale, subsets.types)
rf <- randomForest(subsets.scale) #subsets.scale[,169])
# rf <- randomForest(subsets.types ~ ., data = subsets.scale)
a <- importance(rf)
rownames(a) <- rawname
b <- a[order(a, decreasing = T),]#[1:50]
mydata <- data.frame(names(b),b)
colnames(mydata) <- c("group", "importance")
mydata$group <- factor(mydata$group, levels = mydata$group)
png(filename = paste(args[1], "subsets.unsupervised.rfimportance.png", sep = "."), width = 1500, height = 600)
ggplot(data = mydata, aes(x = group, y = importance)) +
	geom_bar(stat = "identity") + 
	ggtitle(args[1]) +
	theme(axis.text.x = element_text(angle = 60, hjust = 1))
dev.off()

c <- a[order(a, decreasing = T),]#[1:50]
colors <- colorRampPalette(c("blue", "white", "red"))(100)
# data <- scores[, order(a, decreasing = T)[1:50]]
colnames(subsets.scale) <- rawname
data <- subsets.scale[, order(a, decreasing = T)]
data[data>3] <- 3
data[data<-3] <- -3
png(filename = paste(args[1], "unsupervised.topimportance.heatmap.png", sep = "."), width = 1500, height = 1200)
myplot <- pheatmap(data, scale = "none", 
	show_rownames = F, show_colnames = T, color = colors, 
	cluster_cols = F, cluster_rows = F)
dev.off()


rf.1 <- randomForest(subsets.scale, subsets.types) #subsets.scale[,169])
# rf <- randomForest(subsets.types ~ ., data = subsets.scale)
a <- importance(rf.1)
rownames(a) <- rawname
b <- a[order(a, decreasing = T),]#[1:50]
mydata <- data.frame(names(b),b)
colnames(mydata) <- c("group", "importance")
mydata$group <- factor(mydata$group, levels = mydata$group)
png(filename = paste(args[1], "subsets.supervised.rfimportance.png", sep = "."), width = 1500, height = 600)
ggplot(data = mydata, aes(x = group, y = importance)) +
	geom_bar(stat = "identity") + 
	ggtitle(args[1]) +
	theme(axis.text.x = element_text(angle = 60, hjust = 1))
dev.off()

c <- a[order(a, decreasing = T),]#[1:50]
colors <- colorRampPalette(c("blue", "white", "red"))(100)
# data <- scores[, order(a, decreasing = T)[1:50]]
colnames(subsets.scale) <- rawname
data <- subsets.scale[, order(a, decreasing = T)]
data[data>3] <- 3
data[data<-3] <- -3
png(filename = paste(args[1], "supervised.topimportance.heatmap.png", sep = "."), width = 1500, height = 1200)
myplot <- pheatmap(data, scale = "none", 
	show_rownames = F, show_colnames = T, color = colors, 
	cluster_cols = F, cluster_rows = F)
dev.off()

save.image("hkgvstsg.RData")

###################################################

myPCA <- prcomp(scores, scale. = T, center = T)
pca.scores <- as.data.frame(myPCA$x)
png(filename = paste(args[1], "pca.png", sep = "."), width = 1500, height = 1200)
ggplot(data = pca.scores, aes(x = PC1, y = PC2)) + #, label = rownames(pca.scores))) +
	geom_point(aes(colour = types, alpha = 1/10)) +
	ggtitle(args[1])
dev.off()
myPCA.t <- prcomp(t(scores))
pca.scores.t <- as.data.frame(myPCA.t$x)
png(filename = paste(args[1], "pca.t.png", sep = "."), width = 1500, height = 1500)
ggplot(data = pca.scores.t, aes(x = PC1, y = PC2, label = rownames(pca.scores.t))) +
	geom_point(aes(alpha = 1/10)) + geom_text(colour = "tomato", alpha = 0.8, size = 4) +
	ggtitle(args[1])
dev.off()
# feature.select <- order(abs(pca.scores.t[,1]),decreasing=T)[1:50]
# colnames(scores)[feature.select]

set.seed(123)
a.tsne <- Rtsne(scores, check_duplicates = F, dims = 2, perplexity = 30, verbose = TRUE, max_iter = 500)
a.tsne.scale <- Rtsne(scores, check_duplicates = F, dims = 2, perplexity = 30, verbose = TRUE, max_iter = 500, pca_center = T, pca_scale = T)

## Plotting
tsne.scores <- as.data.frame(a.tsne$Y)
colnames(tsne.scores) <- c("tSNE.1", "tSNE.2")
# png(filename = paste(args[1], "subsets.noscale.tsne.png", sep = "."), width = 1500, height = 1200)
pdf(file = paste(args[1], "noscale.tsne.pdf", sep = "."), width = 15, height = 12)
ggplot(data = tsne.scores, aes(x = tSNE.1, y = tSNE.2)) + #, label = rownames(pca.scores))) +
	geom_point(aes(colour = types, alpha = 1/10)) +
	ggtitle(args[1])
dev.off()
tsne.scores <- as.data.frame(a.tsne.scale$Y)
colnames(tsne.scores) <- c("tSNE.1", "tSNE.2")
# png(filename = paste(args[1], "subsets.scale.tsne.png", sep = "."), width = 1500, height = 1200)
pdf(file = paste(args[1], "scale.tsne.pdf", sep = "."), width = 15, height = 12)
ggplot(data = tsne.scores, aes(x = tSNE.1, y = tSNE.2)) + #, label = rownames(pca.scores))) +
	geom_point(aes(colour = types, alpha = 1/10)) +
	ggtitle(args[1])
dev.off()




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
