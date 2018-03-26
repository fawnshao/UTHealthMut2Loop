## conda.R
# .libPaths("/work/04935/shaojf/stampede2/myTools/anaconda2/lib/R/library")
library(data.table)
args <- c("GTEx.log2tpm.median.tsv", "GTEx_sample.tissue.txt", "GTEx.gene.class.byGTExorder")
raw.table <- fread(args[1], sep = "\t", header = T)
tpm <- data.matrix(raw.table[,-1])
rownames(tpm) <- as.matrix(raw.table[,1])
tissues.lab <- colnames(raw.table)[-1]
colnames(tpm) <- tissues.lab
class <- fread(args[3], sep = "\t", header = T)
class.lab <- as.matrix(class[,2])[,1]

#####################
library(ggplot2)
# class.lab <- as.list(class.lab)
# tissues.lab <- as.list(tissues.lab)

myPCA <- prcomp(tpm, scale. = F, center = F)
pca.scores <- as.data.frame(myPCA$x)
p1 <- ggplot(data = pca.scores, aes(x = PC1, y = PC2)) +
	geom_point(aes(colour = class.lab, alpha = 1/10)) +
	ggtitle(args[1])
png(filename = paste(args[1], "pca.png", sep = "."), width = 1500, height = 1200)
print(p1)
dev.off()
pdf(file = paste(args[1], "pca.pdf", sep = "."), width = 15, height = 12)
print(p1)
dev.off()
class.lab.sim <- class.lab
class.lab.sim[class.lab.sim!="hkg" & class.lab.sim!="mixTSG" & class.lab.sim!="Testis" & class.lab.sim!="other"] <- "othersingleTSG"
p11 <- ggplot(data = pca.scores, aes(x = PC1, y = PC2)) +
	geom_point(aes(colour = class.lab.sim, alpha = 1/10)) +
	ggtitle(args[1])
png(filename = paste(args[1], "pca.sim.png", sep = "."), width = 1500, height = 1200)
print(p11)
dev.off()
pdf(file = paste(args[1], "pca.sim.pdf", sep = "."), width = 15, height = 12)
print(p11)
dev.off()


subs <- pca.scores[class.lab!="other",]
subs.class <- class.lab[class.lab!="other"]
p2 <- ggplot(data = subs, aes(x = PC1, y = PC2)) +
	geom_point(aes(colour = subs.class, alpha = 1/10)) +
	ggtitle(args[1])
png(filename = paste(args[1], "pca.sub.png", sep = "."), width = 1500, height = 1200)
print(p2)
dev.off()
pdf(file = paste(args[1], "pca.sub.pdf", sep = "."), width = 15, height = 12)
print(p2)
dev.off()
subs.class.sim <- subs.class
subs.class.sim[subs.class.sim!="hkg" & subs.class.sim!="mixTSG" & subs.class.sim!="Testis"] <- "othersingleTSG"
p3 <- ggplot(data = subs, aes(x = PC1, y = PC2)) +
	geom_point(aes(colour = subs.class.sim, alpha = 1/10)) +
	ggtitle(args[1])
png(filename = paste(args[1], "pca.sub.sim.png", sep = "."), width = 1500, height = 1200)
print(p3)
dev.off()
pdf(file = paste(args[1], "pca.sub.sim.pdf", sep = "."), width = 15, height = 12)
print(p3)
dev.off()
save.image("median.PCA.RData")


#####################
# conda.R
load("median.PCA.RData")
.libPaths("/work/04935/shaojf/stampede2/myTools/anaconda2/lib/R/library")
library(ggplot2)
library(Rtsne)
set.seed(123)
tsne <- Rtsne(tpm, check_duplicates = FALSE, dims = 2, perplexity = 30, verbose = TRUE, max_iter = 500)
# tsne.scale <- Rtsne(tpm, check_duplicates = FALSE, dims = 2, perplexity = 30, verbose = TRUE, max_iter = 500, pca_center = T, pca_scale = T)
t.tsne <- Rtsne(t(tpm), check_duplicates = FALSE, dims = 2, perplexity = 10, verbose = TRUE, max_iter = 500)
# t.tsne.scale <- Rtsne(t(tpm), check_duplicates = FALSE, dims = 2, perplexity = 30, verbose = TRUE, max_iter = 500, pca_center = T, pca_scale = T)
tpm.sub <- tpm[class.lab!="other",]
tsne.sub <- Rtsne(tpm.sub, check_duplicates = FALSE, dims = 2, perplexity = 30, verbose = TRUE, max_iter = 500)
tsne.scale.sub <- Rtsne(tpm.sub, check_duplicates = FALSE, dims = 2, perplexity = 30, verbose = TRUE, max_iter = 500, pca_center = T, pca_scale = T)
t.tsne.sub <- Rtsne(t(tpm.sub), check_duplicates = FALSE, dims = 2, perplexity = 10, verbose = TRUE, max_iter = 500)
t.tsne.scale.sub <- Rtsne(t(tpm.sub), check_duplicates = FALSE, dims = 2, perplexity = 10, verbose = TRUE, max_iter = 500, pca_center = T, pca_scale = T)

# save.image("median.tsne.RData")
save.image("median.tsne.1.RData")
# load("median.tsne.RData")

## Plotting
v.tsne <- as.data.frame(tsne.sub$Y)
colnames(v.tsne) <- c("tSNE.1", "tSNE.2")
# png(filename = paste(args[1], "subsets.noscale.tsne.png", sep = "."), width = 1500, height = 1200)
pdf(file = paste(args[1], "tsne.sub.pdf", sep = "."), width = 15, height = 12)
ggplot(data = v.tsne, aes(x = tSNE.1, y = tSNE.2)) + 
	geom_point(aes(colour = subs.class.sim, alpha = 1/10)) +
	ggtitle(args[1])
dev.off()
v.tsne <- as.data.frame(tsne.scale.sub$Y)
colnames(v.tsne) <- c("tSNE.1", "tSNE.2")
# png(filename = paste(args[1], "subsets.noscale.tsne.png", sep = "."), width = 1500, height = 1200)
pdf(file = paste(args[1], "tsne.scale.sub.pdf", sep = "."), width = 15, height = 12)
ggplot(data = v.tsne, aes(x = tSNE.1, y = tSNE.2)) + 
	geom_point(aes(colour = subs.class.sim, alpha = 1/10)) +
	ggtitle(args[1])
dev.off()


v.tsne <- as.data.frame(t.tsne.sub$Y)
colnames(v.tsne) <- c("tSNE.1", "tSNE.2")
# png(filename = paste(args[1], "subsets.noscale.tsne.png", sep = "."), width = 1500, height = 1200)
pdf(file = paste(args[1], "t.tsne.sub.pdf", sep = "."), width = 15, height = 12)
ggplot(data = v.tsne, aes(x = tSNE.1, y = tSNE.2)) + 
	geom_point(aes(colour = as.vector(unlist(tissues.lab)), alpha = 1/10)) +
	ggtitle(args[1])
dev.off()
v.tsne <- as.data.frame(t.tsne.scale.sub$Y)
colnames(v.tsne) <- c("tSNE.1", "tSNE.2")
# png(filename = paste(args[1], "subsets.noscale.tsne.png", sep = "."), width = 1500, height = 1200)
pdf(file = paste(args[1], "t.tsne.scale.sub.pdf", sep = "."), width = 15, height = 12)
ggplot(data = v.tsne, aes(x = tSNE.1, y = tSNE.2)) + 
	geom_point(aes(colour = as.vector(unlist(tissues.lab)), alpha = 1/10)) +
	ggtitle(args[1])
dev.off()

v.tsne <- as.data.frame(tsne$Y)
colnames(v.tsne) <- c("tSNE.1", "tSNE.2")
# png(filename = paste(args[1], "subsets.noscale.tsne.png", sep = "."), width = 1500, height = 1200)
pdf(file = paste(args[1], "tsne.pdf", sep = "."), width = 15, height = 12)
ggplot(data = v.tsne, aes(x = tSNE.1, y = tSNE.2)) + 
	geom_point(aes(colour = class.lab.sim, alpha = 1/10)) +
	ggtitle(args[1])
dev.off()
v.tsne <- as.data.frame(t.tsne$Y)
colnames(v.tsne) <- c("tSNE.1", "tSNE.2")
# png(filename = paste(args[1], "subsets.noscale.tsne.png", sep = "."), width = 1500, height = 1200)
pdf(file = paste(args[1], "t.tsne.pdf", sep = "."), width = 15, height = 12)
ggplot(data = v.tsne, aes(x = tSNE.1, y = tSNE.2)) + 
	geom_point(aes(colour = as.vector(unlist(tissues.lab)))) + #, alpha = 7/10)) +
	scale_colour_manual(values = rainbow(length(as.vector(unlist(tissues.lab))))) + 
	ggtitle(args[1])
dev.off()
