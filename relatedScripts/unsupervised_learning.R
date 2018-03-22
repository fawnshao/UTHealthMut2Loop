## conda.R
# .libPaths("/work/04935/shaojf/stampede2/myTools/anaconda2/lib/R/library")
library(data.table)
args <- c("onlyhkg.tsg.log2tpm.median.tsv")
raw.table <- fread(args[1], sep = "\t", header = T)
tpm <- data.matrix(raw.table[,-c(1,2)])
rownames(tpm) <- as.matrix(raw.table[,1])
tissues.lab <- colnames(raw.table)[-c(1,2)]
colnames(tpm) <- tissues.lab
class <- as.matrix(raw.table[,2])[,1]

#####################
library(ggplot2)
tpm.hkg <- tpm[class == "hkg",]
myPCA <- prcomp(tpm.hkg, scale. = F, center = F)
myPCA.t <- prcomp(t(tpm.hkg), scale. = F, center = F)
myPCA.scale <- prcomp(tpm.hkg, scale. = T, center = T)

pca.scores <- as.data.frame(myPCA$x)
gene.class <- cut(pca.scores[,1], breaks = quantile(pca.scores[,1], probs = seq(0, 1, 0.1)), 
	include.lowest = TRUE, labels = 1:10)
write.csv(data.frame(rownames(pca.scores),gene.class), paste(args[1], "onlyhkg.pca.csv", sep = "."))

p1 <- ggplot(data = pca.scores, aes(x = PC1, y = PC2)) +
	geom_point(aes(colour = "red")) +
	# geom_text(data = subset(pca.scores, wt > 4 | mpg > 25), aes(wt, mpg, label = name))
	ggtitle(args[1])
png(filename = paste(args[1], "onlyhkg.pca.png", sep = "."), width = 1500, height = 1200)
print(p1)
dev.off()
pdf(file = paste(args[1], "onlyhkg.pca.pdf", sep = "."), width = 15, height = 12)
print(p1)
dev.off()

pca.scores <- as.data.frame(myPCA.t$x)
p2 <- ggplot(data = pca.scores, aes(x = PC1, y = PC2, label = rownames(pca.scores))) +
	geom_point(aes(colour = tissues.lab)) +
	geom_text(aes(colour = tissues.lab), alpha = 0.8, size = 5) +
	theme(legend.position = "none") +
	ggtitle(args[1])
png(filename = paste(args[1], "onlyhkg.pca.t.png", sep = "."), width = 1500, height = 1200)
print(p2)
dev.off()
pdf(file = paste(args[1], "onlyhkg.pca.t.pdf", sep = "."), width = 15, height = 12)
print(p2)
dev.off()

pca.scores <- as.data.frame(myPCA.scale$x)
p3 <- ggplot(data = pca.scores, aes(x = PC1, y = PC2)) +
	geom_point(aes(colour = "red")) +
	ggtitle(args[1])
png(filename = paste(args[1], "onlyhkg.pca.scale.png", sep = "."), width = 1500, height = 1200)
print(p3)
dev.off()
pdf(file = paste(args[1], "onlyhkg.pca.scale.pdf", sep = "."), width = 15, height = 12)
print(p3)
dev.off()
save.image("onlyhkg.pca.RData")



######
load("onlyhkg.pca.RData")
library(Rtsne)
.libPaths("/work/04935/shaojf/stampede2/myTools/anaconda2/lib/R/library")
library(ggplot2)
set.seed(123)
tsne <- Rtsne(tpm.hkg, check_duplicates = FALSE, dims = 2, perplexity = 30, verbose = TRUE, max_iter = 500)
v.tsne <- as.data.frame(tsne$Y)
colnames(v.tsne) <- c("tSNE.1", "tSNE.2")
pdf(file = paste(args[1], "onlyhkg.tsne.pdf", sep = "."), width = 15, height = 12)
ggplot(data = v.tsne, aes(x = tSNE.1, y = tSNE.2)) + 
	geom_point(aes(colour = gene.class, alpha = 1/10)) +
	ggtitle(args[1])
dev.off()

######
final.class <- c(paste(class[1:length(gene.class)], gene.class, sep = ""), class[(length(gene.class) + 1):length(class)])
tsne.all <- Rtsne(tpm, check_duplicates = FALSE, dims = 2, perplexity = 30, verbose = TRUE, max_iter = 500)
v.tsne <- as.data.frame(tsne.all$Y)
colnames(v.tsne) <- c("tSNE.1", "tSNE.2")
pdf(file = paste(args[1], "tsne.10hkg.pdf", sep = "."), width = 15, height = 12)
ggplot(data = v.tsne, aes(x = tSNE.1, y = tSNE.2)) + 
	geom_point(aes(colour = final.class, alpha = 1/10)) +
	ggtitle(args[1])
dev.off()

final.class.sim <- final.class
final.class.sim[class!="hkg" & final.class.sim!="mixTSG" & final.class.sim!="Testis"] <- "othersingleTSG"
pdf(file = paste(args[1], "tsne.10hkg.sim.pdf", sep = "."), width = 15, height = 12)
ggplot(data = v.tsne, aes(x = tSNE.1, y = tSNE.2)) + 
	geom_point(aes(colour = final.class.sim, alpha = 1/10)) +
	ggtitle(args[1])
dev.off()

save.image("onlyhkg_tsg.pca.tsne.RData")


