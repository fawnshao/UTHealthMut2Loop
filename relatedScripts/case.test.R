.libPaths("/work/04935/shaojf/stampede2/myTools/anaconda2/lib/R/library")
library(ggplot2)
library(pheatmap)
library(data.table)
library(Rtsne)
set.seed(123)
args <- "hkg.gm.liver.HiChIP-Loop.count.sim.mat"
input <- fread(args[1], sep = "\t", header = T, na.strings = "/")

Plotheatmap <- function(x, y = annos, z = colors, pre = args[1]){
	png(filename = paste(pre, "pheatmap.png", sep = "."), width = 1000, height = 800)
	# pdf(file = paste(pre, "pheatmap.pdf", sep = "."), width = 12, height = 12)
	myplot <- pheatmap(x, scale = "none", annotation_row = y,
		show_rownames = F, show_colnames = T, color = z,
		cluster_cols = F, cluster_rows = F)
	dev.off()
}
Plotheatmap.cluster <- function(x, y = annos, z = colors, pre = args[1]){
	png(filename = paste(pre, "pheatmap.cluster.png", sep = "."), width = 1000, height = 800)
	# pdf(file = paste(pre, "pheatmap.pdf", sep = "."), width = 12, height = 12)
	myplot <- pheatmap(x, scale = "none", annotation_row = y,
		show_rownames = F, show_colnames = T, color = z,
		cluster_cols = T, cluster_rows = T)
	dev.off()
}
Plotheatmap.cluster.withname <- function(x, y = annos, z = colors, pre = args[1]){
	png(filename = paste(pre, "pheatmap.cluster.png", sep = "."), width = 1000, height = 1800)
	# pdf(file = paste(pre, "pheatmap.pdf", sep = "."), width = 10, height = 16)
	myplot <- pheatmap(x, scale = "none", annotation_row = y,
		show_rownames = T, show_colnames = T, color = z, fontsize_row = 1,
		cluster_cols = T, cluster_rows = T)
	dev.off()
}

scores <- data.matrix(input[,-c(1:2)])
genes <- as.matrix(input[,1])
class <- as.matrix(input[,2])
class.sim <- class
class.sim[1:2647] <- "HKG"
types <- factor(class)
types.sim <- factor(class.sim)
rownames(scores) <- paste(1:nrow(scores), genes, sep = ": ")

annos <- data.frame(class = types)
rownames(annos) <- rownames(scores)
colors <- colorRampPalette(c("white", "blue"))(10)

data.x <- scores
data.x[is.na(data.x)] <- 0
data.x[data.x > 1] <- 1
Plotheatmap(data.x, pre = args[1])
Plotheatmap.cluster(data.x, pre = args[1])

Plotheatmap.cluster(data.x[1:2647,], pre = paste(args[1], "HKG", sep = "."))
Plotheatmap.cluster(data.x[2648:3051,], pre = paste(args[1], "nonHKG", sep = "."))

Plotheatmap.cluster.withname(data.x[1:2647,], pre = paste(args[1], "HKG", sep = "."))
Plotheatmap.cluster.withname(data.x[2648:3051,], pre = paste(args[1], "nonHKG", sep = "."))



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
loop.tsne <- RuntSNE(x = data.x, pre = paste(args[1],"loop", sep = "."))
loop.tsne.plot <- PlottSNE(x = loop.tsne, pre = paste(args[1],"loop", sep = "."))

loop.tsne <- RuntSNE(x = data.x[2648:3051,], pre = paste(args[1],"nonHKGloop", sep = "."))
loop.tsne.plot <- PlottSNE(x = loop.tsne, y = types[2648:3051], z = types.sim[2648:3051], pre = paste(args[1],"nonHKGloop", sep = "."))
