library(ggplot2)
library(pheatmap)
library(data.table)

# args <- c("targets.tsv")
# args <- c("targets.div.tsv")
args <- c("close.targets.div.tsv")
input <- fread(args[1], sep = "\t", header = T, na.strings = "/")
pairs <- paste(as.matrix(input[,1])[,1], as.matrix(input[,2])[,1], sep = ": ")

anno.class <- as.matrix(input[,5])[,1]
anno.class[anno.class!="protein_coding"] <- "others"
anno.class.1 <- anno.class
anno.class <- as.matrix(input[,6])[,1]
anno.class[anno.class!="protein_coding"] <- "others"
anno.class.2 <- anno.class
annos <- data.frame(genes = anno.class.1, gene.neighbors = anno.class.2)
rownames(annos) <- pairs

datax <- as.matrix(input[,-c(1:6)])
datax[datax > 10] <- 10
rownames(datax) <- pairs

colors <- colorRampPalette(c("blue", "white", "red"))(100)

png(filename = paste(args[1], "pairs.pheatmap.png", sep = "."), width = 3000, height = 2000)
pheatmap(as.matrix(datax), scale = "none", annotation_row = annos, fontsize_row = 5,
	show_rownames = T, show_colnames = T, color = colors, # gaps_row = c(2552, 5287),
	cluster_cols = F, cluster_rows = F)
dev.off()

pdf(file = paste(args[1], "pairs.pheatmap.pdf", sep = "."), width = 30, height = 20)
pheatmap(as.matrix(datax), scale = "none", annotation_row = annos, fontsize_row = 5,
	show_rownames = T, show_colnames = T, color = colors, # gaps_row = c(2552, 5287),
	cluster_cols = F, cluster_rows = F)
dev.off()
