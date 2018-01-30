args <- commandArgs(TRUE)
library(pheatmap)
# args <- c("div.gene.pairs.annotation.txt", "div.gene.tpm.gct", "GTEx_sample.tissue.txt")
pairs <- as.matrix(read.table(args[1], sep = "\t"))
exprs <- as.matrix(read.table(args[2], sep = "\t", row.names = 1))
info <- as.matrix(read.table(args[3], sep = "\t"))
tissues <- unique(info[,2])
genes <- rownames(exprs)
res <- matrix(ncol = 4, nrow = nrow(pairs))
ave <- matrix(ncol = 2, nrow = nrow(pairs))
fc <- matrix(ncol = ncol(exprs), nrow = nrow(pairs))
ave.rev.tissues <- matrix(ncol = length(tissues), nrow = nrow(pairs))
ave.fwd.tissues <- matrix(ncol = length(tissues), nrow = nrow(pairs))

for(i in 1:nrow(pairs)){
	x1 <- exprs[genes==pairs[i,1],]
	x2 <- exprs[genes==pairs[i,2],]
	if(length(x1) == length(x2) && length(x1) == ncol(exprs)){
		cor.x <- cor.test(x1, x2, method = "pearson")
		cor.y <- cor.test(x1, x2, method = "spearman")
		res[i,] <- c(cor.x$estimate[[1]], cor.x$p.value, 
			cor.y$estimate[[1]], cor.y$p.value)
		fc[i,] <- (x2 + 1) / (x1 + 1)
	}
	ave[i,] <- c(mean(x1), mean(x2))

	for(j in 1:length(tissues)){
		ave.rev.tissues[i,j] <- mean(x1[info[,2]==tissues[j]])
		ave.fwd.tissues[i,j] <- mean(x2[info[,2]==tissues[j]])
	}
}
out <- data.frame(pairs, ave, res)
colnames(out) <- c("gene.rev","gene.fwd",
	"gene.rev.type","gene.fwd.type",
	"ave.gene.rev","ave.gene.fwd",
	"pearson.cor", "pearson.pvalues", 
	"spearman.cor", "spearman.pvalues")
write.table(out, file = paste(args[1], "cor.tsv", sep = "."), sep = "\t", row.names = FALSE, quote = FALSE)

colnames(fc) <- colnames(exprs)
gene.pair <- paste(pairs[,2], pairs[,1], sep = "/")
write.table(data.frame(gene.pair, fc), file = paste(args[1], "fc.tsv", sep = "."), sep = "\t", row.names = FALSE, quote = FALSE)

colnames(ave.rev.tissues) <- tissues
colnames(ave.fwd.tissues) <- tissues
write.table(data.frame(pairs[,1], ave.rev.tissues), file = paste(args[1], "rev.mean.tsv", sep = "."), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(pairs[,2], ave.fwd.tissues), file = paste(args[1], "fwd.mean.tsv", sep = "."), sep = "\t", row.names = FALSE, quote = FALSE)
fc.mean <- ave.rev.tissues / ave.fwd.tissues
write.table(data.frame(gene.pair, fc.mean), file = paste(args[1], "fc.mean.tsv", sep = "."), sep = "\t", row.names = FALSE, quote = FALSE)


x <- data.frame(pairs[,3:4], ave.rev.tissues, ave.fwd.tissues)
rownames(x) <- paste(pairs[,1], pairs[,2], sep = "<->")
nacount <- apply(x[,-c(1:2)], 1, function(x){sum(is.na(x))})
# infcount <- apply(x, 1, function(x){sum(is.infinite(x))})
data <- x[nacount < length(tissues),-c(1:2)]
anno <- as.matrix(x[nacount < length(tissues), c(1:2)])
anno[anno != "protein_coding"] <- "ncRNA"
colnames(anno) <- c("reverse", "forward")
# row_annos <- data.frame(Group = factor(paste(anno[,1], anno[,2], sep = "<->")))
row_annos <- data.frame(forward = factor(anno[,2]), reverse = factor(anno[,1]))
rownames(row_annos) <- rownames(data)
ann_colors = list(
    forward = c(ncRNA = "#1B9E77", protein_coding = "#D95F02"),
    reverse = c(ncRNA = "#1B9E77", protein_coding = "#D95F02")
)
# infcount <- apply(data, 1, function(x){sum(is.infinite(x))})
# data <- data[infcount == 0,]
output <- "div.gene.pairs.fc.mean"
# colors <- colorRampPalette(c("blue", "white", "red"))(100)

# breaklists <- c(seq(0,8,by=0.1),seq(8.5,10,by=0.5),seq(11,15,by=1))
breaklists <- c(seq(0,8,by=0.1),seq(9,15,by=1))
colorn <- length(breaklists)
colors <- colorRampPalette(c("blue", "white", "red"))(colorn)

dist_methods <- c("correlation", "euclidean", "maximum", "manhattan", "canberra", "minkowski")
clust_methods <- c("ward.D", "ward.D2", "complete", "average", "mcquitty", "median", "centroid")
i <- 2
j <- 2
png(filename = paste(output, dist_methods[i], clust_methods[j], "png", sep = "."), width = 2000, height = 1000)
p1 <- pheatmap(log(data + 1, 2), scale = "none", show_rownames = F, show_colnames = F, 
         color = colors, cluster_cols = F, clustering_distance_rows = dist_methods[i], 
         clustering_method = clust_methods[j], breaks = breaklists,
         annotation_row = row_annos, 
         cutree_row = 8, annotation_colors = ann_colors)
dev.off()

pairscluster <- cutree(p1$tree_row, k = 8)
write.table(data.frame(pairscluster, anno, rownames(data), data), file = paste(args[1], "fc.mean.pheatmap.tsv", sep = "."), sep = "\t", row.names = FALSE, quote = FALSE)

save.image("div.gene.RData")