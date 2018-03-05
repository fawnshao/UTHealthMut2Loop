library(data.table)
library(ggplot2)
library(pheatmap)
args <- c("housekeepinggene.GTRD.mat", 
           "tissuespecificgene.GTRD.mat",
           "GTEx.GTRD.mat","GTRD.maxcells.GTEx.hkg.tsg.txt")
# args <- c("housekeepinggene.GTRD.count.mat", 
#           "tissuespecificgene.GTRD.count.mat",
#           "GTEx.GTRD.count.mat","GTRD.count.GTEx.hkg.tsg.txt")
# hkg <- read.table(args[1], sep = "\t", header = T, row.names = 1)
hkg <- fread(args[1], sep = "\t", header = T)
tsg <- fread(args[2], sep = "\t", header = T)
all <- fread(args[3], sep = "\t", header = T)
anno <- as.matrix(read.table(args[4], sep = "\t", header = T, row.names = 1))

hkg.score <- hkg[,-1]
tsg.score <- tsg[,-1]
all.score <- all[,-1]
# row_annos <- data.frame(type = factor(as.matrix(anno[,2])))
rownames(hkg.score) <- as.matrix(hkg[,1])
rownames(tsg.score) <- as.matrix(tsg[,1])
rownames(all.score) <- as.matrix(all[,1])
# rownames(row_annos) <- as.matrix(anno[,1])
# dim(); some HKG or TSG is overlapped with none of the peaks. 
# so the final line count is different from orignals
colorn <- 10
colors <- colorRampPalette(c("blue", "white", "red"))(colorn)

# data <- log2(data.matrix(hkg.score) + 1)
data <- data.matrix(hkg.score)
data[data > 0 & data < 2] <- 1
data[data > 1] <- 2
png(filename = "hkg.GTRD.png", width = 1500, height = 1200)
myplot <- pheatmap(data, scale = "none", 
	show_rownames = F, show_colnames = F, color = colors, 
	cluster_cols = F, cluster_rows = T)
dev.off()
cluster <- cutree(myplot$tree_row, k = 5)
write.table(data.frame(cluster[myplot$tree_row$order], 
	rownames(data)[myplot$tree_row$order],
	data[myplot$tree_row$order, ]), 
	file = "hkg.GTRD.tsv", 
	sep = "\t", row.names = FALSE, quote = FALSE)
# a <- data

# data <- log2(data.matrix(tsg.score) + 1)
data <- data.matrix(tsg.score)
data[data > 0 & data < 2] <- 1
data[data > 1] <- 2
png(filename = "tsg.GTRD.png", width = 1500, height = 1200)
myplot <- pheatmap(data, scale = "none", 
	show_rownames = F, show_colnames = F, color = colors, 
	cluster_cols = F, cluster_rows = T)
dev.off()
cluster <- cutree(myplot$tree_row, k = 5)
write.table(data.frame(cluster[myplot$tree_row$order], 
	rownames(data)[myplot$tree_row$order],
	data[myplot$tree_row$order, ]), 
	file = "tsg.GTRD.tsv", 
	sep = "\t", row.names = FALSE, quote = FALSE)
# b <- data

######### different columns?
# data <- rbind(a,b)

# data <- rbind(data.matrix(hkg.score), data.matrix(tsg.score))
# data[data > 0 & data < 2] <- 1
# data[data > 1] <- 2
# row_annos <- data.frame(type = factor(c(rep("hkg",nrow(hkg.score)),rep("tsg",nrow(tsg.score)))))
# rownames(row_annos) <- rownames(data)
# subrow <- sample(1:nrow(all.score), 3000)
# data <- data.matrix(all.score)[subrow,]
# row_annos <- data.frame(type = factor(anno[subrow,1]))
data <- data.matrix(all.score)
row_annos <- data.frame(type = factor(anno[,1]))
rownames(row_annos) <- rownames(data)
ann_colors <- list(type = c(hkg = "#1B9E77", tsg = "#D95F02", other = "azure4"))
data[data > 0 & data < 2] <- 1
data[data > 1] <- 2
png(filename = "hkg_tsg.GTRD.1.png", width = 1500, height = 1200)
myplot <- pheatmap(data, scale = "none", 
	show_rownames = F, show_colnames = F, color = colors, 
	cluster_cols = F, cluster_rows = T,
	annotation_row = row_annos, annotation_colors = ann_colors,
	cutree_row = 5)
dev.off()
cluster <- cutree(myplot$tree_row, k = 5)
write.table(data.frame(cluster[myplot$tree_row$order], 
	rownames(data)[myplot$tree_row$order],
	data[myplot$tree_row$order, ]), 
	file = "hkg_tsg.GTRD.1.tsv", 
	sep = "\t", row.names = FALSE, quote = FALSE)
# Error in col[as.numeric(cut(x, breaks = breaks, include.lowest = T))] : 
#   object of type 'closure' is not subsettable



### PCA
data <- data.matrix(all.score)
types <- factor(anno[,1])
# myPCA <- prcomp(data, scale. = F, center = F)
myPCA <- prcomp(data, scale. = T, center = T)
pca.scores <- as.data.frame(myPCA$x)
# png(filename = "hkg_tsg.GTRD.PCA.png", width = 1500, height = 1200)
# # plot(myPCA$x[,1:2], col = factor(anno[,1]))
# # biplot(myPCA)
# ggplot(data = pca.scores[anno[,1]!="other",], aes(x = PC1, y = PC2)) + #, label = rownames(pca.scores))) +
# 	geom_point(aes(colour = factor(anno[anno[,1]!="other",1]), alpha = 1/10)) +
# 	# geom_hline(yintercept = 0, colour = "gray65") +
# 	# geom_vline(xintercept = 0, colour = "gray65") +
# 	# geom_text(colour = "tomato", alpha = 0.8, size = 4) +
# 	ggtitle("PCA plot of GTRD")
# dev.off()
png(filename = "hkg_tsg.GTRD.maxcells.PCA.png", width = 1500, height = 1200)
ggplot(data = pca.scores[anno[,1]!="other",], aes(x = PC1, y = PC2)) + 
	geom_point(aes(colour = factor(anno[anno[,1]!="other",1]), alpha = 1/10)) +
	ggtitle("PCA plot of GTRD max cell count")
dev.off()
myPCA <- prcomp(t(data), scale. = T, center = T)
pca.scores <- as.data.frame(myPCA$x)
png(filename = "hkg_tsg.GTRD.maxcells.PCA.t.png", width = 1500, height = 1200)
ggplot(data = pca.scores, aes(x = PC1, y = PC2, label = rownames(pca.scores))) + 
	geom_point(aes(alpha = 1/10)) +
	geom_text(colour = "tomato", alpha = 0.8, size = 4) +
	ggtitle("PCA plot of GTRD max cell count (t)")
dev.off()

data <- data.matrix(hkg.score)
myPCA <- prcomp(t(data), scale. = T, center = T)
pca.scores <- as.data.frame(myPCA$x)
png(filename = "hkg.GTRD.maxcells.PCA.t.png", width = 1500, height = 1200)
ggplot(data = pca.scores, aes(x = PC1, y = PC2, label = rownames(pca.scores))) + 
	geom_point(aes(alpha = 1/10)) +
	geom_text(colour = "tomato", alpha = 0.8, size = 4) +
	ggtitle("PCA plot of HKG GTRD max cell count (t)")
dev.off()
data <- data.matrix(tsg.score)
myPCA <- prcomp(t(data), scale. = T, center = T)
pca.scores <- as.data.frame(myPCA$x)
png(filename = "tsg.GTRD.maxcells.PCA.t.png", width = 1500, height = 1200)
ggplot(data = pca.scores, aes(x = PC1, y = PC2, label = rownames(pca.scores))) + 
	geom_point(aes(alpha = 1/10)) +
	geom_text(colour = "tomato", alpha = 0.8, size = 4) +
	ggtitle("PCA plot of TSG GTRD max cell count (t)")
dev.off()




# save.image("GTRD.RData")

#####################
args <- c("housekeepinggene.GTRD.count.txt", 
          "tissuespecificgene.GTRD.count.txt") #, 
          # "GTEx.100way.tab")
hkg <- read.table(args[1], sep = "\t", header = F)
tsg <- read.table(args[2], sep = "\t", header = F)
# all <- fread(args[3], sep = "\t", header = F)

makematrix <- function(x) {
	col.lists <- unique(x[,1])
	row.lists <- unique(x[,2])
	out.matrix <- matrix(nrow = length(row.lists), ncol = length(col.lists))
	for(i in 1:nrow(out.matrix)){
		for(j in 1:ncol(out.matrix)){
			a <- x[x[,1]==col.lists[j] & x[,2]==row.lists[i],3]
			if(length(a) == 1){
				out.matrix[i,j] <- a
			}
		}
	}
	colnames(out.matrix) <- col.lists
	rownames(out.matrix) <- row.lists
	return(out.matrix)
}
hkg.score <- makematrix(hkg)
hkg.score <- makematrix(hkg)
write.csv(hkg.score, "hkg.GTRD.occurence.csv")
write.csv(tsg.score, "tsg.GTRD.occurence.csv")
# all.dis <- unique(all[,c(1,6)])

colorn <- 10
colors <- colorRampPalette(c("white", "blue"))(colorn)

data <- hkg.score
png(filename = "hkg.GTRD.peaks.png", width = 1500, height = 1200)
myplot <- pheatmap(data, scale = "none", 
	show_rownames = F, show_colnames = F, color = colors, 
	cluster_cols = F, cluster_rows = T)
dev.off()
cluster <- cutree(myplot$tree_row, k = 5)
write.table(data.frame(cluster[myplot$tree_row$order], 
	rownames(data)[myplot$tree_row$order],
	data[myplot$tree_row$order, ]), 
	file = "hkg.GTRD.peaks.tsv", 
	sep = "\t", row.names = FALSE, quote = FALSE)

data <- tsg.score
png(filename = "tsg.GTRD.peaks.png", width = 1500, height = 1200)
myplot <- pheatmap(data, scale = "none", 
	show_rownames = F, show_colnames = F, color = colors, 
	cluster_cols = F, cluster_rows = T)
dev.off()
cluster <- cutree(myplot$tree_row, k = 5)
write.table(data.frame(cluster[myplot$tree_row$order], 
	rownames(data)[myplot$tree_row$order],
	data[myplot$tree_row$order, ]), 
	file = "tsg.GTRD.peaks.tsv", 
	sep = "\t", row.names = FALSE, quote = FALSE)
save.image("GTRD.RData")
