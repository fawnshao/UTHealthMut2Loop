library(data.table)
library(pheatmap)
args <- c("housekeepinggene.byGTEx.v1.1.raw.tsv", "GTEx_sample.tissue.txt")

## read in data
# tpm <- as.matrix(read.table(args[1], sep = "\t", row.names = 1, header = T))
raw.table <- fread(args[1], sep = "\t", header = T)

info <- as.matrix(read.table(args[2], sep = "\t"))
tissues <- unique(info[,2])
tissues <- tissues[-c(7,24,25,53)]
tissues.count <- table(info[,2])
tissues.count <- tissues.count[-c(7,24,25,53)]

all.stats <- as.matrix(raw.table[,2:5])
tissues.mean <- as.matrix(raw.table[,6:54])
tissues.sd <- as.matrix(raw.table[,55:103])
tissues.largefc <- as.matrix(raw.table[,104:152])
tissues.null <- as.matrix(raw.table[,153:201])

rownames(tissues.mean) <- as.matrix(raw.table[,1])
rownames(tissues.sd) <- as.matrix(raw.table[,1])
rownames(tissues.largefc) <- as.matrix(raw.table[,1])
rownames(tissues.null) <- as.matrix(raw.table[,1])
colnames(tissues.mean) <- tissues
colnames(tissues.sd) <- tissues
colnames(tissues.largefc) <- tissues
colnames(tissues.null) <- tissues

sd.percent <- tissues.sd/tissues.mean
sd.percent.max <- apply(sd.percent, 1, function(x){max(x[is.finite(x)], na.rm = T)})
sd.percent.min <- apply(sd.percent, 1, function(x){min(x[is.finite(x)], na.rm = T)})
length(sd.percent.max[is.finite(sd.percent.max) & sd.percent.max < 0.2])

tissues.largefc.max <- apply(tissues.largefc, 1, function(x){max(x[is.finite(x)], na.rm = T)})
tissues.largefc.percent <- t(t(tissues.largefc)/as.vector(tissues.count))
tissues.largefc.percent.max <- apply(tissues.largefc.percent, 1, function(x){max(x[is.finite(x)], na.rm = T)})
tissues.largefc.percent.min <- apply(tissues.largefc.percent, 1, function(x){min(x[is.finite(x)], na.rm = T)})

all.sd.percent <- all.stats[,2]/all.stats[,1]
length(all.sd.percent[is.finite(all.sd.percent) & all.sd.percent < 0.2])

# housekeepinggene <- raw.table[all.stats[,4] == 0 & is.finite(all.sd.percent) & all.sd.percent < 0.2 & all.stats[,3] == 0 & is.finite(sd.percent.max) & sd.percent.max < 0.2 & tissues.largefc.max < 10, ]
housekeepinggene <- raw.table[all.stats[,4] == 0 & is.finite(all.sd.percent) & all.sd.percent < 0.2 & all.stats[,3] == 0 & is.finite(sd.percent.max) & sd.percent.max < 0.2 & tissues.largefc.percent.max < 0.01, ]
housekeepinggene[grep("GAPDH", as.matrix(housekeepinggene[,1])),1]

tissues.mean.min <- apply(tissues.mean, 1, function(x){min(x[is.finite(x)], na.rm = T)})
tissues.mean.max <- apply(tissues.mean, 1, function(x){max(x[is.finite(x)], na.rm = T)})
tissues.mean.dif.max <- tissues.mean.max - tissues.mean.min
tissues.mean.dif <- tissues.mean - tissues.mean.min
# tissuespecificgene <- raw.table[tissues.mean.dif.max > 2 & is.finite(sd.percent.min) & sd.percent.min < 0.2 & tissues.largefc.percent.min < 0.01, ]
## check ENSG00000204983.8|PRSS1 Pancreas
# tissuespecificgene[grep("PRSS1",as.matrix(tissuespecificgene[,1])),1]
# raw.table[grep("PRSS1",as.matrix(raw.table[,1])),1:5][10,]
tissuespecificgene <- data.frame()
tissueflags <- c()
# a <- raw.table[tissues.mean.dif.max > 2 & is.finite(all.sd.percent) & all.sd.percent > 0.5 & all.stats[,3] > 0 & is.finite(sd.percent.max) & sd.percent.max > 1, ]
# should be sd.percent.min
for(i in 1:nrow(raw.table)){
	if(tissues.mean.dif.max[i] > 2 && is.finite(all.sd.percent[i]) && all.sd.percent[i] > 0.5 && all.stats[i,3] > 0 && is.finite(sd.percent.min[i]) && sd.percent.min[i] < 0.2){
		flag <- 0
		tflags <- c()
		for(j in 1:length(tissues)){
			if(sd.percent[i,j] < 0.2 && tissues.largefc.percent[i,j] < 0.01 && tissues.mean.dif[i,j] > 2){
				flag <- 1
				tflags <- c(tflags, tissues[j])
			}
		}
		if(flag == 1){
			tissuespecificgene <- rbind.data.frame(tissuespecificgene, raw.table[i,])
			tissueflags <- c(tissueflags, paste(tflags, collapse=","))
		}
	}
}

write.table(housekeepinggene, file = "HKG.v1.2.by.percentage.tsv", 
	sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(tissueflags, tissuespecificgene), file = "TSG.v1.2.by.percentage.tsv", 
	sep = "\t", row.names = FALSE, quote = FALSE)

hkg.mean <- as.matrix(housekeepinggene[,6:54])
tsg.mean <- as.matrix(tissuespecificgene[,6:54])
hkg.sd <- as.matrix(housekeepinggene[,55:103])
tsg.sd <- as.matrix(tissuespecificgene[,55:103])
rownames(hkg.mean) <- as.matrix(housekeepinggene[,1])
colnames(hkg.mean) <- tissues
rownames(tsg.mean) <- as.matrix(tissuespecificgene[,1])
colnames(tsg.mean) <- tissues
rownames(hkg.sd) <- as.matrix(housekeepinggene[,1])
colnames(hkg.sd) <- tissues
rownames(tsg.sd) <- as.matrix(tissuespecificgene[,1])
colnames(tsg.sd) <- tissues

breaklists <- c(seq(0, 2, by = 0.01),seq(2.1, 5.7, by = 0.1))
colorn <- length(breaklists)
colors <- colorRampPalette(c("blue", "yellow", "red"))(colorn)

breaklists2 <- c(seq(0, 6, by = 0.05),seq(6.1, 16, by = 0.1))
colorn2 <- length(breaklists2)
colors2 <- colorRampPalette(c("blue", "yellow", "red"))(colorn2)

png(filename = "HKG.v1.2.by.percentage.mean.png", width = 1000, height = 1500)
data <- hkg.mean
p1 <- pheatmap(data, scale = "none", show_rownames = F, show_colnames = T, 
         color = colors2, cluster_cols = T, breaks = breaklists2,
         clustering_distance_cols = "euclidean", clustering_distance_rows = "euclidean", 
         clustering_method = "ward.D2"
         )
dev.off()
cluster <- cutree(p1$tree_row, k = 10)
write.table(data.frame(cluster[p1$tree_row$order], 
	rownames(data)[p1$tree_row$order],
	data[p1$tree_row$order, p1$tree_col$order]), 
	file = "HKG.by.percentage.mean.pheatmap.tsv", 
	sep = "\t", row.names = FALSE, quote = FALSE)

png(filename = "HKG.v1.2.by.percentage.sd.png", width = 1000, height = 1500)
data <- hkg.sd
pheatmap(data[p1$tree_row$order,p1$tree_col$order], scale = "none", 
	show_rownames = F, show_colnames = T, color = colors, 
	cluster_cols = F, cluster_rows = F, breaks = breaklists)
dev.off()

png(filename = "TSG.v1.2.by.percentage.mean.png", width = 1000, height = 1500)
data <- tsg.mean
p2 <- pheatmap(data, scale = "none", show_rownames = F, show_colnames = T, 
         color = colors2, cluster_cols = T, breaks = breaklists2,
         clustering_distance_cols = "euclidean", clustering_distance_rows = "euclidean", 
         clustering_method = "ward.D2"
         )
dev.off()
cluster <- cutree(p2$tree_row, k = 10)
write.table(data.frame(cluster[p2$tree_row$order], 
	tissueflags[p2$tree_row$order],
	rownames(data)[p2$tree_row$order],
	data[p2$tree_row$order, p2$tree_col$order]), 
	file = "TSG.by.percentage.mean.pheatmap.tsv", 
	sep = "\t", row.names = FALSE, quote = FALSE)

png(filename = "TSG.v1.2.by.percentage.sd.png", width = 1000, height = 1500)
data <- tsg.sd
pheatmap(data[p2$tree_row$order,p2$tree_col$order], scale = "none", 
	show_rownames = F, show_colnames = T, color = colors, 
	cluster_cols = F, cluster_rows = F, breaks = breaklists)
dev.off()


