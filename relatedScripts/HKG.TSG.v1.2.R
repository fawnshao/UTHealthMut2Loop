library(data.table)
library(pheatmap)
args <- c("GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct", "GTEx_sample.tissue.txt")

## read in data
# tpm <- as.matrix(read.table(args[1], sep = "\t", row.names = 1, header = T))
raw.table <- fread(args[1], sep = "\t", header = T)
tpm <- as.matrix(raw.table[,-1])
rownames(tpm) <- as.matrix(raw.table[,1])
info <- as.matrix(read.table(args[2], sep = "\t"))
tissues <- unique(info[,2])
tissues <- tissues[-c(7,24,25,53)]

## calculate
nullcount <- matrix(ncol = length(tissues), nrow = nrow(tpm))
log2tpm.median <- matrix(ncol = length(tissues), nrow = nrow(tpm))
log2tpm.iqr <- matrix(ncol = length(tissues), nrow = nrow(tpm))
log2tpm.tsscore <- matrix(ncol = length(tissues), nrow = nrow(tpm))
for(i in 1:length(tissues)){
	z <- log2(tpm[, info[,2] == tissues[i]] + 1)
	nullcount[,i] <- apply(z, 1, function(x){length(x[x==0])})
	log2tpm.median[,i] <- apply(z, 1, function(x){quantile(x, probs = 0.5, na.rm = T)})
	log2tpm.iqr[,i] <- apply(z, 1, function(x){quantile(x, probs = 0.75, na.rm = T) - quantile(x, probs = 0.25, na.rm = T)})
	log2tpm.tsscore[,i] <- apply(abs(z - log2tpm.median[,i]) / log2tpm.iqr[,i], 1, function(x){max(x[is.finite(x)], na.rm = T)})
}
log2tpm.varscore <- log2tpm.iqr / log2tpm.median
colnames(nullcount) <- tissues
colnames(log2tpm.median) <- tissues
colnames(log2tpm.iqr) <- tissues
colnames(log2tpm.tsscore) <- tissues
colnames(log2tpm.varscore) <- tissues

log2tpm.median.median <- apply(log2tpm.median, 1, function(x){quantile(x, probs = 0.5, na.rm = T)})
log2tpm.median.iqr <- apply(log2tpm.median, 1, function(x){quantile(x, probs = 0.75, na.rm = T) - quantile(x, probs = 0.25, na.rm = T)})
log2tpm.median.tsscore <- apply(abs(z - log2tpm.median[,i]) / log2tpm.iqr[,i], 1, function(x){max(x[is.finite(x)], na.rm = T)})
log2tpm.median.varscore <- log2tpm.median.iqr / log2tpm.median.median
nullcount.sum <- apply(nullcount, 1, sum)
log2tpm.varscore.max <- apply(log2tpm.varscore, 1, function(x) {max(x[is.finite(x)], na.rm = T)})
log2tpm.tsscore.max <- apply(log2tpm.tsscore, 1, function(x) {max(x[is.finite(x)], na.rm = T)})

all.stats <- data.frame(rownames(tpm), log2tpm.median.median, log2tpm.median.iqr, 
	log2tpm.median.tsscore, log2tpm.median.varscore,
	nullcount.sum, log2tpm.tsscore.max, log2tpm.varscore.max)
rownames(all.stats) <- rownames(tpm)
housekeepinggene <- all.stats[nullcount.sum == 0 & is.finite(log2tpm.median.varscore) & log2tpm.median.varscore < 0.2 & is.finite(log2tpm.median.tsscore) & log2tpm.median.tsscore < 2 & is.finite(log2tpm.varscore.max) & log2tpm.varscore.max < 0.2 & is.finite(log2tpm.tsscore.max) & log2tpm.tsscore.max < 2, ]
housekeepinggene[grep("GAPDH", as.matrix(housekeepinggene[,1])),1]

# tissuespecificgene <- raw.table[tissues.mean.dif.max > 2 & is.finite(sd.percent.min) & sd.percent.min < 0.2 & tissues.largefc.percent.min < 0.01, ]
## check ENSG00000204983.8|PRSS1 Pancreas
# tissuespecificgene[grep("PRSS1",as.matrix(tissuespecificgene[,1])),1]
# raw.table[grep("PRSS1",as.matrix(raw.table[,1])),1:5][10,]
tissuespecificgene <- data.frame()
tissueflags <- c()
# a <- raw.table[tissues.mean.dif.max > 2 & is.finite(all.sd.percent) & all.sd.percent > 0.5 & all.stats[,3] > 0 & is.finite(sd.percent.max) & sd.percent.max > 1, ]
# should be sd.percent.min
for(i in 1:nrow(all.stats)){
	if(is.finite(log2tpm.median.varscore[i]) & log2tpm.median.varscore[i] > 0.5 & is.finite(log2tpm.median.tsscore[i]) & log2tpm.median.tsscore[i] > 2){
		flag <- 0
		tflags <- c()
		for(j in 1:length(tissues)){
			if(nullcount[i,j]==0 & is.finite(log2tpm.varscore[i,j]) & log2tpm.varscore[i,j] < 0.2 & is.finite(log2tpm.tsscore[i,j]) & log2tpm.tsscore[i,j] < 2){
				flag <- 1
				tflags <- c(tflags, tissues[j])
			}
		}
		if(flag == 1){
			tissuespecificgene <- rbind.data.frame(tissuespecificgene, all.stats[i,])
			tissueflags <- c(tissueflags, paste(tflags, collapse=","))
		}
	}
}

write.table(housekeepinggene, file = "HKG.v1.2.by.percentage.tsv", 
	sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(tissueflags, tissuespecificgene), file = "TSG.v1.2.by.percentage.tsv", 
	sep = "\t", row.names = FALSE, quote = FALSE)

write.table(data.frame(rownames(tpm), nullcount), file = "v1.2.nullcount.tsv", 
	sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(rownames(tpm), log2tpm.median), file = "v1.2.log2tpm.median.tsv", 
	sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(rownames(tpm), log2tpm.iqr), file = "v1.2.log2tpm.iqr.tsv", 
	sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(rownames(tpm), log2tpm.varscore), file = "v1.2.log2tpm.varscore.tsv", 
	sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(rownames(tpm), log2tpm.tsscore), file = "v1.2.log2tpm.tsscore.tsv", 
	sep = "\t", row.names = FALSE, quote = FALSE)

breaklists <- c(seq(0, 2, by = 0.01),seq(2.1, 5.7, by = 0.1))
colorn <- length(breaklists)
colors <- colorRampPalette(c("blue", "yellow", "red"))(colorn)

breaklists2 <- c(seq(0, 6, by = 0.05),seq(6.1, 16, by = 0.1))
colorn2 <- length(breaklists2)
colors2 <- colorRampPalette(c("blue", "yellow", "red"))(colorn2)

png(filename = "HKG.v1.2.by.percentage.median.png", width = 1000, height = 1500)
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
	file = "HKG.v1.2.by.percentage.median.pheatmap.tsv", 
	sep = "\t", row.names = FALSE, quote = FALSE)

png(filename = "HKG.v1.2.by.percentage.iqr.png", width = 1000, height = 1500)
data <- hkg.sd
pheatmap(data[p1$tree_row$order,p1$tree_col$order], scale = "none", 
	show_rownames = F, show_colnames = T, color = colors, 
	cluster_cols = F, cluster_rows = F, breaks = breaklists)
dev.off()

png(filename = "TSG.v1.2.by.percentage.median.png", width = 1000, height = 1500)
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

png(filename = "TSG.v1.2.by.percentage.iqr.png", width = 1000, height = 1500)
data <- tsg.sd
pheatmap(data[p2$tree_row$order,p2$tree_col$order], scale = "none", 
	show_rownames = F, show_colnames = T, color = colors, 
	cluster_cols = F, cluster_rows = F, breaks = breaklists)
dev.off()
save.image("v1.2.iqr.RData")

