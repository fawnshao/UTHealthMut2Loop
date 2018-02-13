library(data.table)
library(pheatmap)
args <- c("GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct", "GTEx_sample.tissue.txt")

## read in data
# tpm <- as.matrix(read.table(args[1], sep = "\t", row.names = 1, header = T))
print("Reading Data")
raw.table <- fread(args[1], sep = "\t", header = T)
tpm <- as.matrix(raw.table[,-1])
rownames(tpm) <- as.matrix(raw.table[,1])
info <- as.matrix(read.table(args[2], sep = "\t"))
tissues <- unique(info[,2])
tissues <- tissues[-c(7,24,25,53)]

## calculate
print("Calculating nullcount, median, IQR and Score for each tissue")
nullcount <- matrix(ncol = length(tissues), nrow = nrow(tpm))
log2tpm.median <- matrix(ncol = length(tissues), nrow = nrow(tpm))
log2tpm.iqr <- matrix(ncol = length(tissues), nrow = nrow(tpm))
log2tpm.Sscore <- matrix(ncol = length(tissues), nrow = nrow(tpm))
for(i in 1:length(tissues)){
	z <- log2(tpm[, info[,2] == tissues[i]] + 1)
	nullcount[,i] <- apply(z, 1, function(x) {length(x[x==0])})
	log2tpm.median[,i] <- apply(z, 1, function(x) {quantile(x, probs = 0.5, na.rm = T)})
	log2tpm.iqr[,i] <- apply(z, 1, function(x) {quantile(x, probs = 0.75, na.rm = T) - quantile(x, probs = 0.25, na.rm = T)})
	log2tpm.Sscore[,i] <- apply(abs(z - log2tpm.median[,i]) / log2tpm.iqr[,i], 1, function(x){quantile(x[is.finite(x)], probs = 0.5, na.rm = T)})
}
log2tpm.Vscore <- log2tpm.iqr / log2tpm.median
colnames(nullcount) <- tissues
colnames(log2tpm.median) <- tissues
colnames(log2tpm.iqr) <- tissues
colnames(log2tpm.Sscore) <- tissues
colnames(log2tpm.Vscore) <- tissues

print("Calculating nullcount, median, IQR and across tissues")
log2tpm.median.median <- apply(log2tpm.median, 1, function(x) {quantile(x[is.finite(x)], probs = 0.5, na.rm = T)})
log2tpm.median.iqr <- apply(log2tpm.median, 1, function(x) {quantile(x[is.finite(x)], probs = 0.75, na.rm = T) - quantile(x[is.finite(x)], probs = 0.25, na.rm = T)})

log2tpm.median.tmp.score <- (log2tpm.median - log2tpm.median.median) / log2tpm.median.iqr
log2tpm.median.Sscore <- apply(log2tpm.median.tmp.score, 1, function(x) {max(x[is.finite(x)], na.rm = T)})
log2tpm.median.Sprimescore <- apply(log2tpm.median.tmp.score, 1, function(x) {min(x[is.finite(x)], na.rm = T)})

log2tpm.median.Vscore <- log2tpm.median.iqr / log2tpm.median.median
nullcount.sum <- apply(nullcount, 1, sum)
log2tpm.Vscore.max <- apply(log2tpm.Vscore, 1, function(x) {max(x[is.finite(x)], na.rm = T)})
log2tpm.Sscore.max <- apply(log2tpm.Sscore, 1, function(x) {max(x[is.finite(x)], na.rm = T)})

print("Looking for housekeeping genes")
all.stats <- data.frame(rownames(tpm), log2tpm.median.median, log2tpm.median.iqr, 
	log2tpm.median.Sscore, log2tpm.median.Sprimescore, log2tpm.median.Vscore,
	nullcount.sum, log2tpm.Sscore.max, log2tpm.Vscore.max)
rownames(all.stats) <- rownames(tpm)
housekeepinggene <- all.stats[nullcount.sum == 0 & is.finite(log2tpm.Vscore.max) & log2tpm.Vscore.max < 0.2 & is.finite(log2tpm.Sscore.max) & log2tpm.Sscore.max < 2 & is.finite(log2tpm.median.Vscore) & log2tpm.median.Vscore < 0.2 & is.finite(log2tpm.median.Sscore) & log2tpm.median.Sscore < 2 & is.finite(log2tpm.median.Sprimescore) & log2tpm.median.Sprimescore > -2, ]
housekeepinggene.median <- log2tpm.median[nullcount.sum == 0 & is.finite(log2tpm.Vscore.max) & log2tpm.Vscore.max < 0.2 & is.finite(log2tpm.Sscore.max) & log2tpm.Sscore.max < 2 & is.finite(log2tpm.median.Vscore) & log2tpm.median.Vscore < 0.2 & is.finite(log2tpm.median.Sscore) & log2tpm.median.Sscore < 2 & is.finite(log2tpm.median.Sprimescore) & log2tpm.median.Sprimescore > -2, ]
housekeepinggene.iqr <- log2tpm.iqr[nullcount.sum == 0 & is.finite(log2tpm.Vscore.max) & log2tpm.Vscore.max < 0.2 & is.finite(log2tpm.Sscore.max) & log2tpm.Sscore.max < 2 & is.finite(log2tpm.median.Vscore) & log2tpm.median.Vscore < 0.2 & is.finite(log2tpm.median.Sscore) & log2tpm.median.Sscore < 2 & is.finite(log2tpm.median.Sprimescore) & log2tpm.median.Sprimescore > -2, ]
# housekeepinggene[grep("GAPDH", as.matrix(housekeepinggene[,1])), 1]

## check ENSG00000204983.8|PRSS1 Pancreas
# tissuespecificgene[grep("PRSS1",as.matrix(tissuespecificgene[,1])),1]
# raw.table[grep("PRSS1",as.matrix(raw.table[,1])),1:5][10,]
print("Looking for tissue specific genes")
tissuespecificgene <- data.frame()
tissuehigh <- c()
tissuelow <- c()
rindex <- c()
for(i in 1:nrow(all.stats)){
	if((is.finite(log2tpm.median.Vscore[i]) && log2tpm.median.Vscore[i] > 0.5) && ( (is.finite(log2tpm.median.Sscore[i]) && log2tpm.median.Sscore[i] > 2) || (is.finite(log2tpm.median.Sprimescore[i]) && log2tpm.median.Sprimescore[i] < -2) )){
		flag <- 0
		hflags <- c()
		lflags <- c()
		for(j in 1:length(tissues)){
			if(nullcount[i,j] == 0 && is.finite(log2tpm.Vscore[i,j]) && log2tpm.Vscore[i,j] < 0.2 && is.finite(log2tpm.Sscore[i,j]) && log2tpm.Sscore[i,j] < 2 && is.finite(log2tpm.median.tmp.score[i,j]) && log2tpm.median.tmp.score[i,j] > 2){
				flag <- 1
				hflags <- c(hflags, tissues[j])
			}
			if(nullcount[i,j] == 0 && is.finite(log2tpm.Vscore[i,j]) && log2tpm.Vscore[i,j] < 0.2 && is.finite(log2tpm.Sscore[i,j]) && log2tpm.Sscore[i,j] < 2 && is.finite(log2tpm.median.tmp.score[i,j]) && log2tpm.median.tmp.score[i,j] < -2){
				flag <- 1
				lflags <- c(lflags, tissues[j])
			}
		}
		if(flag == 1){
			tissuespecificgene <- rbind.data.frame(tissuespecificgene, all.stats[i,])
			tissuehigh <- c(tissuehigh, paste(hflags, collapse = ","))
			tissuelow <- c(tissuelow, paste(lflags, collapse = ","))
			rindex <- c(rindex, i)
		}
	}
}
tissuespecificgene.median <- log2tpm.median[rindex, ]
tissuespecificgene.iqr <- log2tpm.iqr[rindex, ]

print("Writing output files")
write.table(housekeepinggene, file = "HKG.v1.2.by.percentage.tsv", 
	sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(tissuespecificgene, tissuehigh, tissuelow), file = "TSG.v1.2.by.percentage.tsv", 
	sep = "\t", row.names = FALSE, quote = FALSE)

write.table(data.frame(rownames(tpm), nullcount), file = "v1.2.nullcount.tsv", 
	sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(rownames(tpm), log2tpm.median), file = "v1.2.log2tpm.median.tsv", 
	sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(rownames(tpm), log2tpm.iqr), file = "v1.2.log2tpm.iqr.tsv", 
	sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(rownames(tpm), log2tpm.Vscore), file = "v1.2.log2tpm.Vscore.tsv", 
	sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(rownames(tpm), log2tpm.Sscore), file = "v1.2.log2tpm.Sscore.tsv", 
	sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(rownames(tpm), all.stats), file = "v1.2.allstats.tsv", 
	sep = "\t", row.names = FALSE, quote = FALSE)

print("Saving RData")
save.image("v1.2.iqr.RData")

print("Clustering")
breaklists <- c(seq(0, 2, by = 0.01),seq(2.1, 5.7, by = 0.1))
colorn <- length(breaklists)
colors <- colorRampPalette(c("blue", "yellow", "red"))(colorn)

breaklists2 <- c(seq(0, 6, by = 0.05),seq(6.1, 16, by = 0.1))
colorn2 <- length(breaklists2)
colors2 <- colorRampPalette(c("blue", "yellow", "red"))(colorn2)

png(filename = "HKG.v1.2.median.png", width = 1000, height = 1500)
data <- housekeepinggene.median
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
	file = "HKG.v1.2.median.pheatmap.tsv", 
	sep = "\t", row.names = FALSE, quote = FALSE)

png(filename = "HKG.v1.2.iqr.png", width = 1000, height = 1500)
data <- housekeepinggene.iqr
pheatmap(data[p1$tree_row$order,p1$tree_col$order], scale = "none", 
	show_rownames = F, show_colnames = T, color = colors, 
	cluster_cols = F, cluster_rows = F, breaks = breaklists)
dev.off()

png(filename = "TSG.v1.2.median.png", width = 1000, height = 1500)
data <- tissuespecificgene.median
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
	file = "TSG.v1.2.median.pheatmap.tsv", 
	sep = "\t", row.names = FALSE, quote = FALSE)

png(filename = "TSG.v1.2.iqr.png", width = 1000, height = 1500)
data <- housekeepinggene.iqr
pheatmap(data[p2$tree_row$order,p2$tree_col$order], scale = "none", 
	show_rownames = F, show_colnames = T, color = colors, 
	cluster_cols = F, cluster_rows = F, breaks = breaklists)
dev.off()


