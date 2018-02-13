library(data.table)
library(pheatmap)
args <- c("GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct", "GTEx_sample.tissue.txt")

### some functions ###
###+++###	
#Function require a vector with expression of one gene in different tissues/samples.
#Mean is calculated taking in account tissues with 0 expression. 2+0+4=2
fmean <- function(x)
	{
		if(!all(is.na(x)))
	 	{
	 		res <- mean(x, na.rm = TRUE)
	 	} else {
	 		res <- NA
	 	}
	 	return(res)
	}
###***###***###	

###+++###	
#Function require a vector with expression of one gene in different tissues/samples.
#Max is calculated taking in account tissues with 0 expression. 2+0+4=2
fmax <- function(x)
	{
		if(!all(is.na(x)))
	 	{
	 		res <- max(x, na.rm = TRUE)
	 	} else {
	 		res <- NA
	 	}
	 	return(res)
	}
###***###***###	

###+++###
#Function require a vector with expression of one gene in different tissues/samples.
#If expression for one tissue is not known, gene specificity for this gene is NA
#Minimum 2 tissues/samples
fTau <- function(x)
{
	if(all(!is.na(x)))
 	{
 		if(min(x, na.rm = TRUE) >= 0)
		{
 			mx <- fmax(x)
 			if(mx != 0 && !is.na(mx))
 			{
 				x <- (1-(x/mx))
 				res <- sum(x, na.rm = TRUE)
 				res <- res/(length(x)-1)
 			} else {
 				res <- 0
 			}
 		} else {
 		res <- NA
 		#print("Expression values have to be positive!")
 		} 
 	} else {
 		res <- NA
 		#print("No data for this gene avalable.")
 	} 
 	return(res)
}
###***###***###

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
print("Calculating nullcount, mean and Tau for each tissue")
nullcount <- matrix(ncol = length(tissues), nrow = nrow(tpm))
log2tpm.mean <- matrix(ncol = length(tissues), nrow = nrow(tpm))
log2tpm.Tau <- matrix(ncol = length(tissues), nrow = nrow(tpm))
for(i in 1:length(tissues)){
	z <- log2(tpm[, info[,2] == tissues[i]] + 1)
	nullcount[,i] <- apply(z, 1, function(x) {length(x[x==0])})
	log2tpm.mean[,i] <- apply(z, 1, fmean)
	log2tpm.Tau[,i] <- apply(z, 1, fTau)
}
colnames(nullcount) <- tissues
colnames(log2tpm.mean) <- tissues
colnames(log2tpm.Tau) <- tissues

print("Calculating nullcount, mean and Tau across tissues")
log2tpm.mean.mean <- apply(log2tpm.mean, 1, fmean)
log2tpm.mean.Tau <- apply(log2tpm.mean, 1, fTau)

nullcount.sum <- apply(nullcount, 1, sum)
log2tpm.Tau.max <- apply(log2tpm.Tau, 1, function(x) {max(x[is.finite(x)], na.rm = T)})
log2tpm.Tau.min <- apply(log2tpm.Tau, 1, function(x) {min(x[is.finite(x)], na.rm = T)})

print("Saving RData")
save.image("v1.3.Tau.RData")

print("Looking for housekeeping genes")
all.stats <- data.frame(rownames(tpm), log2tpm.mean.mean, log2tpm.mean.Tau, 
	nullcount.sum, log2tpm.Tau.max, log2tpm.Tau.min, 
	log2tpm.mean, log2tpm.Tau)
rownames(all.stats) <- rownames(tpm)
housekeepinggene <- all.stats[nullcount.sum == 0 & !is.na(log2tpm.Tau.max) & log2tpm.Tau.max < 0.15 & !is.na(log2tpm.mean.Tau) & log2tpm.mean.Tau < 0.15, ]
housekeepinggene.mean <- log2tpm.mean[nullcount.sum == 0 & !is.na(log2tpm.Tau.max) & log2tpm.Tau.max < 0.15 & !is.na(log2tpm.mean.Tau) & log2tpm.mean.Tau < 0.15, ]
housekeepinggene.Tau <- log2tpm.Tau[nullcount.sum == 0 & !is.na(log2tpm.Tau.max) & log2tpm.Tau.max < 0.15 & !is.na(log2tpm.mean.Tau) & log2tpm.mean.Tau < 0.15, ]
rownames(housekeepinggene.mean) <- housekeepinggene[,1]
rownames(housekeepinggene.Tau) <- housekeepinggene[,1]
# housekeepinggene[grep("GAPDH", as.matrix(housekeepinggene[,1])), 1]

## check ENSG00000204983.8|PRSS1 Pancreas
# tissuespecificgene[grep("PRSS1",as.matrix(tissuespecificgene[,1])),1]
# raw.table[grep("PRSS1",as.matrix(raw.table[,1])),1:5][10,]
print("Looking for tissue specific genes")
tissuespecificgene <- data.frame()
tissueflags <- c()
rindex <- c()
for(i in 1:nrow(all.stats)){
	if(!is.na(log2tpm.mean.Tau[i]) && log2tpm.mean.Tau[i] > 0.85 && !is.na(log2tpm.Tau.min[i]) && log2tpm.Tau.min[i] < 0.15 && min(nullcount[i,], na.rm = T) == 0){
		flag <- 0
		tflags <- c()
		for(j in 1:length(tissues)){
			temp <- log2tpm.mean[i,j]/fmax(log2tpm.mean[i,])
			if(nullcount[i,j] == 0 && !is.na(log2tpm.Tau[i,j]) && log2tpm.Tau[i,j] < 0.15 &&  !is.na(temp) && temp > 0.85){
				flag <- 1
				tflags <- c(tflags, tissues[j])
			}
		}
		if(flag == 1){
			tissuespecificgene <- rbind.data.frame(tissuespecificgene, all.stats[i,])
			tissueflags <- c(tissueflags, paste(tflags, collapse = ","))
			rindex <- c(rindex, i)
		}
	}
}
tissuespecificgene.mean <- log2tpm.mean[rindex, ]
tissuespecificgene.Tau <- log2tpm.Tau[rindex, ]
rownames(tissuespecificgene.mean) <- tissuespecificgene[,1]
rownames(tissuespecificgene.Tau) <- tissuespecificgene[,1]

print("Writing output files")
write.table(housekeepinggene, file = "HKG.v1.3.tsv", 
	sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(tissuespecificgene, tissueflags), file = "TSG.v1.3.tsv", 
	sep = "\t", row.names = FALSE, quote = FALSE)

write.table(data.frame(rownames(tpm), nullcount), file = "v1.3.nullcount.tsv", 
	sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(rownames(tpm), log2tpm.mean), file = "v1.3.log2tpm.mean.tsv", 
	sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(rownames(tpm), log2tpm.Tau), file = "v1.3.log2tpm.Tau.tsv", 
	sep = "\t", row.names = FALSE, quote = FALSE)
write.table(all.stats, file = "v1.3.allstats.tsv", 
	sep = "\t", row.names = FALSE, quote = FALSE)

print("Clustering")
breaklists <- c(seq(0, 1, by = 0.01))
colorn <- length(breaklists)
colors <- colorRampPalette(c("blue", "yellow", "red"))(colorn)

breaklists2 <- c(seq(0, 6, by = 0.05),seq(6.1, 16, by = 0.1))
colorn2 <- length(breaklists2)
colors2 <- colorRampPalette(c("blue", "yellow", "red"))(colorn2)

png(filename = "HKG.v1.3.mean.png", width = 1000, height = 1500)
data <- housekeepinggene.mean
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
	file = "HKG.v1.3.mean.pheatmap.tsv", 
	sep = "\t", row.names = FALSE, quote = FALSE)

png(filename = "HKG.v1.3.Tau.png", width = 1000, height = 1500)
data <- housekeepinggene.Tau
pheatmap(data[p1$tree_row$order,p1$tree_col$order], scale = "none", 
	show_rownames = F, show_colnames = T, color = colors, 
	cluster_cols = F, cluster_rows = F, breaks = breaklists)
dev.off()

png(filename = "TSG.v1.3.mean.png", width = 1000, height = 1500)
data <- tissuespecificgene.mean
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
	file = "TSG.v1.3.mean.pheatmap.tsv", 
	sep = "\t", row.names = FALSE, quote = FALSE)

png(filename = "TSG.v1.3.Tau.png", width = 1000, height = 1500)
data <- housekeepinggene.Tau
pheatmap(data[p2$tree_row$order,p2$tree_col$order], scale = "none", 
	show_rownames = F, show_colnames = T, color = colors, 
	cluster_cols = F, cluster_rows = F, breaks = breaklists)
dev.off()


