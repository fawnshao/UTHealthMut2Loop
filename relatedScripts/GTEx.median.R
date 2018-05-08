library(data.table)
library(pheatmap)

args <- c("GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct", 
	"GTEx_sample.tissue.txt", "v1.5", "0.25", "0.5", "0.1", "3")
# expression tissuename outputpre tissueTau sampleTau expressioncutoff[log2(tpm+1)]

### some functions ###
###+++###
fmedian <- function(x)
	{
		if(!all(is.na(x)))
	 	{
	 		res <- median(x, na.rm = TRUE)
	 	} else {
	 		res <- NA
	 	}
	 	return(res)
	}
###***###***###	

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
# Read 56202 rows and 11689 (of 11689) columns from 2.622 GB file in 00:06:16
tpm <- as.matrix(raw.table[,-1])
rownames(tpm) <- as.matrix(raw.table[,1])
info <- as.matrix(read.table(args[2], sep = "\t"))
tissues <- unique(info[,2])
tissues <- tissues[-c(7,24,25,53)]

outputpre <- args[3]
tissueTau <- as.numeric(args[4])
sampleTau <- as.numeric(args[5])
nullexpression <- as.numeric(args[6])
############
medianthreshold <- as.numeric(args[7])
############

## calculate
print("Calculating nullcount, median/mean and Tau for each tissue")
nullcount <- matrix(ncol = length(tissues), nrow = nrow(tpm))
log2tpm.median <- matrix(ncol = length(tissues), nrow = nrow(tpm))
log2tpm.mean <- matrix(ncol = length(tissues), nrow = nrow(tpm))
log2tpm.Tau <- matrix(ncol = length(tissues), nrow = nrow(tpm))
for(i in 1:length(tissues)){
	z <- log2(tpm[, info[,2] == tissues[i]] + 1)
	nullcount[,i] <- apply(z, 1, function(x) {length(x[x < nullexpression])})
	log2tpm.median[,i] <- apply(z, 1, fmedian)
	log2tpm.mean[,i] <- apply(z, 1, fmean)
	log2tpm.Tau[,i] <- apply(z, 1, fTau)
}
colnames(nullcount) <- tissues
colnames(log2tpm.median) <- tissues
colnames(log2tpm.mean) <- tissues
colnames(log2tpm.Tau) <- tissues

print("Calculating nullcount, mean and Tau across tissues")
log2tpm.median.mean <- apply(log2tpm.median, 1, fmean)
log2tpm.median.Tau <- apply(log2tpm.median, 1, fTau)
log2tpm.mean.mean <- apply(log2tpm.mean, 1, fmean)
log2tpm.mean.Tau <- apply(log2tpm.mean, 1, fTau)

nullcount.sum <- apply(nullcount, 1, sum)
log2tpm.Tau.max <- apply(log2tpm.Tau, 1, function(x) {max(x[is.finite(x)], na.rm = T)})
log2tpm.Tau.min <- apply(log2tpm.Tau, 1, function(x) {min(x[is.finite(x)], na.rm = T)})

### it takes ~30min until this step

# print("Saving RData")
# save.image("v1.3.Tau.RData")
# ### it takes ~50min until this step

print("Writing raw output files")
all.stats <- data.frame(rownames(tpm), log2tpm.mean.mean, log2tpm.mean.Tau, 
	log2tpm.median.mean, log2tpm.median.Tau, 
	nullcount.sum, log2tpm.Tau.max, log2tpm.Tau.min, 
	log2tpm.mean, log2tpm.median, log2tpm.Tau)
rownames(all.stats) <- rownames(tpm)
write.table(data.frame(rownames(tpm), nullcount), file = paste(outputpre, "nullcount.tsv", sep = "."), 
	sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(rownames(tpm), log2tpm.mean), file = paste(outputpre, "log2tpm.mean.tsv", sep = "."), 
	sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(rownames(tpm), log2tpm.median), file = paste(outputpre, "log2tpm.median.tsv", sep = "."), 
	sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(rownames(tpm), log2tpm.Tau), file = paste(outputpre, "log2tpm.Tau.tsv", sep = "."), 
	sep = "\t", row.names = FALSE, quote = FALSE)
write.table(all.stats, file = paste(outputpre, "allstats.tsv", sep = "."), 
	sep = "\t", row.names = FALSE, quote = FALSE)

#load("v1.3.Tau.RData")
#library(pheatmap)
print("Looking for housekeeping genes")
housekeepinggene <- all.stats[nullcount.sum == 0 & !is.na(log2tpm.Tau.max) & log2tpm.Tau.max < sampleTau & !is.na(log2tpm.median.Tau) & log2tpm.median.Tau < tissueTau, ]
housekeepinggene.median <- log2tpm.median[nullcount.sum == 0 & !is.na(log2tpm.Tau.max) & log2tpm.Tau.max < sampleTau & !is.na(log2tpm.median.Tau) & log2tpm.median.Tau < tissueTau, ]
housekeepinggene.Tau <- log2tpm.Tau[nullcount.sum == 0 & !is.na(log2tpm.Tau.max) & log2tpm.Tau.max < sampleTau & !is.na(log2tpm.median.Tau) & log2tpm.median.Tau < tissueTau, ]
rownames(housekeepinggene.median) <- housekeepinggene[,1]
rownames(housekeepinggene.Tau) <- housekeepinggene[,1]
# housekeepinggene[grep("GAPDH", as.matrix(housekeepinggene[,1])), 1]

## check ENSG00000204983.8|PRSS1 Pancreas
# tissuespecificgene[grep("PRSS1",as.matrix(tissuespecificgene[,1])),1]
# raw.table[grep("PRSS1",as.matrix(raw.table[,1])),1:5][10,]
# ENSG00000153266.8|FEZF2	Brain
# ENSG00000055957.6|ITIH1	Liver
# ENSG00000136574.13|GATA4	Artery - Coronary,Heart - Atrial Appendage,Heart - Left Ventricle,Ovary,Testis
print("Looking for tissue specific genes")
tissuespecificgene <- data.frame()
tissueflags <- c()
rindex <- c()
for(i in 1:nrow(all.stats)){
	if(!is.na(log2tpm.median.Tau[i]) && log2tpm.median.Tau[i] > 1 - tissueTau && !is.na(log2tpm.Tau.min[i]) && log2tpm.Tau.min[i] < sampleTau && min(nullcount[i,], na.rm = T) == 0){
		flag <- 0
		tflags <- c()
		for(j in 1:length(tissues)){
			temp <- log2tpm.median[i,j]/fmax(log2tpm.median[i,])
			if(nullcount[i,j] == 0 && !is.na(log2tpm.Tau[i,j]) && log2tpm.Tau[i,j] < sampleTau &&  !is.na(temp) && temp > 1 - tissueTau){
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
tissuespecificgene.median <- log2tpm.median[rindex, ]
tissuespecificgene.Tau <- log2tpm.Tau[rindex, ]
rownames(tissuespecificgene.median) <- tissuespecificgene[,1]
rownames(tissuespecificgene.Tau) <- tissuespecificgene[,1]

print("Writing HKG/TSG output files")
write.table(housekeepinggene, file = paste(outputpre, "HKG.tsv", sep = "."), 
	sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(tissuespecificgene, tissueflags), file = paste(outputpre, "TSG.tsv", sep = "."), 
	sep = "\t", row.names = FALSE, quote = FALSE)

print("Clustering")
breaklists <- c(seq(0, 1, by = 0.01))
colorn <- length(breaklists)
colors <- colorRampPalette(c("blue", "yellow", "red"))(colorn)

breaklists2 <- c(seq(0, 6, by = 0.05),seq(6.1, 16, by = 0.1))
colorn2 <- length(breaklists2)
colors2 <- colorRampPalette(c("blue", "yellow", "red"))(colorn2)

png(filename = paste(outputpre, "HKG.median.png", sep = "."), width = 1000, height = 1500)
data <- data.frame(housekeepinggene.median, housekeepinggene[,4])
p1 <- pheatmap(data, scale = "none", show_rownames = F, show_colnames = T, 
         color = colors2, cluster_cols = T, breaks = breaklists2,
         clustering_distance_cols = "euclidean", clustering_distance_rows = "euclidean", 
         clustering_method = "ward.D2"
         )
dev.off()
cluster <- cutree(p1$tree_row, k = 5)
write.table(data.frame(cluster[p1$tree_row$order], 
	rownames(data)[p1$tree_row$order],
	data[p1$tree_row$order, p1$tree_col$order]), 
	file = paste(outputpre, "HKG.median.pheatmap.tsv", sep = "."), 
	sep = "\t", row.names = FALSE, quote = FALSE)

png(filename = paste(outputpre, "HKG.Tau.png", sep = "."), width = 1000, height = 1500)
data <- data.frame(housekeepinggene.Tau,housekeepinggene[,5])
pheatmap(data[p1$tree_row$order,p1$tree_col$order], scale = "none", 
	show_rownames = F, show_colnames = T, color = colors, 
	cluster_cols = F, cluster_rows = F, breaks = breaklists)
dev.off()
write.table(data.frame(cluster[p1$tree_row$order], 
	rownames(data)[p1$tree_row$order],
	data[p1$tree_row$order, p1$tree_col$order]), 
	file = paste(outputpre, "HKG.Tau.pheatmap.tsv", sep = "."), 
	sep = "\t", row.names = FALSE, quote = FALSE)

png(filename = paste(outputpre, "TSG.median.png", sep = "."), width = 1000, height = 1500)
data <- data.frame(tissuespecificgene.median, tissuespecificgene[,4])
p2 <- pheatmap(data, scale = "none", show_rownames = F, show_colnames = T, 
         color = colors2, cluster_cols = T, breaks = breaklists2,
         clustering_distance_cols = "euclidean", clustering_distance_rows = "euclidean", 
         clustering_method = "ward.D2"
         )
dev.off()
cluster <- cutree(p2$tree_row, k = 5)
write.table(data.frame(cluster[p2$tree_row$order], 
	tissueflags[p2$tree_row$order],
	rownames(data)[p2$tree_row$order],
	data[p2$tree_row$order, p2$tree_col$order]), 
	file = paste(outputpre, "TSG.median.pheatmap.tsv", sep = "."), 
	sep = "\t", row.names = FALSE, quote = FALSE)

png(filename = paste(outputpre, "TSG.Tau.png", sep = "."), width = 1000, height = 1500)
data <- data.frame(tissuespecificgene.Tau, tissuespecificgene[,5])
pheatmap(data[p2$tree_row$order,p2$tree_col$order], scale = "none", 
	show_rownames = F, show_colnames = T, color = colors, 
	cluster_cols = F, cluster_rows = F, breaks = breaklists)
dev.off()
write.table(data.frame(cluster[p2$tree_row$order], 
	rownames(data)[p2$tree_row$order],
	data[p2$tree_row$order, p2$tree_col$order]), 
	file = paste(outputpre, "TSG.Tau.pheatmap.tsv", sep = "."), 
	sep = "\t", row.names = FALSE, quote = FALSE)

