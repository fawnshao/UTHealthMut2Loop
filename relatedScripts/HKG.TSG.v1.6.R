library(data.table)
library(pheatmap)

args <- c("GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct", 
	"GTEx_sample.tissue.txt", "v1.6", "0.25", "0.5", "0.1", "10")
# expression tissuename outputpre tissueTau sampleTau expressioncutoff[raw(tpm+1)]
# use the raw TPM value

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
rawtpm.median <- matrix(ncol = length(tissues), nrow = nrow(tpm))
rawtpm.mean <- matrix(ncol = length(tissues), nrow = nrow(tpm))
rawtpm.Tau <- matrix(ncol = length(tissues), nrow = nrow(tpm))
for(i in 1:length(tissues)){
	z <- tpm[, info[,2] == tissues[i]]
	nullcount[,i] <- apply(z, 1, function(x) {length(x[x < nullexpression])})
	rawtpm.median[,i] <- apply(z, 1, fmedian)
	rawtpm.mean[,i] <- apply(z, 1, fmean)
	rawtpm.Tau[,i] <- apply(z, 1, fTau)
}
colnames(nullcount) <- tissues
colnames(rawtpm.median) <- tissues
colnames(rawtpm.mean) <- tissues
colnames(rawtpm.Tau) <- tissues

print("Calculating nullcount, mean and Tau across tissues")
rawtpm.median.mean <- apply(rawtpm.median, 1, fmean)
rawtpm.median.Tau <- apply(rawtpm.median, 1, fTau)
rawtpm.mean.mean <- apply(rawtpm.mean, 1, fmean)
rawtpm.mean.Tau <- apply(rawtpm.mean, 1, fTau)

nullcount.sum <- apply(nullcount, 1, sum)
rawtpm.Tau.max <- apply(rawtpm.Tau, 1, function(x) {max(x[is.finite(x)], na.rm = T)})
rawtpm.Tau.min <- apply(rawtpm.Tau, 1, function(x) {min(x[is.finite(x)], na.rm = T)})

### it takes ~30min until this step

# print("Saving RData")
# save.image("v1.3.Tau.RData")
# ### it takes ~50min until this step

print("Writing raw output files")
all.stats <- data.frame(rownames(tpm), rawtpm.mean.mean, rawtpm.mean.Tau, 
	rawtpm.median.mean, rawtpm.median.Tau, 
	nullcount.sum, rawtpm.Tau.max, rawtpm.Tau.min, 
	rawtpm.mean, rawtpm.median, rawtpm.Tau)
rownames(all.stats) <- rownames(tpm)
write.table(data.frame(rownames(tpm), nullcount), file = paste(outputpre, "nullcount.tsv", sep = "."), 
	sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(rownames(tpm), rawtpm.mean), file = paste(outputpre, "rawtpm.mean.tsv", sep = "."), 
	sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(rownames(tpm), rawtpm.median), file = paste(outputpre, "rawtpm.median.tsv", sep = "."), 
	sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(rownames(tpm), rawtpm.Tau), file = paste(outputpre, "rawtpm.Tau.tsv", sep = "."), 
	sep = "\t", row.names = FALSE, quote = FALSE)
write.table(all.stats, file = paste(outputpre, "allstats.tsv", sep = "."), 
	sep = "\t", row.names = FALSE, quote = FALSE)

#load("v1.3.Tau.RData")
#library(pheatmap)
print("Looking for housekeeping genes")
rawtpm.median.min <- apply(rawtpm.median, 1, function(x){min(x, na.rm = T)})
housekeepinggene <- all.stats[nullcount.sum == 0 & !is.na(rawtpm.Tau.max) & rawtpm.Tau.max < sampleTau & !is.na(rawtpm.median.Tau) & rawtpm.median.Tau < tissueTau & rawtpm.median.min > medianthreshold, ]
housekeepinggene.median <- rawtpm.median[nullcount.sum == 0 & !is.na(rawtpm.Tau.max) & rawtpm.Tau.max < sampleTau & !is.na(rawtpm.median.Tau) & rawtpm.median.Tau < tissueTau & rawtpm.median.min > medianthreshold, ]
housekeepinggene.Tau <- rawtpm.Tau[nullcount.sum == 0 & !is.na(rawtpm.Tau.max) & rawtpm.Tau.max < sampleTau & !is.na(rawtpm.median.Tau) & rawtpm.median.Tau < tissueTau & rawtpm.median.min > medianthreshold, ]
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
	if(!is.na(rawtpm.median.Tau[i]) && rawtpm.median.Tau[i] > 1 - tissueTau && !is.na(rawtpm.Tau.min[i]) && rawtpm.Tau.min[i] < sampleTau && min(nullcount[i,], na.rm = T) == 0){
		flag <- 0
		tflags <- c()
		for(j in 1:length(tissues)){
			temp <- rawtpm.median[i,j]/fmax(rawtpm.median[i,])
			if(nullcount[i,j] == 0 && !is.na(rawtpm.Tau[i,j]) && rawtpm.Tau[i,j] < sampleTau &&  !is.na(temp) && temp > 1 - tissueTau && rawtpm.median[i,j] > medianthreshold){
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
tissuespecificgene.median <- rawtpm.median[rindex, ]
tissuespecificgene.Tau <- rawtpm.Tau[rindex, ]
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

