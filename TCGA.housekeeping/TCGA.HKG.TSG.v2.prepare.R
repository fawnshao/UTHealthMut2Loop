####################
# input data is log2(RSEM + 1)
####################
library(data.table)
# library(pheatmap)
################## some mathematics ##################
# IQR = Q3 - Q1
# outlier < Q1 - 1.5*IQR or outlier > Q3 + 1.5*IQR

# 1 <= i <= n
# z = (xi - mean(x)) / sd(x)
# sd = sqrt(sum((xi - mean(x))^2) / sqrt(n - 1)) (sample sd)
# se = sd / sqrt(n)

# http://www.statisticshowto.com/measures-variation/
# https://people.richland.edu/james/lecture/m170/ch03-var.html
# Different Measures of Variation
# The Range (RANGE = MAXIMUM - MINIMUM)
# Quartiles
# Interquartile Range (IQR)
# Variance (sd)

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5444245/
# Tau = sum(1 - xi/max(x)) / (n - 1)
# https://www.ncbi.nlm.nih.gov/pubmed/15388519
# genes with midrange profiles having 0.15<tau<0.85 were found to constitute >50% of all expression patterns

# https://www.qiagen.com/us/resources/faq?id=06a192c2-e72d-42e8-9b40-3171e1eb4cb8&lang=en
# How much RNA does a typical mammalian cell contain?
# FAQ ID -2946
# The RNA content and RNA make up of a cell depend very much on its developmental stage and the type of cell. To estimate the approximate yield of RNA that can be expected from your starting material, we usually calculate that a typical mammalian cell contains 10–30 pg total RNA.

# The majority of RNA molecules are tRNAs and rRNAs. mRNA accounts for only 1–5% of the total cellular RNA although the actual amount depends on the cell type and physiological state. Approximately 360,000 mRNA molecules are present in a single mammalian cell, made up of approximately 12,000 different transcripts with a typical length of around 2 kb. Some mRNAs comprise 3% of the mRNA pool whereas others account for less than 0.1%. These rare or low-abundance mRNAs may have a copy number of only 5–15 molecules per cell.

######################################################
# args <- c("GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct", "GTEx_sample.tissue.txt", "v2", "0.1")
args <- c("tcga_RSEM_gene_tpm", "tcga_RSEM_gene_tpm.samples.sim", "v2.T", "0.1")
# args <- c("tcga_RSEM_gene_tpm", "tcga_RSEM_gene_tpm.samples.sim", "v2.N", "0.1")
# expression tissuename outputpre tissueTau sampleTau expressioncutoff[log2(tpm+1)]

### some functions ###
###+++###
fmedian <- function(x) {
	if(!all(is.na(x))) {
		res <- median(x, na.rm = TRUE)
	 	} 
	else {
		res <- NA
	 	}
	return(res)
}
###***###***###	

###+++###
fiqr <- function(x) {
	if(!all(is.na(x))) {
		res <- quantile(x, probs = 0.75, na.rm = T) - quantile(x, probs = 0.25, na.rm = T)
	 	} 
	else {
		res <- NA
	 	}
	return(res)
}
###***###***###	

###+++###
fsd <- function(x) {
	if(!all(is.na(x))) {
		res <- sd(x, na.rm = T)
	 	} 
	else {
		res <- NA
	 	}
	return(res)
}
###***###***###	

###+++###	
#Function require a vector with expression of one gene in different tissues/samples.
#Mean is calculated taking in account tissues with 0 expression. 2+0+4=2
fmean <- function(x){
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
foutlier <- function(x) {
	if(!all(is.na(x))) {
		z <- (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = T)
		res <- length(z[abs(z) > 2])
	 	} 
	else {
		res <- NA
	 	}
	return(res)
}
###***###***###	

###+++###	
#Function require a vector with expression of one gene in different tissues/samples.
#Max is calculated taking in account tissues with 0 expression. 2+0+4=2
fmax <- function(x){
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
####### GETx
# fTau <- function(x){
# 	if(all(!is.na(x)))
#  	{
#  		if(min(x, na.rm = TRUE) >= 0)
# 		{
#  			mx <- fmax(x)
#  			if(mx != 0 && !is.na(mx))
#  			{
#  				x <- (1-(x/mx))
#  				res <- sum(x, na.rm = TRUE)
#  				res <- res/(length(x)-1)
#  			} else {
#  				res <- 0
#  			}
#  		} else {
#  		res <- NA
#  		#print("Expression values have to be positive!")
#  		} 
#  	} else {
#  		res <- NA
#  		#print("No data for this gene avalable.")
#  	} 
#  	return(res)
# }
####### TCGA
fTau <- function(x)
{
	if(all(!is.na(x)))
 	{
 		# if(min(x, na.rm = TRUE) >= 0)
		# {
 			mx <- fmax(x)
 			if(mx != 0 && !is.na(mx))
 			{
 				x <- (1-(x/mx))
 				res <- sum(x, na.rm = TRUE)
 				res <- res/(length(x)-1)
 			} else {
 				res <- 0
 			}
 		# } else {
 		# res <- NA
 		# #print("Expression values have to be positive!")
 		# } 
 	} else {
 		res <- NA
 		#print("No data for this gene avalable.")
 	} 
 	return(res)
}
###***###***###

## read in data
print("Reading Data")
raw.table <- fread(args[1], sep = "\t", header = T)
# Read 56202 rows and 11689 (of 11689) columns from 2.622 GB file in 00:06:16
tpm <- as.matrix(raw.table[,-1])
rownames(tpm) <- as.matrix(raw.table[,1])
info <- as.matrix(read.table(args[2], sep = "\t", header = T))
samplecount <- table(info[,2])
tissues <- names(samplecount[samplecount > 20])
tissues <- tissues[-grep("Normal", tissues)]
# tissues <- tissues[grep("Normal", tissues)]

outputpre <- args[3]
nullexpression <- as.numeric(args[4])
## calculate
print("Calculating nullcount, median/mean and Tau for each tissue")
nullcount <- matrix(ncol = length(tissues), nrow = nrow(tpm))
log2tpm.median <- matrix(ncol = length(tissues), nrow = nrow(tpm))
log2tpm.mean <- matrix(ncol = length(tissues), nrow = nrow(tpm))
log2tpm.sd <- matrix(ncol = length(tissues), nrow = nrow(tpm))
log2tpm.iqr <- matrix(ncol = length(tissues), nrow = nrow(tpm))
log2tpm.outlier <- matrix(ncol = length(tissues), nrow = nrow(tpm))
log2tpm.Tau <- matrix(ncol = length(tissues), nrow = nrow(tpm))
for(i in 1:length(tissues)){
	print(tissues[i])
	# z <- log2(tpm[, info[,2] == tissues[i]] + 1)
	z <- tpm[, info[,2] == tissues[i]]
	nullcount[,i] <- apply(z, 1, function(x) {length(x[x <= nullexpression])})
	log2tpm.median[,i] <- apply(z, 1, fmedian)
	log2tpm.mean[,i] <- apply(z, 1, fmean)
	log2tpm.sd[,i] <- apply(z, 1, fsd)
	log2tpm.iqr[,i] <- apply(z, 1, fiqr)
	log2tpm.outlier[,i] <- apply(z, 1, foutlier)
	log2tpm.Tau[,i] <- apply(z, 1, fTau)
}
colnames(nullcount) <- tissues
colnames(log2tpm.median) <- tissues
colnames(log2tpm.mean) <- tissues
colnames(log2tpm.sd) <- tissues
colnames(log2tpm.iqr) <- tissues
colnames(log2tpm.outlier) <- tissues
colnames(log2tpm.Tau) <- tissues

print("Calculating nullcount, mean and Tau across tissues")
log2tpm.median.mean <- apply(log2tpm.median, 1, fmean)
log2tpm.median.Tau <- apply(log2tpm.median, 1, fTau)

# nullcount.sum <- apply(nullcount, 1, sum)
# log2tpm.Tau.max <- apply(log2tpm.Tau, 1, function(x) {max(x[is.finite(x)], na.rm = T)})
# log2tpm.Tau.min <- apply(log2tpm.Tau, 1, function(x) {min(x[is.finite(x)], na.rm = T)})

### it takes ~30min until this step

# print("Saving RData")
# save.image("v1.3.Tau.RData")
# ### it takes ~50min until this step

print("Writing raw output files")
# write.table(data.frame(rownames(tpm), nullcount), file = paste(outputpre, args[4], "nullcount.tsv", sep = "."), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(rownames(tpm), nullcount), file = paste(outputpre, "nullcount.tsv", sep = "."), 
	sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(rownames(tpm), log2tpm.mean), file = paste(outputpre, "log2tpm.mean.tsv", sep = "."), 
	sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(rownames(tpm), log2tpm.median), file = paste(outputpre, "log2tpm.median.tsv", sep = "."), 
	sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(rownames(tpm), log2tpm.Tau), file = paste(outputpre, "log2tpm.Tau.tsv", sep = "."), 
	sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(rownames(tpm), log2tpm.sd), file = paste(outputpre, "log2tpm.sd.tsv", sep = "."), 
	sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(rownames(tpm), log2tpm.iqr), file = paste(outputpre, "log2tpm.iqr.tsv", sep = "."), 
	sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(rownames(tpm), log2tpm.outlier), file = paste(outputpre, "log2tpm.outlier.tsv", sep = "."), 
	sep = "\t", row.names = FALSE, quote = FALSE)

all.stats <- data.frame(rownames(tpm), log2tpm.median.mean, log2tpm.median.Tau)
write.table(all.stats, file = paste(outputpre, "allstats.tsv", sep = "."), 
	sep = "\t", row.names = FALSE, quote = FALSE)
