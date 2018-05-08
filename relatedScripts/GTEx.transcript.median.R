library(data.table)
library(pheatmap)

args <- c("GTEx.transcript.tpm", "GTEx.transcript.id", "GTEx.transcript.tissue.txt", "GTEx.transcript")
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
# tpm1 <- as.matrix(raw.table[1:50000,-1])
# tpm2 <- as.matrix(raw.table[50001:70000,-1])
# tpm3 <- as.matrix(raw.table[70001:90000,-1])
# tpm4 <- as.matrix(raw.table[90001:110000,-1])
# tpm5 <- as.matrix(raw.table[110001:130000,-1])
# tpm6 <- as.matrix(raw.table[130001:150000,-1])
# tpm7 <- as.matrix(raw.table[150001:196520,-1])
# tpm <- rbind.data.frame(tpm1, tpm2, tpm3, tpm4, tpm5, tpm6, tpm7)
tpm <- data.frame()
for(i in 1:19){
	j <- (i - 1) * 10000 + 1
	k <- i * 10000
	tpm <- rbind.data.frame(tpm, as.matrix(raw.table[j:k,]))
}
tpm <- rbind.data.frame(tpm, as.matrix(raw.table[190001:196520,]))
rownames(tpm) <- as.matrix(read.table(args[2], header = T))
info <- as.matrix(read.table(args[3], sep = "\t"))
tissues <- unique(info[,2])

outputpre <- args[4]

## calculate
print("Calculating median/mean for each tissue")
log2tpm.median <- matrix(ncol = length(tissues), nrow = nrow(tpm))
log2tpm.mean <- matrix(ncol = length(tissues), nrow = nrow(tpm))
for(i in 1:length(tissues)){
	z <- log2(tpm[, info[,2] == tissues[i]] + 1)
	log2tpm.median[,i] <- apply(z, 1, fmedian)
	log2tpm.mean[,i] <- apply(z, 1, fmean)
}
colnames(log2tpm.median) <- tissues
colnames(log2tpm.mean) <- tissues

print("Writing raw output files")
write.table(data.frame(rownames(tpm), log2tpm.mean), file = paste(outputpre, "log2tpm.mean.tsv", sep = "."), 
	sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(rownames(tpm), log2tpm.median), file = paste(outputpre, "log2tpm.median.tsv", sep = "."), 
	sep = "\t", row.names = FALSE, quote = FALSE)


