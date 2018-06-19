library(data.table)
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
fmin <- function(x){
	if(!all(is.na(x)))
 	{
 		res <- min(x, na.rm = TRUE)
 	} else {
 		res <- NA
 	}
 	return(res)
}
###+++###

print("Reading Data")
# outputpre tissuename tissueTau sampleTau
args <- c("v2", "GTEx_sample.tissue.txt", "0.15", "1", "3", "0.1")
outputpre <- args[1]
tau.threshold <- as.numeric(args[3])
sd.threshold <- as.numeric(args[4])
median.threshold <- as.numeric(args[5])
outlier.threshold <- as.numeric(args[6])

all.stats.file <- paste(outputpre, "allstats.tsv", sep = ".")
all.stats <- fread(all.stats.file, sep = "\t", header = T)
nullcount.file <- paste(outputpre, "nullcount.tsv", sep = ".")
log2tpm.median.file <- paste(outputpre, "log2tpm.median.tsv", sep = ".")
log2tpm.sd.file <- paste(outputpre, "log2tpm.sd.tsv", sep = ".")
log2tpm.outlier.file <- paste(outputpre, "log2tpm.outlier.tsv", sep = ".")
nullcount <- as.matrix(fread(nullcount.file, sep = "\t", header = T)[,-1])
log2tpm.median <- as.matrix(fread(log2tpm.median.file, sep = "\t", header = T)[,-1])
log2tpm.sd <- as.matrix(fread(log2tpm.sd.file, sep = "\t", header = T)[,-1])
log2tpm.outlier <- as.matrix(fread(log2tpm.outlier.file, sep = "\t", header = T)[,-1])

info <- as.matrix(read.table(args[2], sep = "\t"))
tissues <- unique(info[,2])
tissues.count <- table(info[,2])
# tissues <- tissues[-c(7,24,25,53)]
tissues <- tissues[-c(7,24,25,31,53)]
tissues.count <- tissues.count[-c(7,24,25,31,53)]

genes <- as.matrix(all.stats[,1])
log2tpm.median.mean <- as.matrix(all.stats[,2])
log2tpm.median.Tau <- as.matrix(all.stats[,3])
rownames(log2tpm.median) <- genes
colnames(log2tpm.median) <- tissues
rownames(nullcount) <- genes
colnames(nullcount) <- tissues
rownames(log2tpm.sd) <- genes
colnames(log2tpm.sd) <- tissues
rownames(log2tpm.outlier) <- genes
colnames(log2tpm.outlier) <- tissues

# all.stats[grep("ISL1", genes),]

print("Searching HKG")
log2tpm.median.min <- apply(log2tpm.median, 1, fmin)
log2tpm.sd.max <- apply(log2tpm.sd, 1, fmax)
outlier.percentage <- t(apply(log2tpm.outlier, 1, function(x){x/tissues.count}))
log2tpm.outlier.max <- apply(outlier.percentage, 1, fmax)
nullcount.sum <- apply(nullcount, 1, sum)
housekeepinggene <- all.stats[nullcount.sum == 0 & !is.na(log2tpm.sd.max) & log2tpm.sd.max < sd.threshold & !is.na(log2tpm.outlier.max) & log2tpm.outlier.max < outlier.threshold & !is.na(log2tpm.median.min) & log2tpm.median.min > median.threshold & !is.na(log2tpm.median.Tau) & log2tpm.median.Tau < tau.threshold, ]
housekeepinggene.median <- log2tpm.median[nullcount.sum == 0 & !is.na(log2tpm.sd.max) & log2tpm.sd.max < sd.threshold & !is.na(log2tpm.outlier.max) & log2tpm.outlier.max < outlier.threshold & !is.na(log2tpm.median.min) & log2tpm.median.min > median.threshold & !is.na(log2tpm.median.Tau) & log2tpm.median.Tau < tau.threshold, ]
housekeepinggene.sd <- log2tpm.sd[nullcount.sum == 0 & !is.na(log2tpm.sd.max) & log2tpm.sd.max < sd.threshold & !is.na(log2tpm.outlier.max) & log2tpm.outlier.max < outlier.threshold & !is.na(log2tpm.median.min) & log2tpm.median.min > median.threshold & !is.na(log2tpm.median.Tau) & log2tpm.median.Tau < tau.threshold, ]
rownames(housekeepinggene.median) <- as.matrix(housekeepinggene[,1])
rownames(housekeepinggene.sd) <- as.matrix(housekeepinggene[,1])

print("Searching TSG")
# all.stats[grep("ENSG00000168757.8|TSPY2 ", genes),]
log2tpm.median.max <- apply(log2tpm.median, 1, fmax)
log2tpm.sd.min <- apply(log2tpm.sd, 1, fmin)
log2tpm.outlier.min <- apply(outlier.percentage, 1, fmin)
nullcount.min <- apply(nullcount, 1, fmin)

tissuespecificgene <- data.frame()
tissueflags <- c()
rindex <- c()
for(i in 1:nrow(all.stats)){
	if(!is.na(nullcount.min[i]) && nullcount.min[i] == 0 && !is.na(log2tpm.sd.max[i]) && log2tpm.sd.min[i] < sd.threshold && !is.na(log2tpm.outlier.min[i]) && log2tpm.outlier.min[i] < outlier.threshold && !is.na(log2tpm.median.Tau[i]) && log2tpm.median.Tau[i] > 1 - tau.threshold && !is.na(log2tpm.median.max[i]) && log2tpm.median.max[i] > median.threshold){
		# rindex <- c(rindex, i)
		flag <- 0
		tflags <- c()
		for(j in 1:length(tissues)){
			temp <- log2tpm.median[i,j] / fmax(log2tpm.median[i,])
			if(nullcount[i,j] == 0 && !is.na(log2tpm.sd[i,j]) && log2tpm.sd[i,j] < sd.threshold && !is.na(outlier.percentage[i,j]) && outlier.percentage[i,j] < outlier.threshold && !is.na(temp) && temp > 1 - tau.threshold && !is.na(log2tpm.median[i,j]) && log2tpm.median[i,j] > median.threshold){
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
tissuespecificgene.sd <- log2tpm.sd[rindex, ]
rownames(tissuespecificgene.median) <- as.matrix(tissuespecificgene[,1])
rownames(tissuespecificgene.sd) <- as.matrix(tissuespecificgene[,1])


print("Writing output files")
write.table(data.frame(housekeepinggene, housekeepinggene.median, housekeepinggene.sd), 
	file = paste(outputpre, "HKG.tsv", sep = "."), 
	sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(tissuespecificgene, tissuespecificgene.median, tissuespecificgene.sd,tissueflags), file = paste(outputpre, "TSG.tsv", sep = "."), 
	sep = "\t", row.names = FALSE, quote = FALSE)
