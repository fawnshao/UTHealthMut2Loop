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
args <- c("v2.1", "GTEx_sample.tissue.age.txt", "0.15", "1.5", "2", "0.25") #, "1")
# args <- commandArgs(TRUE)
# outputpre <- args[1]
tau.threshold <- as.numeric(args[3])
sd.threshold <- as.numeric(args[4])
median.threshold <- as.numeric(args[5])
tau.threshold2 <- as.numeric(args[6])
# outlier.threshold <- as.numeric(args[6])

info <- as.matrix(read.table(args[2], sep = "\t", header = T))
tissues <- unique(info[,2])
tissues.count <- table(info[,2])
# tissues <- tissues[-c(7,24,25,53)]
tissues <- tissues[-c(7,24,25,31,53)]
tissues.count <- tissues.count[-c(7,24,25,31,53)]
tissues <- tissues[-c(7,24,25,31,53)]
AgeGroup <- info[,5]
AgeGroup[info[,5] == "20-29"] <- "20-39"
AgeGroup[info[,5] == "30-39"] <- "20-39"
AgeGroup[info[,5] == "60-69"] <- "60-79"
AgeGroup[info[,5] == "70-79"] <- "60-79"
info <- data.frame(info, AgeGroup = AgeGroup)
ages <- unique(info[,7])
# ages.count <- table(info[,7])

for(k in 1:length(ages)){
	print(ages[k])
	outputpre <- paste(args[1], ages[k], sep = ".")
	all.stats.file <- paste(outputpre, "allstats.tsv", sep = ".")
	all.stats <- fread(all.stats.file, sep = "\t", header = T)
	nullcount.file <- paste(outputpre, "nullcount.tsv", sep = ".")
	# nullcount.file <- paste(outputpre, "3.nullcount.tsv", sep = ".")
	log2tpm.median.file <- paste(outputpre, "log2tpm.median.tsv", sep = ".")
	log2tpm.sd.file <- paste(outputpre, "log2tpm.sd.tsv", sep = ".")
	log2tpm.outlier.file <- paste(outputpre, "log2tpm.outlier.tsv", sep = ".")
	nullcount <- as.matrix(fread(nullcount.file, sep = "\t", header = T)[,-1])
	log2tpm.median <- as.matrix(fread(log2tpm.median.file, sep = "\t", header = T)[,-1])
	log2tpm.sd <- as.matrix(fread(log2tpm.sd.file, sep = "\t", header = T)[,-1])
	log2tpm.outlier <- as.matrix(fread(log2tpm.outlier.file, sep = "\t", header = T)[,-1])

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
	# outlier.percentage <- t(apply(log2tpm.outlier, 1, function(x){x/tissues.count}))
	# log2tpm.outlier.max <- apply(outlier.percentage, 1, fmax)
	nullcount.sum <- apply(nullcount, 1, sum)
	# housekeepinggene <- all.stats[nullcount.sum == 0 & !is.na(log2tpm.sd.max) & log2tpm.sd.max < sd.threshold & !is.na(log2tpm.outlier.max) & log2tpm.outlier.max < outlier.threshold & !is.na(log2tpm.median.min) & log2tpm.median.min > median.threshold & !is.na(log2tpm.median.Tau) & log2tpm.median.Tau < tau.threshold, ]
	# housekeepinggene.median <- log2tpm.median[nullcount.sum == 0 & !is.na(log2tpm.sd.max) & log2tpm.sd.max < sd.threshold & !is.na(log2tpm.outlier.max) & log2tpm.outlier.max < outlier.threshold & !is.na(log2tpm.median.min) & log2tpm.median.min > median.threshold & !is.na(log2tpm.median.Tau) & log2tpm.median.Tau < tau.threshold, ]
	# housekeepinggene.sd <- log2tpm.sd[nullcount.sum == 0 & !is.na(log2tpm.sd.max) & log2tpm.sd.max < sd.threshold & !is.na(log2tpm.outlier.max) & log2tpm.outlier.max < outlier.threshold & !is.na(log2tpm.median.min) & log2tpm.median.min > median.threshold & !is.na(log2tpm.median.Tau) & log2tpm.median.Tau < tau.threshold, ]
	housekeepinggene <- all.stats[nullcount.sum == 0 & !is.na(log2tpm.sd.max) & log2tpm.sd.max < sd.threshold & !is.na(log2tpm.median.min) & log2tpm.median.min > median.threshold & !is.na(log2tpm.median.Tau) & log2tpm.median.Tau < tau.threshold, ]
	housekeepinggene.median <- log2tpm.median[nullcount.sum == 0 & !is.na(log2tpm.sd.max) & log2tpm.sd.max < sd.threshold & !is.na(log2tpm.median.min) & log2tpm.median.min > median.threshold & !is.na(log2tpm.median.Tau) & log2tpm.median.Tau < tau.threshold, ]
	housekeepinggene.sd <- log2tpm.sd[nullcount.sum == 0 & !is.na(log2tpm.sd.max) & log2tpm.sd.max < sd.threshold & !is.na(log2tpm.median.min) & log2tpm.median.min > median.threshold & !is.na(log2tpm.median.Tau) & log2tpm.median.Tau < tau.threshold, ]
	rownames(housekeepinggene.median) <- as.matrix(housekeepinggene[,1])
	rownames(housekeepinggene.sd) <- as.matrix(housekeepinggene[,1])

	housekeepinggene2 <- all.stats[nullcount.sum == 0 & !is.na(log2tpm.sd.max) & log2tpm.sd.max < sd.threshold & !is.na(log2tpm.median.min) & log2tpm.median.min > median.threshold & !is.na(log2tpm.median.Tau) & log2tpm.median.Tau >= tau.threshold & log2tpm.median.Tau < tau.threshold2, ]
	housekeepinggene2.median <- log2tpm.median[nullcount.sum == 0 & !is.na(log2tpm.sd.max) & log2tpm.sd.max < sd.threshold & !is.na(log2tpm.median.min) & log2tpm.median.min > median.threshold & !is.na(log2tpm.median.Tau) & log2tpm.median.Tau >= tau.threshold & log2tpm.median.Tau < tau.threshold2, ]
	housekeepinggene2.sd <- log2tpm.sd[nullcount.sum == 0 & !is.na(log2tpm.sd.max) & log2tpm.sd.max < sd.threshold & !is.na(log2tpm.median.min) & log2tpm.median.min > median.threshold & !is.na(log2tpm.median.Tau) & log2tpm.median.Tau >= tau.threshold & log2tpm.median.Tau < tau.threshold2, ]
	rownames(housekeepinggene2.median) <- as.matrix(housekeepinggene2[,1])
	rownames(housekeepinggene2.sd) <- as.matrix(housekeepinggene2[,1])

	print("Searching TSG")
	# all.stats[grep("ENSG00000168757.8|TSPY2 ", genes),]
	log2tpm.median.max <- apply(log2tpm.median, 1, fmax)
	log2tpm.sd.min <- apply(log2tpm.sd, 1, fmin)
	# log2tpm.outlier.min <- apply(outlier.percentage, 1, fmin)
	nullcount.min <- apply(nullcount, 1, fmin)

	tissuespecificgene <- data.frame()
	tissueflags <- c()
	rindex <- c()
	for(i in 1:nrow(all.stats)){
		# if(!is.na(nullcount.min[i]) && nullcount.min[i] == 0 && !is.na(log2tpm.sd.max[i]) && log2tpm.sd.min[i] < sd.threshold && !is.na(log2tpm.outlier.min[i]) && log2tpm.outlier.min[i] < outlier.threshold && !is.na(log2tpm.median.Tau[i]) && log2tpm.median.Tau[i] > 1 - tau.threshold && !is.na(log2tpm.median.max[i]) && log2tpm.median.max[i] > median.threshold){
		if(!is.na(nullcount.min[i]) && nullcount.min[i] == 0 && !is.na(log2tpm.sd.max[i]) && log2tpm.sd.min[i] < sd.threshold && !is.na(log2tpm.median.Tau[i]) && log2tpm.median.Tau[i] > 1 - tau.threshold && !is.na(log2tpm.median.max[i]) && log2tpm.median.max[i] > median.threshold){
			# rindex <- c(rindex, i)
			flag <- 0
			tflags <- c()
			for(j in 1:length(tissues)){
				temp <- log2tpm.median[i,j] / fmax(log2tpm.median[i,])
				# if(nullcount[i,j] == 0 && !is.na(log2tpm.sd[i,j]) && log2tpm.sd[i,j] < sd.threshold && !is.na(outlier.percentage[i,j]) && outlier.percentage[i,j] < outlier.threshold && !is.na(temp) && temp > 1 - tau.threshold && !is.na(log2tpm.median[i,j]) && log2tpm.median[i,j] > median.threshold){
				if(nullcount[i,j] == 0 && !is.na(log2tpm.sd[i,j]) && log2tpm.sd[i,j] < sd.threshold && !is.na(temp) && temp > 1 - tau.threshold && !is.na(log2tpm.median[i,j]) && log2tpm.median[i,j] > median.threshold){
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

	tissuespecificgene2 <- data.frame()
	tissueflags2 <- c()
	rindex2 <- c()
	for(i in 1:nrow(all.stats)){
		if(!is.na(nullcount.min[i]) && nullcount.min[i] == 0 && !is.na(log2tpm.sd.max[i]) && log2tpm.sd.min[i] < sd.threshold && !is.na(log2tpm.median.Tau[i]) && log2tpm.median.Tau[i] <= 1 - tau.threshold && log2tpm.median.Tau[i] > 1 - tau.threshold2 && !is.na(log2tpm.median.max[i]) && log2tpm.median.max[i] > median.threshold){
			# rindex <- c(rindex, i)
			flag <- 0
			tflags <- c()
			for(j in 1:length(tissues)){
				temp <- log2tpm.median[i,j] / fmax(log2tpm.median[i,])
				if(nullcount[i,j] == 0 && !is.na(log2tpm.sd[i,j]) && log2tpm.sd[i,j] < sd.threshold && !is.na(temp) && temp <= 1 - tau.threshold && temp > 1 - tau.threshold2 && !is.na(log2tpm.median[i,j]) && log2tpm.median[i,j] > median.threshold){
					flag <- 1
					tflags <- c(tflags, tissues[j])
				}
			}
			if(flag == 1){
				tissuespecificgene2 <- rbind.data.frame(tissuespecificgene2, all.stats[i,])
				tissueflags2 <- c(tissueflags2, paste(tflags, collapse = ","))
				rindex2 <- c(rindex2, i)
			}
		}
	}
	tissuespecificgene2.median <- log2tpm.median[rindex2, ]
	tissuespecificgene2.sd <- log2tpm.sd[rindex2, ]
	rownames(tissuespecificgene2.median) <- as.matrix(tissuespecificgene2[,1])
	rownames(tissuespecificgene2.sd) <- as.matrix(tissuespecificgene2[,1])


	print("Writing output files")
	write.table(data.frame(housekeepinggene, housekeepinggene.median, housekeepinggene.sd), 
		file = paste(outputpre, "HKG.tsv", sep = "."), 
		sep = "\t", row.names = FALSE, quote = FALSE)
	write.table(data.frame(tissuespecificgene, tissuespecificgene.median, tissuespecificgene.sd,tissueflags), file = paste(outputpre, "TSG.tsv", sep = "."), 
		sep = "\t", row.names = FALSE, quote = FALSE)
	write.table(data.frame(housekeepinggene2, housekeepinggene2.median, housekeepinggene2.sd), 
		file = paste(outputpre, "HKG.2.tsv", sep = "."), 
		sep = "\t", row.names = FALSE, quote = FALSE)
	write.table(data.frame(tissuespecificgene2, tissuespecificgene2.median, tissuespecificgene2.sd,tissueflags2), file = paste(outputpre, "TSG.2.tsv", sep = "."), 
		sep = "\t", row.names = FALSE, quote = FALSE)
}
