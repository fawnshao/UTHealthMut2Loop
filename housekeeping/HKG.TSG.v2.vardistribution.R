library(data.table)
library(pheatmap)
library(ggplot2)
library(reshape2)
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

print("Reading Data")
# outputpre tissuename tissueTau sampleTau
args <- c("v2", "GTEx_sample.tissue.txt")
outputpre <- args[1]
all.stats.file <- paste(outputpre, "allstats.tsv", sep = ".")
all.stats <- fread(all.stats.file, sep = "\t", header = T)
nullcount.file <- paste(outputpre, "nullcount.tsv", sep = ".")
log2tpm.mean.file <- paste(outputpre, "log2tpm.mean.tsv", sep = ".")
log2tpm.median.file <- paste(outputpre, "log2tpm.median.tsv", sep = ".")
log2tpm.Tau.file <- paste(outputpre, "log2tpm.Tau.tsv", sep = ".") 
log2tpm.sd.file <- paste(outputpre, "log2tpm.sd.tsv", sep = ".")
log2tpm.iqr.file <- paste(outputpre, "log2tpm.iqr.tsv", sep = ".")
log2tpm.outlier.file <- paste(outputpre, "log2tpm.outlier.tsv", sep = ".")
nullcount <- as.matrix(fread(nullcount.file, sep = "\t", header = T)[,-1])
log2tpm.mean <- as.matrix(fread(log2tpm.mean.file, sep = "\t", header = T)[,-1])
log2tpm.median <- as.matrix(fread(log2tpm.median.file, sep = "\t", header = T)[,-1])
log2tpm.Tau <- as.matrix(fread(log2tpm.Tau.file, sep = "\t", header = T)[,-1])
log2tpm.sd <- as.matrix(fread(log2tpm.sd.file, sep = "\t", header = T)[,-1])
log2tpm.iqr <- as.matrix(fread(log2tpm.iqr.file, sep = "\t", header = T)[,-1])
log2tpm.outlier <- as.matrix(fread(log2tpm.outlier.file, sep = "\t", header = T)[,-1])

info <- as.matrix(read.table(args[2], sep = "\t"))
tissues <- unique(info[,2])
tissues <- tissues[-c(7,24,25,53)]

genes <- as.matrix(all.stats[,1])
log2tpm.median.mean <- as.matrix(all.stats[,2])
log2tpm.median.Tau <- as.matrix(all.stats[,3])
rownames(log2tpm.mean) <- genes
colnames(log2tpm.mean) <- tissues
rownames(log2tpm.median) <- genes
colnames(log2tpm.median) <- tissues
rownames(log2tpm.Tau) <- genes
colnames(log2tpm.Tau) <- tissues
rownames(nullcount) <- genes
colnames(nullcount) <- tissues
rownames(log2tpm.sd) <- genes
colnames(log2tpm.sd) <- tissues
rownames(log2tpm.iqr) <- genes
colnames(log2tpm.iqr) <- tissues
rownames(log2tpm.outlier) <- genes
colnames(log2tpm.outlier) <- tissues

# all.stats[grep("ISL1", genes),]

print("Plotting distribution for mean and Tau")

data <- melt(log2tpm.Tau)
colnames(data) <- c("gene", "tissue", "Tau")
pdf(file = paste(args[1], "Tau.pdf", sep = "."), width = 20, height = 10)
ggplot(data, aes(x = tissue, y = Tau, fill = tissue)) + 
	geom_violin() + 
	theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1))
dev.off()

data <- melt(log2tpm.outlier)
colnames(data) <- c("gene", "tissue", "outlier")
pdf(file = paste(args[1], "outliercounts.pdf", sep = "."), width = 20, height = 10)
ggplot(data, aes(x = outlier, y = Tau, fill = tissue)) + 
	geom_violin() + 
	theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1))
dev.off()

data <- melt(log2tpm.sd)
colnames(data) <- c("gene", "tissue", "sd")
pdf(file = paste(args[1], "sd.pdf", sep = "."), width = 20, height = 10)
ggplot(data, aes(x = sd, y = Tau, fill = tissue)) + 
	geom_violin() + 
	theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1))
dev.off()

data <- melt(log2tpm.iqr)
colnames(data) <- c("gene", "tissue", "iqr")
pdf(file = paste(args[1], "iqr.pdf", sep = "."), width = 20, height = 10)
ggplot(data, aes(x = iqr, y = Tau, fill = tissue)) + 
	geom_violin() + 
	theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1))
dev.off()

data <- melt(log2tpm.median.Tau)
pdf(file = paste(args[1], "tissueTau.pdf", sep = "."), width = 20, height = 10)
ggplot(data, aes(x = value)) + 
	geom_histogram(bins = 50) + 
	theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1))
dev.off()

data <- melt(log2tpm.median)
colnames(data) <- c("gene", "tissue", "mean")
pdf(file = paste(args[1], "median.pdf", sep = "."), width = 20, height = 10)
ggplot(data, aes(x = tissue, y = mean, fill = tissue)) + 
	geom_violin() + # scale_y_continuous(limits = c(0,6)) +
	theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1))
dev.off()

data <- melt(log2tpm.median.mean)
pdf(file = paste(args[1], "tissuemean.pdf", sep = "."), width = 20, height = 10)
ggplot(data, aes(x = value)) + 
	geom_histogram(bins = 1000) + # scale_x_continuous(limits = c(0,6)) +
	theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1))
dev.off()

data <- data.frame(log2tpm.median.mean, log2tpm.median.Tau)
colnames(data) <- c("mean", "Tau")
pdf(file = paste(args[1], "tissuemeanVStissueTau.pdf", sep = "."), width = 10, height = 10)
ggplot(data, aes(x = Tau, y = mean)) + 
	geom_point()+ 
	theme(legend.position = "none") 
dev.off()

for(i in 1:length(tissues)){
	data <- data.frame(log2tpm.median[,i], log2tpm.Tau[,i])
	colnames(data) <- c("median", "Tau")
	p <- ggplot(data, aes(x = Tau, y = median)) + geom_point() + ggtitle(tissues[i])
	pdf(file = paste(args[1], i, "medianVSTau.pdf", sep = "."), width = 10, height = 10)
	print(p)
	dev.off()
}

for(i in 1:length(tissues)){
	data <- data.frame(log2tpm.median[,i], log2tpm.sd[,i])
	colnames(data) <- c("median", "sd")
	p <- ggplot(data, aes(x = sd, y = median)) + geom_point() + ggtitle(tissues[i])
	pdf(file = paste(args[1], i, "medianVSsd.pdf", sep = "."), width = 10, height = 10)
	print(p)
	dev.off()
}

for(i in 1:length(tissues)){
	data <- data.frame(log2tpm.median[,i], log2tpm.mean[,i])
	colnames(data) <- c("median", "mean")
	p <- ggplot(data, aes(x = mean, y = median)) + geom_point() + ggtitle(tissues[i])
	pdf(file = paste(args[1], i, "medianVSmean.pdf", sep = "."), width = 10, height = 10)
	print(p)
	dev.off()
}

for(i in 1:length(tissues)){
	data <- data.frame(log2tpm.median[,i], log2tpm.iqr[,i])
	colnames(data) <- c("median", "iqr")
	p <- ggplot(data, aes(x = iqr, y = median)) + geom_point() + ggtitle(tissues[i])
	pdf(file = paste(args[1], i, "medianVSiqr.pdf", sep = "."), width = 10, height = 10)
	print(p)
	dev.off()
}

for(i in 1:length(tissues)){
	data <- data.frame(log2tpm.median[,i], log2tpm.outlier[,i])
	colnames(data) <- c("median", "outlier")
	p <- ggplot(data, aes(x = outlier, y = median)) + geom_point() + ggtitle(tissues[i])
	pdf(file = paste(args[1], i, "medianVSoutlier.pdf", sep = "."), width = 10, height = 10)
	print(p)
	dev.off()
}

for(i in 1:length(tissues)){
	data <- data.frame(log2tpm.median[,i], nullcount[,i])
	colnames(data) <- c("median", "nullcount")
	p <- ggplot(data, aes(x = nullcount, y = median)) + geom_point() + ggtitle(tissues[i])
	pdf(file = paste(args[1], i, "medianVSnullcount.pdf", sep = "."), width = 10, height = 10)
	print(p)
	dev.off()
}
