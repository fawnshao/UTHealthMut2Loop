library(data.table)
library(pheatmap)
library(ggplot2)
library(reshape2)

print("Reading Data")
args <- c("v1.3.2", "GTEx_sample.tissue.txt", "0.25", "0.25")
all.stats.file <- paste(args[1], "allstats.tsv", sep = ".")
all.stats <- fread(all.stats.file, , sep = "\t", header = T)
# all.stats <- data.frame(rownames(tpm), log2tpm.mean.mean, log2tpm.mean.Tau, 
# 	nullcount.sum, log2tpm.Tau.max, log2tpm.Tau.min, 
# 	log2tpm.mean, log2tpm.Tau)
info <- as.matrix(read.table(args[2], sep = "\t"))
tissues <- unique(info[,2])
tissues <- tissues[-c(7,24,25,53)]

tissueTau <- as.numeric(args[3])
sampleTau <- as.numeric(args[4])

genes <- as.matrix(all.stats[,1])
log2tpm.mean.mean <- as.matrix(all.stats[,2])
log2tpm.mean.Tau <- as.matrix(all.stats[,3])
nullcount.sum <- as.matrix(all.stats[,4])
log2tpm.Tau.max <- as.matrix(all.stats[,5])
log2tpm.Tau.min <- as.matrix(all.stats[,6])
log2tpm.mean <- as.matrix(all.stats[,7:55])
log2tpm.Tau <- as.matrix(all.stats[,56:104])
rownames(log2tpm.mean) <- genes
colnames(log2tpm.mean) <- tissues
rownames(log2tpm.Tau) <- genes
colnames(log2tpm.Tau) <- tissues


print("Plotting distribution for mean and Tau")

data <- melt(log2tpm.Tau)
colnames(data) <- c("gene", "tissue", "Tau")
png(filename = paste(args[1], "Tau.png", sep = "."), width = 2000, height = 1000)
ggplot(data, aes(x = tissue, y = Tau, fill = tissue)) + 
geom_violin() + 
theme(legend.position = "none") 
dev.off()

data <- melt(log2tpm.mean.Tau)
png(filename = paste(args[1], "tissueTau.png", sep = "."), width = 2000, height = 1000)
ggplot(data, aes(x = value)) + 
geom_histogram(bins = 1000) + 
theme(legend.position = "none") 
dev.off()

data <- melt(log2tpm.mean[log2tpm.mean.mean > 0.1,])
colnames(data) <- c("gene", "tissue", "mean")
png(filename = paste(args[1], "mean.png", sep = "."), width = 2000, height = 1000)
ggplot(data, aes(x = tissue, y = mean, fill = tissue)) + 
geom_violin() + scale_y_continuous(limits = c(0,6)) +
theme(legend.position = "none") 
dev.off()

data <- melt(log2tpm.mean.mean)
png(filename = paste(args[1], "tissuemean.png", sep = "."), width = 2000, height = 1000)
ggplot(data, aes(x = value)) + 
geom_histogram(bins = 1000) + scale_x_continuous(limits = c(0,6)) +
theme(legend.position = "none") 
dev.off()

data <- data.frame(log2tpm.mean.mean, log2tpm.mean.Tau)
colnames(data) <- c("mean", "Tau")
png(filename = paste(args[1], "tissuemeanVStissueTau.png", sep = "."), width = 1500, height = 1500)
ggplot(data, aes(x = Tau, y = mean)) + 
geom_point()+ 
theme(legend.position = "none") 
dev.off()

# scatterplot <- function(input, fileout, name){
# 	require(ggplot2)
# 	png(filename = fileout, width = 1500, height = 1500)
# 	ggplot(input, aes(x = Tau, y = mean)) + geom_point() + ggtitle(name)
# 	dev.off()
# }
for(i in 1:length(tissues)){
	data <- data.frame(log2tpm.mean[,i], log2tpm.Tau[,i])
	colnames(data) <- c("mean", "Tau")
	# scatterplot(data, paste(args[1], i, "meanVSTau.png", sep = "."), tissues[i])
	png(filename = paste(args[1], i, "meanVSTau.png", sep = "."), width = 1500, height = 1500)
	ggplot(data, aes(x = Tau, y = mean)) + geom_point() + ggtitle(tissues[i])
	dev.off()
}
for(i in 1:length(tissues)){
	print(paste(tissues[i], nrow(log2tpm.Tau[log2tpm.Tau[,i] < 0.15,]), sep = ": "))
}
# all.stats[grep("GAPDH", as.matrix(all.stats[,1])), 1]
# all.stats[Name.Description == "ENSG00000111640.10|GAPDH", 1]
# as.numeric(all.stats[ grep("\\|ACTB$",genes), ])

print("Looking for housekeeping genes")
housekeepinggene <- all.stats[nullcount.sum == 0 & !is.na(log2tpm.Tau.max) & log2tpm.Tau.max < sampleTau & !is.na(log2tpm.mean.Tau) & log2tpm.mean.Tau < tissueTau, ]
housekeepinggene.mean <- log2tpm.mean[nullcount.sum == 0 & !is.na(log2tpm.Tau.max) & log2tpm.Tau.max < sampleTau & !is.na(log2tpm.mean.Tau) & log2tpm.mean.Tau < tissueTau, ]
housekeepinggene.Tau <- log2tpm.Tau[nullcount.sum == 0 & !is.na(log2tpm.Tau.max) & log2tpm.Tau.max < sampleTau & !is.na(log2tpm.mean.Tau) & log2tpm.mean.Tau < tissueTau, ]
rownames(housekeepinggene.mean) <- as.matrix(housekeepinggene[,1])
rownames(housekeepinggene.Tau) <- as.matrix(housekeepinggene[,1])
# housekeepinggene[grep("GAPDH", as.matrix(housekeepinggene[,1])), 1]

## check ENSG00000204983.8|PRSS1 Pancreas
# tissuespecificgene[grep("PRSS1",as.matrix(tissuespecificgene[,1])),1]
# raw.table[grep("PRSS1",as.matrix(raw.table[,1])),1:5][10,]
print("Looking for tissue specific genes")
tissuespecificgene <- data.frame()
tissueflags <- c()
rindex <- c()
for(i in 1:nrow(all.stats)){
	if(!is.na(log2tpm.mean.Tau[i]) && log2tpm.mean.Tau[i] > 1 - tissueTau && !is.na(log2tpm.Tau.min[i]) && log2tpm.Tau.min[i] < sampleTau && min(nullcount[i,], na.rm = T) == 0){
		flag <- 0
		tflags <- c()
		for(j in 1:length(tissues)){
			temp <- log2tpm.mean[i,j]/fmax(log2tpm.mean[i,])
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
tissuespecificgene.mean <- log2tpm.mean[rindex, ]
tissuespecificgene.Tau <- log2tpm.Tau[rindex, ]
rownames(tissuespecificgene.mean) <- tissuespecificgene[,1]
rownames(tissuespecificgene.Tau) <- tissuespecificgene[,1]
