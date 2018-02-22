library(data.table)
library(pheatmap)
library(ggplot2)
library(reshape2)
# v1.5: add an arbitary cut off for HKG/TSG
# save.image("v1.4.vis.RData")
# load("v1.4.vis.RData")
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
args <- c("v1.5", "GTEx_sample.tissue.txt", "0.25", "0.5", "3")
# args <- c("v1.6", "GTEx_sample.tissue.txt", "0.25", "0.5", "10")
all.stats.file <- paste(args[1], "allstats.tsv", sep = ".")
all.stats <- fread(all.stats.file, sep = "\t", header = T)
nullcount.file <- paste(args[1], "nullcount.tsv", sep = ".")
nullcount <- as.matrix(fread(nullcount.file, sep = "\t", header = T)[,-1])
# all.stats <- data.frame(rownames(tpm), rawtpm.mean.mean, rawtpm.mean.Tau, 
# 	rawtpm.median.mean, rawtpm.median.Tau, 
# 	nullcount.sum, rawtpm.Tau.max, rawtpm.Tau.min, 
# 	rawtpm.mean, rawtpm.median, rawtpm.Tau)
info <- as.matrix(read.table(args[2], sep = "\t"))
tissues <- unique(info[,2])
tissues <- tissues[-c(7,24,25,53)]

outputpre <- args[1]
tissueTau <- as.numeric(args[3])
sampleTau <- as.numeric(args[4])
############
medianthreshold <- as.numeric(args[5])
############

genes <- as.matrix(all.stats[,1])
rawtpm.mean.mean <- as.matrix(all.stats[,2])
rawtpm.mean.Tau <- as.matrix(all.stats[,3])
rawtpm.median.mean <- as.matrix(all.stats[,4])
rawtpm.median.Tau <- as.matrix(all.stats[,5])
nullcount.sum <- as.matrix(all.stats[,6])
rawtpm.Tau.max <- as.matrix(all.stats[,7])
rawtpm.Tau.min <- as.matrix(all.stats[,8])
rawtpm.mean <- as.matrix(all.stats[,9:57])
rawtpm.median <- as.matrix(all.stats[,58:106])
rawtpm.Tau <- as.matrix(all.stats[,107:155])
rownames(rawtpm.mean) <- genes
colnames(rawtpm.mean) <- tissues
rownames(rawtpm.median) <- genes
colnames(rawtpm.median) <- tissues
rownames(rawtpm.Tau) <- genes
colnames(rawtpm.Tau) <- tissues
rownames(nullcount) <- genes
colnames(nullcount) <- tissues

# all.stats[grep("ISL1", genes),]

print("Plotting distribution for mean and Tau")

data <- melt(rawtpm.Tau)
colnames(data) <- c("gene", "tissue", "Tau")
png(filename = paste(args[1], "Tau.png", sep = "."), width = 2000, height = 1000)
ggplot(data, aes(x = tissue, y = Tau, fill = tissue)) + 
geom_violin() + 
theme(legend.position = "none") 
dev.off()

data <- melt(rawtpm.median.Tau)
png(filename = paste(args[1], "tissueTau.png", sep = "."), width = 2000, height = 1000)
ggplot(data, aes(x = value)) + 
geom_histogram(bins = 50) + 
theme(legend.position = "none") 
dev.off()

data <- melt(rawtpm.median[rawtpm.median.mean > 1,])
# data <- melt(rawtpm.median)
colnames(data) <- c("gene", "tissue", "mean")
png(filename = paste(args[1], "median.png", sep = "."), width = 2000, height = 1000)
ggplot(data, aes(x = tissue, y = mean, fill = tissue)) + 
geom_violin() + scale_y_continuous(limits = c(0,6)) +
theme(legend.position = "none") 
dev.off()

data <- melt(rawtpm.median.mean)
png(filename = paste(args[1], "tissuemean.png", sep = "."), width = 2000, height = 1000)
ggplot(data, aes(x = value)) + 
geom_histogram(bins = 1000) + scale_x_continuous(limits = c(0,6)) +
theme(legend.position = "none") 
dev.off()

data <- data.frame(rawtpm.median.mean, rawtpm.median.Tau)
colnames(data) <- c("mean", "Tau")
png(filename = paste(args[1], "tissuemeanVStissueTau.png", sep = "."), width = 1500, height = 1500)
ggplot(data, aes(x = Tau, y = mean)) + 
geom_point()+ 
theme(legend.position = "none") 
dev.off()

scatterplot <- function(x){
	require(ggplot2)
	print(tissues[x])
	data <- data.frame(rawtpm.median[,x], rawtpm.Tau[,x])
	colnames(data) <- c("median", "Tau")
	# png(filename = paste(args[1], i, "medianVSTau.png", sep = "."), width = 1500, height = 1500)
	p <- ggplot(data, aes(x = Tau, y = median)) + geom_point() + ggtitle(tissues[x])
	# dev.off()
	# return(p)
}
for(i in 1:length(tissues)){
	# data <- data.frame(rawtpm.median[,i], rawtpm.Tau[,i])
	# colnames(data) <- c("median", "Tau")
	# scatterplot(data, paste(args[1], i, "meanVSTau.png", sep = "."), tissues[i])
	png(filename = paste(args[1], i, "medianVSTau.png", sep = "."), width = 1500, height = 1500)
	# ggplot(data, aes(x = Tau, y = median)) + geom_point() + ggtitle(tissues[i])
	print(scatterplot(i))
	dev.off()
}

scatterplot2 <- function(x){
	require(ggplot2)
	print(tissues[x])
	data <- data.frame(rawtpm.mean[,x], rawtpm.Tau[,x])
	colnames(data) <- c("mean", "Tau")
	p <- ggplot(data, aes(x = Tau, y = mean)) + geom_point() + ggtitle(tissues[x])
}
for(i in 1:length(tissues)){
	png(filename = paste(args[1], i, "meanVSTau.png", sep = "."), width = 1500, height = 1500)
	print(scatterplot2(i))
	dev.off()
}

for(i in 1:length(tissues)){
	print(paste(tissues[i], nrow(rawtpm.Tau[rawtpm.Tau[,i] < 0.15,]), sep = ": "))
}
# all.stats[grep("GAPDH", as.matrix(all.stats[,1])), 1]
# all.stats[Name.Description == "ENSG00000111640.10|GAPDH", 1]
# as.numeric(all.stats[ grep("\\|ACTB$",genes), ])
# all.stats[rawtpm.median.Tau < 0.1,1:6]

print("Looking for housekeeping genes")
rawtpm.median.min <- apply(rawtpm.median, 1, function(x){min(x, na.rm = T)})
housekeepinggene <- all.stats[nullcount.sum == 0 & !is.na(rawtpm.median.Tau) & rawtpm.median.Tau < tissueTau & rawtpm.median.min > medianthreshold, ]
housekeepinggene <- all.stats[nullcount.sum == 0 & !is.na(rawtpm.Tau.max) & rawtpm.Tau.max < sampleTau & !is.na(rawtpm.median.Tau) & rawtpm.median.Tau < tissueTau & rawtpm.median.min > medianthreshold, ]
housekeepinggene.median <- rawtpm.median[nullcount.sum == 0 & !is.na(rawtpm.Tau.max) & rawtpm.Tau.max < sampleTau & !is.na(rawtpm.median.Tau) & rawtpm.median.Tau < tissueTau & rawtpm.median.min > medianthreshold, ]
housekeepinggene.Tau <- rawtpm.Tau[nullcount.sum == 0 & !is.na(rawtpm.Tau.max) & rawtpm.Tau.max < sampleTau & !is.na(rawtpm.median.Tau) & rawtpm.median.Tau < tissueTau & rawtpm.median.min > medianthreshold, ]
rownames(housekeepinggene.median) <- as.matrix(housekeepinggene[,1])
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
	if(!is.na(rawtpm.median.Tau[i]) && rawtpm.median.Tau[i] > 1 - tissueTau && !is.na(rawtpm.Tau.min[i]) && rawtpm.Tau.min[i] < sampleTau && min(nullcount[i,], na.rm = T) == 0){
		flag <- 0
		tflags <- c()
		for(j in 1:length(tissues)){
			temp <- rawtpm.median[i,j] / fmax(rawtpm.median[i,])
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
rownames(tissuespecificgene.median) <- as.matrix(tissuespecificgene[,1])
rownames(tissuespecificgene.Tau) <- as.matrix(tissuespecificgene[,1])

print("Writing output files")
write.table(housekeepinggene, file = paste(outputpre, "HKG.tsv", sep = "."), 
	sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(tissuespecificgene, tissueflags), file = paste(outputpre, "TSG.tsv", sep = "."), 
	sep = "\t", row.names = FALSE, quote = FALSE)

missed <- setdiff(x = genes[rawtpm.median.Tau < tissueTau], y = as.matrix(housekeepinggene[,1]))
missed.genes <- all.stats[genes %in% missed,]
write.table(missed.genes, file = paste(outputpre, "lowexpr.hkg.tsv", sep = "."), 
	sep = "\t", row.names = FALSE, quote = FALSE)


print("HKG/TSG counts")

data <- melt(data.frame(housekeepinggene.Tau,housekeepinggene[,5]))
colnames(data) <- c("tissue", "Tau")
png(filename = paste(args[1], "HKG.Tau.violin.png", sep = "."), width = 2000, height = 1000)
ggplot(data, aes(x = tissue, y = Tau, fill = tissue)) + 
geom_violin() + 
theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1)) 
dev.off()

data <- melt(data.frame(housekeepinggene.median,housekeepinggene[,4]))
colnames(data) <- c("tissue", "median")
png(filename = paste(args[1], "HKG.median.violin.png", sep = "."), width = 2000, height = 1000)
ggplot(data, aes(x = tissue, y = median, fill = tissue)) + 
geom_violin() + 
theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1)) 
dev.off()

data <- melt(data.frame(tissuespecificgene.Tau,tissuespecificgene[,5]))
colnames(data) <- c("tissue", "Tau")
png(filename = paste(args[1], "TSG.Tau.violin.png", sep = "."), width = 2000, height = 1000)
ggplot(data, aes(x = tissue, y = Tau, fill = tissue)) + 
geom_violin() + 
theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1)) 
dev.off()

data <- melt(data.frame(tissuespecificgene.median,tissuespecificgene[,4]))
colnames(data) <- c("tissue", "median")
png(filename = paste(args[1], "TSG.median.violin.png", sep = "."), width = 2000, height = 1000)
ggplot(data, aes(x = tissue, y = median, fill = tissue)) + 
geom_violin() + scale_y_continuous(limits = c(0,5)) +
theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1)) 
dev.off()

#### piecharts
tissueflags.count <- table(tissueflags)
lbls <- paste(names(tissueflags.count), ": ", tissueflags.count, sep="")
png(filename = paste(args[1], "TSG.piechart.png", sep = "."), width = 1000, height = 1000)
pie(tissueflags.count, labels = lbls)
dev.off()
lbls <- paste(names(tissueflags.count[tissueflags.count>30]), 
	": ", tissueflags.count[tissueflags.count>30], sep="")
png(filename = paste(args[1], "TSG.piechart.top.png", sep = "."), width = 1000, height = 1000)
pie(tissueflags.count[tissueflags.count>30], labels = lbls)
dev.off()
lbls <- paste(names(tissueflags.count[tissueflags.count>50]), 
	": ", tissueflags.count[tissueflags.count>50], sep="")
png(filename = paste(args[1], "TSG.piechart.top2.png", sep = "."), width = 1000, height = 1000)
pie(tissueflags.count[tissueflags.count>50], labels = lbls, radius = 0.8)
dev.off()

### mean vs Tau
data <- housekeepinggene[,4:5]
colnames(data) <- c("mean", "Tau")
png(filename = paste(args[1], "HKG.tissuemeanVStissueTau.png", sep = "."), width = 1500, height = 1500)
ggplot(data, aes(x = Tau, y = mean)) + 
geom_point() + scale_x_continuous(limits = c(0,1)) +
theme(legend.position = "none") 
dev.off()
data <- tissuespecificgene[,4:5]
colnames(data) <- c("mean", "Tau")
png(filename = paste(args[1], "TSG.tissuemeanVStissueTau.png", sep = "."), width = 1500, height = 1500)
ggplot(data, aes(x = Tau, y = mean)) + 
geom_point() + scale_x_continuous(limits = c(0,1)) +
theme(legend.position = "none") 
dev.off()
data <- rbind(housekeepinggene[,4:5], tissuespecificgene[,4:5])
colnames(data) <- c("mean", "Tau")
png(filename = paste(args[1], "HKG.TSG.tissuemeanVStissueTau.png", sep = "."), width = 1500, height = 1500)
ggplot(data, aes(x = Tau, y = mean)) + 
geom_point() + scale_x_continuous(limits = c(0,1)) +
theme(legend.position = "none") 
dev.off()
data1 <- data.frame(rawtpm.median.mean, rawtpm.median.Tau)
data2 <- housekeepinggene[,4:5]
data3 <- tissuespecificgene[,4:5]
colnames(data1) <- c("mean", "Tau")
colnames(data2) <- c("mean", "Tau")
colnames(data3) <- c("mean", "Tau")
png(filename = paste(args[1], "HKG.TSG.tissuemeanVStissueTau.1.png", sep = "."), width = 1500, height = 1500)
ggplot() + 
geom_point(data = data1, aes(x = Tau, y = mean), size = 0.3) + 
geom_point(data = data2, aes(x = Tau, y = mean), color = "red", size = 0.5) +
geom_point(data = data3, aes(x = Tau, y = mean), color = "blue", size = 0.5) +
scale_x_continuous(limits = c(0,1)) +
theme(legend.position = "none") 
dev.off()

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



##################
### need some pre-done inputs
##################
breaklists3 <- c(seq(0, 15, by = 0.05),seq(15.1, 16, by = 0.1))
colorn3 <- length(breaklists3)
colors3 <- colorRampPalette(c("blue", "yellow", "red"))(colorn3)

selected.tsg.file <- paste(args[1], "TSG.onlytissue.proteincoding.genes", sep = ".")
selected.tsg <- as.matrix(read.table(selected.tsg.file))[,1]
# selected.tsg.median <- data.frame(tissuespecificgene[rownames(tissuespecificgene.median) %in% selected.tsg, 4], 
# 	tissuespecificgene.median[rownames(tissuespecificgene.median) %in% selected.tsg,])
selected.tsg.median <- tissuespecificgene.median[rownames(tissuespecificgene.median) %in% selected.tsg,]
data <- selected.tsg.median
png(filename = paste(outputpre, "TSG.onlytissue.proteincoding.median.png", sep = "."), width = 1000, height = 1500)
p3 <- pheatmap(data, scale = "none", show_rownames = F, show_colnames = T, 
         color = colors3, cluster_cols = T, breaks = breaklists3,
         clustering_distance_cols = "euclidean", clustering_distance_rows = "euclidean", 
         clustering_method = "ward.D"
         )
dev.off()
cluster <- cutree(p3$tree_row, k = 40)
write.table(data.frame(cluster[p3$tree_row$order], 
	rownames(data)[p3$tree_row$order],
	data[p3$tree_row$order, p3$tree_col$order]), 
	file = paste(outputpre, "TSG.onlytissue.proteincoding.median.pheatmap.tsv", sep = "."), 
	sep = "\t", row.names = FALSE, quote = FALSE)

selected.tsg.file <- paste(args[1], "TSG.proteincoding.genes", sep = ".")
selected.tsg <- as.matrix(read.table(selected.tsg.file))[,1]
selected.tsg.median <- tissuespecificgene.median[rownames(tissuespecificgene.median) %in% selected.tsg,]
data <- selected.tsg.median
png(filename = paste(outputpre, "TSG.proteincoding.median.png", sep = "."), width = 1000, height = 1500)
p3 <- pheatmap(data, scale = "none", show_rownames = F, show_colnames = T, 
         color = colors3, cluster_cols = T, breaks = breaklists3,
         clustering_distance_cols = "euclidean", clustering_distance_rows = "euclidean", 
         clustering_method = "ward.D"
         )
dev.off()
cluster <- cutree(p3$tree_row, k = 40)
write.table(data.frame(cluster[p3$tree_row$order], 
	rownames(data)[p3$tree_row$order],
	data[p3$tree_row$order, p3$tree_col$order]), 
	file = paste(outputpre, "TSG.proteincoding.median.pheatmap.tsv", sep = "."), 
	sep = "\t", row.names = FALSE, quote = FALSE)

##
selected.hkg.file <- paste(args[1], "HKG.proteincoding.noMT.genes", sep = ".")
selected.hkg <- as.matrix(read.table(selected.hkg.file))[,1]
selected.hkg.median <- data.frame(housekeepinggene[rownames(housekeepinggene.median) %in% selected.hkg, 4], 
	housekeepinggene.median[rownames(housekeepinggene.median) %in% selected.hkg,])
selected.hkg.Tau <- data.frame(housekeepinggene[rownames(housekeepinggene.median) %in% selected.hkg, 5], 
	housekeepinggene.Tau[rownames(housekeepinggene.median) %in% selected.hkg,])
data <- selected.hkg.median
png(filename = paste(outputpre, "HKG.proteincoding.noMT.median.png", sep = "."), width = 1000, height = 1500)
p4 <- pheatmap(data, scale = "none", show_rownames = F, show_colnames = T, 
         color = colors3, cluster_cols = T, breaks = breaklists3,
         clustering_distance_cols = "euclidean", clustering_distance_rows = "euclidean", 
         clustering_method = "ward.D"
         )
dev.off()
cluster <- cutree(p4$tree_row, k = 5)
write.table(data.frame(cluster[p4$tree_row$order], 
	rownames(data)[p4$tree_row$order],
	data[p4$tree_row$order, p4$tree_col$order]), 
	file = paste(outputpre, "HKG.proteincoding.noMT.median.pheatmap.tsv", sep = "."), 
	sep = "\t", row.names = FALSE, quote = FALSE)
png(filename = paste(outputpre, "HKG.proteincoding.noMT.Tau.png", sep = "."), width = 1000, height = 1500)
data <- selected.hkg.Tau
pheatmap(data[p4$tree_row$order,p4$tree_col$order], scale = "none", 
	show_rownames = F, show_colnames = T, color = colors, 
	cluster_cols = F, cluster_rows = F, breaks = breaklists)
dev.off()
write.table(data.frame(cluster[p4$tree_row$order], 
	rownames(data)[p4$tree_row$order],
	data[p4$tree_row$order, p4$tree_col$order]), 
	file = paste(outputpre, "HKG.proteincoding.noMT.Tau.pheatmap.tsv", sep = "."), 
	sep = "\t", row.names = FALSE, quote = FALSE)
