library(data.table)
library(pheatmap)
library(ggplot2)
library(reshape2)

print("Reading Data")
args <- c("v2", "GTEx_sample.tissue.txt")
outputpre <- args[1]

housekeepinggene.file <- paste(outputpre, "HKG.tsv", sep = ".")
housekeepinggene <- fread(housekeepinggene.file, sep = "\t", header = T)
tissuespecificgene.file <- paste(outputpre, "TSG.tsv", sep = ".")
tissuespecificgene <- fread(tissuespecificgene.file, sep = "\t", header = T)

info <- as.matrix(read.table(args[2], sep = "\t"))
tissues <- unique(info[,2])
tissues <- tissues[-c(7,24,25,31,53)]

housekeepinggene.median <- data.matrix(housekeepinggene[,4:(3+length(tissues))])
tissuespecificgene.median <- data.matrix(tissuespecificgene[,4:(3+length(tissues))])

rownames(housekeepinggene.median) <- as.matrix(housekeepinggene[,1])
rownames(tissuespecificgene.median) <- as.matrix(tissuespecificgene[,1])
colnames(housekeepinggene.median) <- tissues
colnames(tissuespecificgene.median) <- tissues

housekeepinggene.median.srt <- housekeepinggene.median[order(housekeepinggene[,2], decreasing = T),]
tissuespecificgene.median.srt <- tissuespecificgene.median[order(tissuespecificgene[,100], decreasing = F),]
datax <- rbind.data.frame(housekeepinggene.median.srt, tissuespecificgene.median.srt)
annosR <- data.frame(Type = c(rep("HKG", nrow(housekeepinggene.median)),as.matrix(tissuespecificgene[order(tissuespecificgene[,100], decreasing = F), 100])))
rownames(annosR) <- rownames(datax)

print("Plotting")
colors <- colorRampPalette(c("blue", "white", "red"))(100)
datax[datax > 10] <- 10
png(filename = paste(args[1], "pheatmap.png", sep = "."), width = 2500, height = 1500)
myplot <- pheatmap(datax, scale = "none", annotation_row = annosR,
	show_rownames = F, show_colnames = T, color = colors, 
	cluster_cols = F, cluster_rows = F)
dev.off()


