# args <- commandArgs(TRUE)
# library(pheatmap)
library(data.table)
args <- c("GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct", "GTEx_sample.tissue.txt")

## read in data
# tpm <- as.matrix(read.table(args[1], sep = "\t", row.names = 1, header = T))
raw.table <- fread(args[1], sep = "\t", header = T)
tpm <- as.matrix(raw.table[,-1])
rownames(tpm) <- as.matrix(raw.table[,1])
info <- as.matrix(read.table(args[2], sep = "\t"))
tissues <- unique(info[,2])
tissues <- tissues[-c(7,24,25,53)]

## calculate
nullcount <- matrix(ncol = length(tissues), nrow = nrow(tpm))
log2tpm.mean <- matrix(ncol = length(tissues), nrow = nrow(tpm))
log2tpm.sd <- matrix(ncol = length(tissues), nrow = nrow(tpm))
largefc.count <- matrix(ncol = length(tissues), nrow = nrow(tpm))
for(i in 1:length(tissues)){
	z <- log2(tpm[, info[,2] == tissues[i]] + 1)
	nullcount[,i] <- apply(z, 1, function(x){length(x[x==0])})
	log2tpm.mean[,i] <- apply(z, 1, mean)
	log2tpm.sd[,i] <- apply(z, 1, sd)
	fc.all <- z/log2tpm.mean[,i]
	largefc.count[,i] <- apply(fc.all, 1, function(x){length(x[x > 2 | x < 0.5])})
}
colnames(log2tpm.mean) <- tissues
colnames(log2tpm.sd) <- tissues
colnames(largefc.count) <- tissues

log2tpm.mean.mean <- apply(log2tpm.mean, 1, mean)
log2tpm.mean.sd <- apply(log2tpm.mean, 1, sd)
fc.all <- log2tpm.mean/log2tpm.mean.mean
largefc.count.all <- apply(fc.all, 1, function(x){length(x[x > 2 | x < 0.5])})
nullcount.sum <- apply(nullcount, 1, sum)

## output
merge.sd <- data.frame(log2tpm.mean.mean, log2tpm.mean.sd, largefc.count.all, nullcount.sum,
	log2tpm.mean, log2tpm.sd, largefc.count, nullcount)
rownames(merge.sd) <- rownames(tpm)
log2tpm.sd.max <- apply(log2tpm.sd, 1, function(x){max(x[is.finite(x)], na.rm = T)})
largefc.count.max <- apply(largefc.count, 1, function(x){max(x[is.finite(x)], na.rm = T)})
housekeepinggene <- merge.sd[largefc.count.all == 0 & largefc.count.max == 0 & log2tpm.mean.sd < 1 & log2tpm.sd.max < 1,]
write.table(data.frame(rownames(housekeepinggene), housekeepinggene), 
	file = "housekeepinggene.byGTEx.v1.1.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

housekeepinggene <- merge.sd[largefc.count.all == 0 & largefc.count.max < 10 & log2tpm.mean.sd < 1.5 & log2tpm.sd.max < 1.5,]
write.table(data.frame(rownames(housekeepinggene), housekeepinggene), 
	file = "housekeepinggene.byGTEx.v1.1.loose.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

write.table(data.frame(rownames(merge.sd), merge.sd), 
	file = "housekeepinggene.byGTEx.v1.1.raw.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

save.image("HKG.v1.1.RData")


args <- c("housekeepinggene.byGTEx.v1.1.raw.sd","housekeepinggene.byGTEx.v1.1.raw.mean")
breaklists <- c(seq(0, 2, by = 0.01),seq(2.1, 5.7, by = 0.1))
colorn <- length(breaklists)
colors <- colorRampPalette(c("blue", "yellow", "red"))(colorn)
dist_methods <- c("correlation", "euclidean", "maximum", "manhattan", "canberra", "minkowski")
clust_methods <- c("ward.D", "ward.D2", "complete", "average", "mcquitty", "median", "centroid")
i <- 2
j <- 2

data <- data.frame(log2tpm.mean.sd,log2tpm.sd)
rownames(data) <- rownames(tpm)
png(filename = paste(args[1], dist_methods[i], clust_methods[j], "png", sep = "."), width = 1500, height = 1500)
p1 <- pheatmap(data, scale = "none", show_rownames = F, show_colnames = T, 
         color = colors, cluster_cols = T, breaks = breaklists,
         clustering_distance_cols = dist_methods[i], clustering_distance_rows = dist_methods[i], 
         clustering_method = clust_methods[j],
         cutree_row = 10)
dev.off()
cluster <- cutree(p1$tree_row, k = 10)
write.table(data.frame(cluster[p1$tree_row$order], 
	rownames(data)[p1$tree_row$order],
	data[p1$tree_row$order,p1$tree_col$order]), 
	file = paste(args[1], "pheatmap.tsv", sep = "."), 
	sep = "\t", row.names = FALSE, quote = FALSE)

png(filename = paste(args[1], dist_methods[i], clust_methods[j], "notree.png", sep = "."), width = 1000, height = 1500)
pheatmap(data[p1$tree_row$order,p1$tree_col$order], scale = "none", show_rownames = F, show_colnames = T, 
         color = colors, cluster_cols = F, cluster_rows = F, breaks = breaklists)
dev.off()

# breaklists <- c(seq(0, 5, by = 0.1),seq(5.5, 15, by = 0.5))
breaklists <- c(seq(0, 6, by = 0.05),seq(6.1, 16, by = 0.1))
colorn <- length(breaklists)
colors <- colorRampPalette(c("blue", "yellow", "red"))(colorn)
data <- data.frame(log2tpm.mean.mean,log2tpm.mean)
rownames(data) <- rownames(tpm)
png(filename = paste(args[2], dist_methods[i], clust_methods[j], "notree.png", sep = "."), width = 1000, height = 1500)
pheatmap(data[p1$tree_row$order,p1$tree_col$order], scale = "none", show_rownames = F, show_colnames = T, 
         color = colors, cluster_cols = F, cluster_rows = F, breaks = breaklists)
dev.off()






#### visulization
library(pheatmap)
args <- c("housekeepinggene.byGTEx.v1.1.raw.sd.tsv", 
	"housekeepinggene.byGTEx.v1.1.loose.sd.tsv")
breaklists <- c(seq(0, 2, by = 0.01),seq(2.1, 5.7, by = 0.1))
colorn <- length(breaklists)
colors <- colorRampPalette(c("blue", "yellow", "red"))(colorn)
dist_methods <- c("correlation", "euclidean", "maximum", "manhattan", "canberra", "minkowski")
clust_methods <- c("ward.D", "ward.D2", "complete", "average", "mcquitty", "median", "centroid")
i <- 2
j <- 2

data <- as.matrix(read.table(args[1], sep = "\t", row.names = 1, header = T))
png(filename = paste(args[1], dist_methods[i], clust_methods[j], "png", sep = "."), width = 1500, height = 1500)
p1 <- pheatmap(data, scale = "none", show_rownames = F, show_colnames = T, 
         color = colors, cluster_cols = T, breaks = breaklists,
         clustering_distance_cols = dist_methods[i], clustering_distance_rows = dist_methods[i], 
         clustering_method = clust_methods[j],
         cutree_row = 10)
dev.off()
cluster <- cutree(p1$tree_row, k = 10)
write.table(data.frame(cluster[p1$tree_row$order], 
	rownames(data)[p1$tree_row$order],
	data[p1$tree_row$order,p1$tree_col$order]), 
	file = paste(args[1], "pheatmap.tsv", sep = "."), 
	sep = "\t", row.names = FALSE, quote = FALSE)
save.image("HKG.v1.1.heatmap.RData")

# source("https://bioconductor.org/biocLite.R")
# biocLite("karyoploteR")
library(karyoploteR)
hk.input <- read.table("HKG.mean.sd.txt", sep = "\t")
all.input <- read.table("all.mean.sd.txt", sep = "\t")

tad <- read.table("hESC_domains_hg19.bed", sep = "\t")
colnames(tad) <- c("chr","start", "end")
tad.gr <- makeGRangesFromDataFrame(df = tad)

tad2 <- read.table("IMR90_domains_hg19.bed", sep = "\t")
colnames(tad2) <- c("chr","start", "end")
tad.gr2 <- makeGRangesFromDataFrame(df = tad2)

# kp <- plotKaryotype(plot.type = 1, chromosomes = c("chr1", "chr2", "chr3"))
pdf(file = "HKG.v1.1.byGTEx.pdf", width = 15, height = 12)
kp <- plotKaryotype(genome="hg19", plot.type = 2)
kpDataBackground(kp, r0 = 0, r1 = 0.8, data.panel = 1)
kpDataBackground(kp, r0 = 0, r1 = 0.5, data.panel = 2)
input <- hk.input
# kpSegments(kp, chr = as.matrix(input$V1), x0 = input$V2, x1 = input$V3, 
#            y0 = input$V5-input$V7, y1 = input$V5+input$V7, 
#            col="blue", ymin = 0, ymax = 16, r0 = 0, r1 = 0.8)
kpRect(kp, chr = as.matrix(input$V1), x0 = input$V2, x1 = input$V3, 
       y0 = input$V5-input$V7, y1 = input$V5+input$V7, 
       border=NA, col="blue", ymin = 0, ymax = 16, r0 = 0, r1 = 0.8)
kpPlotRegions(kp, r0 = 0, r1 = 0.2, data.panel = 2, data = tad.gr, col = "burlywood4")
kpPlotRegions(kp, r0 = 0.3, r1 = 0.5, data.panel = 2, data = tad.gr2, col = "coral")
dev.off()

# pdf(file = "housekeepinggene.byGTEx.karyoploteR.hESC_domains.pdf", width = 12, height = 10)
# pdf(file = "all.hESC_domains.rect.pdf", width = 12, height = 10)
pdf(file = "HKG.v1.1.byGTEx.vsALL.pdf", width = 15, height = 12)
kp <- plotKaryotype(genome="hg19", plot.type = 2)
kpDataBackground(kp, r0 = 0, r1 = 0.8, data.panel = 1)
kpDataBackground(kp, r0 = 0, r1 = 0.5, data.panel = 2)
# kpSegments(kp, chr = input$V1, x0 = input$V2, x1 = input$V3, 
#            y0 = input$V5-input$V7, y1 = input$V5+input$V7, 
#            col="blue", ymin = 0, ymax = 16, r0 = 0, r1 = 0.8)
input <- all.input
kpRect(kp, chr = as.matrix(input$V1), x0 = input$V2, x1 = input$V3, 
       y0 = input$V5-input$V7, y1 = input$V5+input$V7, 
       border=NA, col="green", ymin = 0, ymax = 16, r0 = 0, r1 = 0.8)
input <- hk.input
kpRect(kp, chr = as.matrix(input$V1), x0 = input$V2, x1 = input$V3, 
       y0 = input$V5-input$V7, y1 = input$V5+input$V7, 
       border=NA, col="blue", ymin = 0, ymax = 16, r0 = 0, r1 = 0.8)
kpPlotRegions(kp, r0 = 0, r1 = 0.2, data.panel = 2, data = tad.gr, col = "burlywood4")
kpPlotRegions(kp, r0 = 0.3, r1 = 0.5, data.panel = 2, data = tad.gr2, col = "coral")
dev.off()
# kpLines()


# ##
# x <- c(rep(1,70),rep(0.0001,26))
# p <- x/(sum(x))
# -sum(p*log2(p))

# x <- c(rep(1,30),rep(0.0001,66))
# p <- x/(sum(x))
# -sum(p*log2(p))