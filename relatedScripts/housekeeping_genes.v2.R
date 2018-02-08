# args <- commandArgs(TRUE)
library(pheatmap)
# library(data.table)
args <- c("HPA.tissue.tsv", "HPA.tissue", "HPA.normal.tissue.protein.tsv")
tpm <- as.matrix(read.table(args[1], sep = "\t", row.names = 1, header = T))
# log2tpm <- log2(tpm + 1.001)
# entropy <- apply(log2tpm, 1, function(x){p=x/sum(x);-sum(p*log2(p))})
nullcount <- apply(tpm, 1, function(x){length(x[x==0])})
tpm <- tpm[nullcount==0,]
dim(tpm) # [1] 11902    37
tpm <- log(tpm,2)
range(tpm) # [1] -3.321928 16.277242
ave.all <- apply(tpm, 1, mean)
sd.all <- apply(tpm, 1, sd)
fc.all <- log(tpm/ave.all,2)
largefccount <- apply(fc.all, 1, function(x){length(x[abs(x)>1])})
housekeepinggene <- tpm[largefccount==0 & sd.all<1.5,]
write.table(data.frame(rownames(housekeepinggene), 
	ave.all[largefccount==0 & sd.all<1.5], sd.all[largefccount==0 & sd.all<1.5], 
	housekeepinggene), 
	file = paste("housekeepinggene", args[2],"tsv", sep = "."), sep = "\t", row.names = FALSE, quote = FALSE)
save.image("HPA.RData")

protein <- as.matrix(read.table(args[3], sep = "\t", row.names = 1, header = T))
# tissueprotein <- colSums(protein)
notdetect <- apply(protein, 2, function(x){length(x[x==0])})
dist_methods <- c("correlation", "euclidean", "maximum", "manhattan", "canberra", "minkowski")
clust_methods <- c("ward.D", "ward.D2", "complete", "average", "mcquitty", "median", "centroid")
i <- 2
j <- 2
colors <- colorRampPalette(c("grey", "blue", "yellow", "red"))(4)

# png(filename = paste(args[2], dist_methods[i], clust_methods[j], "protein.png", sep = "."), width = 1800, height = 2000)
protein <- protein[,notdetect < nrow(protein)/2]
png(filename = paste(args[2], dist_methods[i], clust_methods[j], "protein.subset.tissue.png", sep = "."), width = 1800, height = 2000)
p1 <- pheatmap(protein, scale = "none", show_rownames = F, show_colnames = T, 
         color = colors, cluster_cols = T, clustering_distance_rows = dist_methods[i], 
         clustering_method = clust_methods[j],
         cutree_row = 5)
dev.off()
cluster <- cutree(p1$tree_row, k = 5)
write.table(data.frame(cluster[p1$tree_row$order], 
	rownames(protein)[p1$tree_row$order],
	protein[p1$tree_row$order,p1$tree_col$order]), 
	file = paste(args[2], "protein.subset.tissue.pheatmap.tsv", sep = "."), 
	sep = "\t", row.names = FALSE, quote = FALSE)
save.image("HPA.protein.RData")



#### visulization
source("https://bioconductor.org/biocLite.R")
biocLite("karyoploteR")
library(karyoploteR)
hk.input <- read.table("housekeepinggene.byGTEx.gene.meanlog2.sdlog2.txt", sep = "\t")
all.input <- read.table("all.mean.sd.txt", sep = "\t")

tad <- read.table("hESC_domains_hg19.bed", sep = "\t")
colnames(tad) <- c("chr","start", "end")
tad.gr <- makeGRangesFromDataFrame(df = tad)

# kp <- plotKaryotype(plot.type = 1, chromosomes = c("chr1", "chr2", "chr3"))
pdf(file = "housekeepinggene.byGTEx.karyoploteR.pdf", width = 12, height = 10)
kp <- plotKaryotype(plot.type = 1)
# kpDataBackground(kp, data.panel=1)
kpDataBackground(kp, r0 = 0, r1 = 0.8, data.panel = 1)
kpSegments(kp, chr = input$V1, x0 = input$V2, x1 = input$V3, 
           y0 = input$V5-input$V7, y1 = input$V5+input$V7, 
           col="blue", ymin = 0, ymax = 16, r0 = 0, r1 = 0.8)
dev.off()

# pdf(file = "housekeepinggene.byGTEx.karyoploteR.hESC_domains.pdf", width = 12, height = 10)
# pdf(file = "all.hESC_domains.rect.pdf", width = 12, height = 10)
pdf(file = "all.hk.hESC_domains.rect.pdf", width = 12, height = 10)
kp <- plotKaryotype(plot.type = 2)
kpDataBackground(kp, r0 = 0, r1 = 0.8, data.panel = 1)
kpDataBackground(kp, r0 = 0, r1 = 0.3, data.panel = 2)
# kpSegments(kp, chr = input$V1, x0 = input$V2, x1 = input$V3, 
#            y0 = input$V5-input$V7, y1 = input$V5+input$V7, 
#            col="blue", ymin = 0, ymax = 16, r0 = 0, r1 = 0.8)
input <- all.input
kpRect(kp, chr = input$V1, x0 = input$V2, x1 = input$V3, 
       y0 = input$V5-input$V7, y1 = input$V5+input$V7, 
       border=NA, col="green", ymin = 0, ymax = 16, r0 = 0, r1 = 0.8)
input <- hk.input
kpRect(kp, chr = input$V1, x0 = input$V2, x1 = input$V3, 
       y0 = input$V5-input$V7, y1 = input$V5+input$V7, 
       border=NA, col="blue", ymin = 0, ymax = 16, r0 = 0, r1 = 0.8)
kpPlotRegions(kp, r0 = 0, r1 = 0.3, data.panel = 2, data = tad.gr)
dev.off()
# kpLines()
