# args <- commandArgs(TRUE)
# library(pheatmap)
# library(data.table)
args <- c("GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct", "GTEx_sample.tissue.txt")
# exprs <- fread(args[1], header = T, sep = "\t")
# tpm <- exprs[,-1]
# rownames(tpm) <- as.matrix(exprs[,1])

exprs <- as.matrix(read.table(args[1], sep = "\t", row.names = 1, header = T))
info <- as.matrix(read.table(args[2], sep = "\t"))
tissues <- unique(info[,2])
# load("GTEx.RData")
dim(exprs) # [1] 56202 11281
tpm <- exprs[,info[,2]!="Whole Blood"]
dim(tpm) # [1] 56202 11281
nullcount <- apply(tpm, 1, function(x){length(x[x==0])})
tpm <- tpm[nullcount==0,]
dim(tpm) # [1] 11522 11281
tpm <- log(tpm,2)
range(tpm) # [1] -9.341917 17.679343
ave.all <- apply(tpm, 1, mean)
sd.all <- apply(tpm, 1, sd)
# fc.all <- apply(tpm, 1, function(x){log(x/mean(x),2)}) # need to t(fc.all)
fc.all <- log(tpm/ave.all,2)
largefccount <- apply(fc.all, 1, function(x){length(x[abs(x)>1])})
housekeepinggene <- tpm[largefccount==0 & sd.all<1.5,]
write.table(data.frame(rownames(housekeepinggene), 
	ave.all[largefccount==0 & sd.all<1.5], sd.all[largefccount==0 & sd.all<1.5], 
	housekeepinggene), 
	file = "housekeepinggene.byGTEx.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
housekeepinggene <- tpm[largefccount<100 & sd.all<1.5,]
write.table(data.frame(rownames(housekeepinggene), 
	ave.all[largefccount<100 & sd.all<1.5], sd.all[largefccount<100 & sd.all<1.5], 
	housekeepinggene), 
	file = "housekeepinggene.byGTEx.loose.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
# load("HKG.RData")
write.table(data.frame(rownames(tpm), 
	largefccount, ave.all, sd.all), 
	file = "housekeepinggene.byGTEx.raw.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
save.image("HKG.RData")




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


##
x <- c(rep(1,70),rep(0.0001,26))
p <- x/(sum(x))
-sum(p*log2(p))

x <- c(rep(1,30),rep(0.0001,66))
p <- x/(sum(x))
-sum(p*log2(p))