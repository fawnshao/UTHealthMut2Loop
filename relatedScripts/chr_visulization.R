# cut -f 1-2 HKG.v1.1.by.percentage.tsv | sed -n '2,$p' | perl ~/myScripts/add_any_2files_together.pl gencode.v19.gene.bed /dev/stdin 3 0 | awk -vOFS="\t" '{print $3,$4,$5,$1,$2,$8}' | bedtools sort -i - > housekeepinggene.byGTEx.gene.bed
# cut -f 1,3 HKG.v1.1.by.percentage.tsv | perl ~/myScripts/add_any_2files_together.pl /dev/stdin housekeepinggene.byGTEx.gene.bed 0 3 | cut -f 1-6,8 > HKG.mean.sd.txt
# cut -f 2-3 TSG.v1.1.by.percentage.tsv | sed -n '2,$p' | perl ~/myScripts/add_any_2files_together.pl gencode.v19.gene.bed /dev/stdin 3 0 | awk -vOFS="\t" '{print $3,$4,$5,$1,$2,$8}' | bedtools sort -i - > tissuespecificgene.byGTEx.gene.bed
# cut -f 2,4 TSG.v1.1.by.percentage.tsv | perl ~/myScripts/add_any_2files_together.pl /dev/stdin tissuespecificgene.byGTEx.gene.bed 0 3 | cut -f 1-6,8 > TSG.mean.sd.txt
# cut -f 1-3 ../housekeepinggene.byGTEx.v1.1.raw.tsv| perl ~/myScripts/add_any_2files_together.pl /dev/stdin gencode.v19.gene.bed 0 3 | awk -vOFS="\t" '{print $1,$2,$3,$4,$8,$6,$9}' | grep -v "/" > all.mean.sd.txt

library(karyoploteR)
hk.input <- read.table("HKG.mean.sd.txt", sep = "\t")
tsg.input <- read.table("TSG.mean.sd.txt", sep = "\t")
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

pdf(file = "TSG.v1.1.byGTEx.pdf", width = 15, height = 12)
kp <- plotKaryotype(genome="hg19", plot.type = 2)
kpDataBackground(kp, r0 = 0, r1 = 0.8, data.panel = 1)
kpDataBackground(kp, r0 = 0, r1 = 0.5, data.panel = 2)
input <- tsg.input
# kpSegments(kp, chr = as.matrix(input$V1), x0 = input$V2, x1 = input$V3, 
#            y0 = input$V5-input$V7, y1 = input$V5+input$V7, 
#            col="blue", ymin = 0, ymax = 16, r0 = 0, r1 = 0.8)
kpRect(kp, chr = as.matrix(input$V1), x0 = input$V2, x1 = input$V3, 
       y0 = input$V5-input$V7, y1 = input$V5+input$V7, 
       border=NA, col="sienna", ymin = 0, ymax = 16, r0 = 0, r1 = 0.8)
kpPlotRegions(kp, r0 = 0, r1 = 0.2, data.panel = 2, data = tad.gr, col = "burlywood4")
kpPlotRegions(kp, r0 = 0.3, r1 = 0.5, data.panel = 2, data = tad.gr2, col = "coral")
dev.off()

# pdf(file = "housekeepinggene.byGTEx.karyoploteR.hESC_domains.pdf", width = 12, height = 10)
# pdf(file = "all.hESC_domains.rect.pdf", width = 12, height = 10)
pdf(file = "HKG.TSG.v1.1.byGTEx.vsALL.pdf", width = 15, height = 12)
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
input <- tsg.input
kpRect(kp, chr = as.matrix(input$V1), x0 = input$V2, x1 = input$V3, 
       y0 = input$V5-input$V7, y1 = input$V5+input$V7, 
       border=NA, col="sienna", ymin = 0, ymax = 16, r0 = 0, r1 = 0.8)
kpPlotRegions(kp, r0 = 0, r1 = 0.2, data.panel = 2, data = tad.gr, col = "burlywood4")
kpPlotRegions(kp, r0 = 0.3, r1 = 0.5, data.panel = 2, data = tad.gr2, col = "coral")
dev.off()
