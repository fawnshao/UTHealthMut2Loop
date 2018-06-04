#!/bin/sh
myperl=/home1/04935/shaojf/myScripts/add_any_2files_together.pl
homermatrix=/home1/04935/shaojf/myTools/UTHealthMut2Loop/relatedScripts/make.matrix.from.homer.motif.file.pl
sumcpg=/home1/04935/shaojf/myTools/UTHealthMut2Loop/relatedScripts/sum_cpg_repeat_length.pl
# rclone sync /home1/04935/shaojf/stampede2/housekeeping_genes/both.pc.and.nc.genes/appris.tss/ mygoogle:hkg_tsg/both.pc.and.nc.genes/appris.tss/
bedtools intersect -wo -a hg19.cpgIslandExt.withCount.bed -b <(cut -f 3- hkg.tsg.srtbyPCA.transcript.bed | awk -vOFS="\t" '{a=$2-1000;b=$2+1000;if($6=="-"){a=$3-1000;b=$3+1000}if(a<0){a=0;}print $1,a,b,$4,"1000",$6}') > hkg.tsg.srtbyPCA.transcript.cpgIslandExt
perl $sumcpg hkg.tsg.srtbyPCA.transcript.cpgIslandExt > hkg.tsg.srtbyPCA.transcript.cpgIslandExt.length
for f in *.promoter.vert.motifs.txt
do
	pre=`echo $f | awk -F"." '{print $1}'`
	perl $homermatrix <(cut -f 1,22- $f | sed '1s?/Homer Distance From Peak(sequence,strand,conservation)??g') > $pre.promoter.vert.motifs.mat
done
# head_line HKG.promoter.vert.motifs.mat | grep -w -e "TATA-Box" -e "NRF" -e "Klf9" -e "ETS" -e "YY1" -e "E2A(bHLH),near_PU.1" | awk '{print $1}' | tr "\n" ","
cat *.promoter.vert.motifs.mat | cut -f 1,60,67,69,70,71,72,73,76,80,81,82,83,84,85,86,87,88,89,108,163,208,210,255,256,257,293,294,307,334 > hkg.tsg.srtbyPCA.transcript.candidates.motifs.mat
# head_line hkg.tsg.srtbyPCA.transcript.candidates.motifs.mat | cut -f 2 -d " " | cut -f 1 -d "/" | tr "\n" " "
echo "gene|transcript cpgIslandLength E2A(bHLH),near_PU.1 EHF(ETS) ELF1(ETS) ELF3(ETS) ELF5(ETS) Elk1(ETS) Elk4(ETS) ERG(ETS) ETS1(ETS) Ets1-distal(ETS) ETS:E-box(ETS,bHLH) ETS(ETS) ETS:RUNX(ETS,Runt) ETV1(ETS) Etv2(ETS) EWS:ERG-fusion(ETS) EWS:FLI1-fusion(ETS) Fli1(ETS) GABPA(ETS) Klf9(Zf) NRF1(NRF) NRF(NRF) PU.1:IRF8(ETS:IRF) PU.1-IRF(ETS:IRF) PU.1(ETS) SPDEF(ETS) SpiB(ETS) TATA-Box(TBP) YY1(Zf) HiChIP.GM12878 HiChIP.K562 HiChIP.Naive HiChIP.Th17 HiChIP.Treg" | tr " " "\t"> hkg.tsg.srtbyPCA.transcript.candidates.mat
paste <(cut -f 3 hkg.tsg.srtbyPCA.transcript.class | tail -n +2 | perl $myperl hkg.tsg.srtbyPCA.transcript.cpgIslandExt.length /dev/stdin 0 0 | cut -f 1,3) <(cut -f 3 hkg.tsg.srtbyPCA.transcript.class | tail -n +2 | perl $myperl hkg.tsg.srtbyPCA.transcript.candidates.motifs.mat /dev/stdin 0 0 | cut -f 3-) <(cut -f 3 hkg.tsg.srtbyPCA.transcript.class | tail -n +2 | perl $myperl hkg.tsg.srtbyPCA.transcript.HiChIP-Loop.count /dev/stdin 0 0 | cut -f 3-) | sed 's?/?0?g' >> hkg.tsg.srtbyPCA.transcript.candidates.mat
####################################################
# library(data.table)
# library(pheatmap)
# args <- c("hkg.tsg.srtbyPCA.transcript.candidates.mat", "hkg.tsg.srtbyPCA.transcript.class")
# input <- fread(args[1], sep = "\t", header = T, na.strings = "/")
# class <- fread(args[2], sep = "\t", header = T, na.strings = "/")
# scores <- data.matrix(input[,-1])
# rownames(scores) <- as.matrix(input[,1])
# annos <- class[,2]
# rownames(annos) <- rownames(scores)
# annos[Type == "hkg1" | Type == "hkg2" | Type == "hkg3" | Type == "hkg4"] <- "HKG"
# annos[Type != "HKG" & Type != "mixTSG"] <- "singleTSG"
# # scores[scores[,1] > 0, 1] <- 1
# myscale <- function(x){
# 	th <- as.vector(quantile(x, probs = 0.95))
# 	x[x > th] <- th
# 	return(x)
# }
# scores <- apply(scores, 2, myscale)
# nullcount <- apply(scores, 2, function(x){length(x[x==0])})
# range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}
# data.x <- apply(scores, 2, range01)
# colors <- colorRampPalette(c("white", "blue"))(10)

# png(filename = paste(args[1], "pheatmap.png", sep = "."), width = 1000, height = 1200)
# myplot <- pheatmap(data.x, scale = "none", annotation_row = annos,
# 	show_rownames = F, show_colnames = T, color = colors, 
# 	cluster_cols = T, cluster_rows = F)
# dev.off()
# png(filename = paste(args[1], "pheatmap.1.png", sep = "."), width = 1000, height = 1200)
# myplot <- pheatmap(data.x, scale = "none", annotation_row = annos,
# 	show_rownames = F, show_colnames = T, color = colors, 
# 	cluster_cols = T, cluster_rows = T)
# dev.off()
# png(filename = paste(args[1], "pheatmap.2.png", sep = "."), width = 1000, height = 1200)
# myplot <- pheatmap(data.x[, nullcount < 8500], scale = "none", annotation_row = annos,
# 	show_rownames = F, show_colnames = T, color = colors, 
# 	cluster_cols = T, cluster_rows = F)
# dev.off()
# png(filename = paste(args[1], "pheatmap.3.png", sep = "."), width = 1000, height = 1200)
# myplot <- pheatmap(data.x[, nullcount < 8500], scale = "none", annotation_row = annos,
# 	show_rownames = F, show_colnames = T, color = colors, 
# 	cluster_cols = T, cluster_rows = T)
# dev.off()
# png(filename = paste(args[1], "pheatmap.4.png", sep = "."), width = 1000, height = 1200)
# myplot <- pheatmap(data.x[, c(1,11,13,15,21:23,29,31:35)], scale = "none", annotation_row = annos,
# 	show_rownames = F, show_colnames = T, color = colors, 
# 	cluster_cols = T, cluster_rows = F)
# dev.off()
# png(filename = paste(args[1], "pheatmap.5.png", sep = "."), width = 1000, height = 1200)
# myplot <- pheatmap(data.x[, c(1,11,13,15,21:23,29,31:35)], scale = "none", annotation_row = annos,
# 	show_rownames = F, show_colnames = T, color = colors, 
# 	cluster_cols = T, cluster_rows = T)
# dev.off()

# x <- data.x
# x[x > 0] <- 1
# png(filename = paste(args[1], "pheatmap.bin.png", sep = "."), width = 1000, height = 1200)
# myplot <- pheatmap(x[, c(1,11,13,15,21:23,29,31:35)], scale = "none", annotation_row = annos,
# 	show_rownames = F, show_colnames = T, color = colors, 
# 	cluster_cols = T, cluster_rows = T)
# dev.off()
# annos.sim <- annos[1:8106,]
# rownames(annos.sim) <- rownames(scores)[1:8106]
# png(filename = paste(args[1], "pheatmap.bin.2type.png", sep = "."), width = 1000, height = 1200)
# myplot <- pheatmap(x[1:8106, c(1,11,13,15,21:23,29,31:35)], scale = "none", annotation_row = annos.sim,
# 	show_rownames = F, show_colnames = T, color = colors, 
# 	cluster_cols = T, cluster_rows = T)
# dev.off()
# png(filename = paste(args[1], "pheatmap.bin.2type.1.png", sep = "."), width = 1000, height = 1200)
# myplot <- pheatmap(x[1:8106, c(1,11,13,15,21:23,29,31:35)], scale = "none", annotation_row = annos.sim,
# 	show_rownames = F, show_colnames = T, color = colors, 
# 	cluster_cols = T, cluster_rows = F)
# dev.off()
# png(filename = paste(args[1], "pheatmap.bin.2type.2.png", sep = "."), width = 1000, height = 1200)
# myplot <- pheatmap(x[1:8106, c(1,11,13,21:23,29,31:35)], scale = "none", annotation_row = annos.sim,
# 	show_rownames = F, show_colnames = T, color = colors, 
# 	cluster_cols = T, cluster_rows = T)
# dev.off()
# png(filename = paste(args[1], "pheatmap.bin.2type.3.png", sep = "."), width = 1000, height = 1200)
# myplot <- pheatmap(x[1:8106, c(1,11,13,21:23,29,31:35)], scale = "none", annotation_row = annos.sim,
# 	show_rownames = F, show_colnames = T, color = colors, 
# 	cluster_cols = T, cluster_rows = F)
# dev.off()
# png(filename = paste(args[1], "pheatmap.bin.HKG.png", sep = "."), width = 1000, height = 1200)
# myplot <- pheatmap(x[1:4181, c(1,13,21:23,29,31:35)], scale = "none", 
# 	show_rownames = F, show_colnames = T, color = colors, 
# 	cluster_cols = T, cluster_rows = T)
# dev.off()

# tmp <- x[1:8106,31:35]
# a <- rowSums(tmp)
# a[a > 0] <- 1
# xx <- data.frame(x[1:8106, c(1,13,21,23,29)], a)
# # colnames(xx)[7] <- "loop"
# colnames(xx) <- c("cpgIslandLength", "ETS", "Klf9", "NRF", "TATA-Box", "loop")
# png(filename = paste(args[1], "pheatmap.bin.2type.sim.png", sep = "."), width = 1000, height = 1200)
# myplot <- pheatmap(xx, scale = "none", 
# 	show_rownames = F, show_colnames = T, color = colors, annotation_row = annos.sim,
# 	cluster_cols = T, cluster_rows = T)
# dev.off()
# png(filename = paste(args[1], "pheatmap.bin.2type.sim.1.png", sep = "."), width = 1000, height = 1200)
# myplot <- pheatmap(xx[order(xx[,1], xx[,5], xx[,6], xx[,2], xx[,3], xx[,4], decreasing = T), c(1,5,6,2,3,4)], scale = "none", 
# 	show_rownames = F, show_colnames = T, color = colors, annotation_row = annos.sim,
# 	cluster_cols = F, cluster_rows = F)
# dev.off()
# png(filename = paste(args[1], "pheatmap.bin.2type.sim.2.png", sep = "."), width = 1000, height = 1200)
# myplot <- pheatmap(xx[order(xx[,6], xx[,1], xx[,5], xx[,2], xx[,3], xx[,4], decreasing = T), c(6,1,5,2,3,4)], scale = "none", 
# 	show_rownames = F, show_colnames = T, color = colors, annotation_row = annos.sim,
# 	cluster_cols = F, cluster_rows = F)
# dev.off()


# tmp <- x[1:4181,31:35]
# a <- rowSums(tmp)
# a[a > 0] <- 1
# xx <- data.frame(x[1:4181, c(1,13,21,23,29)], a)
# # colnames(xx)[7] <- "loop"
# colnames(xx) <- c("cpgIslandLength", "ETS", "Klf9", "NRF", "TATA-Box", "loop")
# png(filename = paste(args[1], "pheatmap.bin.HKG.sim.png", sep = "."), width = 1000, height = 1200)
# myplot <- pheatmap(xx, scale = "none", 
# 	show_rownames = F, show_colnames = T, color = colors, 
# 	cluster_cols = T, cluster_rows = T)
# dev.off()
# png(filename = paste(args[1], "pheatmap.bin.HKG.sim.1.png", sep = "."), width = 1000, height = 1200)
# myplot <- pheatmap(xx[order(xx[,1], xx[,5], xx[,6], xx[,2], xx[,3], xx[,4], decreasing = T), c(1,5,6,2,3,4)], scale = "none", 
# 	show_rownames = F, show_colnames = T, color = colors, 
# 	cluster_cols = F, cluster_rows = F)
# dev.off()
# png(filename = paste(args[1], "pheatmap.bin.HKG.sim.2.png", sep = "."), width = 1000, height = 1200)
# myplot <- pheatmap(xx[order(xx[,6], xx[,1], xx[,5], xx[,2], xx[,3], xx[,4], decreasing = T), c(6,1,5,2,3,4)], scale = "none", 
# 	show_rownames = F, show_colnames = T, color = colors, 
# 	cluster_cols = F, cluster_rows = F)
# dev.off()

# tmp <- data.x[1:4181,31:35]
# a <- rowSums(tmp)/5
# xx <- data.frame(data.x[1:4181, c(1,13,21,23,29)], a)
# colnames(xx) <- c("cpgIslandLength", "ETS", "Klf9", "NRF", "TATA-Box", "loop")
# png(filename = paste(args[1], "pheatmap.HKG.sim.png", sep = "."), width = 1000, height = 1200)
# myplot <- pheatmap(xx, scale = "none", 
# 	show_rownames = F, show_colnames = T, color = colors, 
# 	cluster_cols = T, cluster_rows = T)
# dev.off()
# png(filename = paste(args[1], "pheatmap.HKG.sim.1.png", sep = "."), width = 1000, height = 1200)
# myplot <- pheatmap(xx[order(xx[,1], xx[,5], xx[,6], xx[,2], xx[,4], xx[,3], decreasing = T), c(1,5,6,2,4,3)], scale = "none", 
# 	show_rownames = F, show_colnames = T, color = colors, 
# 	cluster_cols = F, cluster_rows = F)
# dev.off()
####################################################

####################################################
# cd looping.motifs
##### homer automatically add 1 for start in bed
f=all.looped.enhancers.vert.motifs.txt
pre=`echo $f | sed 's/.vert.motifs.txt//'`
# cut -f 22- $f | head -1 | sed 's?/Homer Distance From Peak(sequence,strand,conservation)??g' | awk '{print "ID\t"$0}'> $pre.vert.motifs.mat
# perl $homermatrix <(cut -f 2-4,22- $f | tail -n +2 | awk -F"\t" -vOFS="\t" '{$2=$2-1;print $0}' | sed 's/chr//;s/\t/:/;s/\t/-/') >> $pre.vert.motifs.mat
perl $homermatrix <(cut -f 2-4,22- $f | sed '1s?/Homer Distance From Peak(sequence,strand,conservation)??g' | awk -F"\t" -vOFS="\t" '{$2=$2-1;print $0}' | sed 's/chr//;s/\t/:/;s/\t/-/') > $pre.vert.motifs.mat
# head_line all.looped.enhancers.vert.motifs.mat | grep -w -e "TATA-Box" -e "NRF" -e "Klf9" -e "ETS" -e "YY1" -e "E2A(bHLH),near_PU.1" -e "CTCF(" | awk '{print $1}' | tr "\n" ","
cut -f 1,45,60,67,69,70,71,72,73,76,80,81,82,83,84,85,86,87,88,89,108,163,208,210,255,256,257,293,294,307,334 all.looped.enhancers.vert.motifs.mat > all.looped.enhancers.candidates.motifs.mat

paste <(head -1 hkg.tsg.srtbyPCA.transcript.enhancers.mat) <(head -1 all.looped.enhancers.candidates.motifs.mat | cut -f 2- | tr "\t" "\n" | cut -d"/" -f1 | awk '{print "E."$1}' | tr "\n" "\t"  | sed 's/\t$//') <(head -1 hkg.tsg.srtbyPCA.transcript.candidates.mat | cut -f 2- | tr "\t" "\n" | cut -d"/" -f1 | awk '{print "P."$1}' | tr "\n" "\t"  | sed 's/\t$//') > hkg.tsg.srtbyPCA.transcript.enhancers.loop.ep.motif.mat
paste <(tail -n +2 hkg.tsg.srtbyPCA.transcript.enhancers.mat) <(tail -n +2 hkg.tsg.srtbyPCA.transcript.enhancers.mat | cut -f 1 | cut -f2 -d"%" | perl $myperl all.looped.enhancers.candidates.motifs.mat /dev/stdin 0 0 | cut -f 3-) <(tail -n +2 hkg.tsg.srtbyPCA.transcript.enhancers.mat | cut -f 1 | cut -f1 -d"%" | perl $myperl hkg.tsg.srtbyPCA.transcript.candidates.mat /dev/stdin 0 0 | cut -f 3-) | sed 's?/?0?g' >> hkg.tsg.srtbyPCA.transcript.enhancers.loop.ep.motif.mat
perl $myperl hkg.tsg.srtbyPCA.transcript.class <(cut -f 1 hkg.tsg.srtbyPCA.transcript.enhancers.loop.ep.motif.mat | sed 's/%/\t/') 2 0 | cut -f 1-2,4 | sed 's/\t/%/' > hkg.tsg.srtbyPCA.transcript.enhancers.annotation

####################################################
# library(data.table)
# library(pheatmap)
# args <- c("hkg.tsg.srtbyPCA.transcript.enhancers.loop.ep.motif.mat", "hkg.tsg.srtbyPCA.transcript.enhancers.annotation")
# input <- fread(args[1], sep = "\t", header = T, na.strings = "/")
# class <- fread(args[2], sep = "\t", header = T, na.strings = "/")
# scores <- data.matrix(input[,-1])
# rownames(scores) <- as.matrix(input[,1])
# annos <- class[,2]
# rownames(annos) <- rownames(scores)
# annos[Type == "hkg1" | Type == "hkg2" | Type == "hkg3" | Type == "hkg4"] <- "HKG"
# annos[Type != "HKG" & Type != "mixTSG"] <- "singleTSG"
# # scores[scores[,1] > 0, 1] <- 1
# myscale <- function(x){
# 	th <- as.vector(quantile(x, probs = 0.95))
# 	x[x > th] <- th
# 	return(x)
# }
# scores <- apply(scores, 2, myscale)
# nullcount <- apply(scores, 2, function(x){length(x[x==0])})
# range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}
# data.x <- apply(scores, 2, range01)
# colors <- colorRampPalette(c("white", "blue"))(10)

# x <- data.x
# x[x > 0] <- 1
# png(filename = paste(args[1], "pheatmap.bin.png", sep = "."), width = 1000, height = 1200)
# myplot <- pheatmap(x, scale = "none", annotation_row = annos,
# 	show_rownames = F, show_colnames = T, color = colors, 
# 	cluster_cols = T, cluster_rows = T)
# dev.off()
# png(filename = paste(args[1], "pheatmap.bin.1.png", sep = "."), width = 1000, height = 1200)
# myplot <- pheatmap(x, scale = "none", annotation_row = annos,
# 	show_rownames = F, show_colnames = T, color = colors, 
# 	cluster_cols = F, cluster_rows = T)
# dev.off()
# png(filename = paste(args[1], "pheatmap.bin.2.png", sep = "."), width = 1000, height = 1200)
# myplot <- pheatmap(x[,c(1:5,34,64,36,6,65,35,18,26,28,48,56,58)], scale = "none", annotation_row = annos,
# 	show_rownames = F, show_colnames = T, color = colors, 
# 	cluster_cols = F, cluster_rows = T)
# dev.off()
# png(filename = paste(args[1], "pheatmap.bin.3.png", sep = "."), width = 1000, height = 1200)
# myplot <- pheatmap(x[,c(34,64,36,6,65,35,18,26,28,48,56,58)], scale = "none", annotation_row = annos,
# 	show_rownames = F, show_colnames = T, color = colors, 
# 	cluster_cols = F, cluster_rows = T)
# dev.off()
# y <- x[annos$Type != "mixTSG",]
# annos.sim <- annos[Type != "mixTSG"]
# rownames(annos.sim) <- rownames(y)
# png(filename = paste(args[1], "pheatmap.bin.3.png", sep = "."), width = 800, height = 600)
# myplot <- pheatmap(y[,c(34,64,36,6,65,35,18,26,28,48,56,58)], scale = "none", annotation_row = annos.sim,
# 	show_rownames = F, show_colnames = T, color = colors, fontsize_col = 12,
# 	cluster_cols = F, cluster_rows = T)
# dev.off()
####################################################
