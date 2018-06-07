# rclone sync /home1/04935/shaojf/stampede2/housekeeping_genes/both.pc.and.nc.genes/PANCAN.hkg.tsg/ mygoogle:hkg_tsg/both.pc.and.nc.genes/PANCAN.hkg.tsg/
# find hkg and tsg in TCGA pancancer data 
# use cohort.TCGA.Pan-Cancer.PANCAN log2(RSEM + 1)
myperl=/home1/04935/shaojf/myTools/BioinformaticsDaily/textProcess/add_any_2files_together.pl
# input=tcga_RSEM_gene_tpm
# head -1 PANCAN_clinicalMatrix | cut -f 1,21-22,25 > $input.samples
# perl $myperl <(cut -f 1,21-22,25 PANCAN_clinicalMatrix) <(head_line $input | awk '{print $2}' | tail -n +2) 0 0 | cut -f 1,3- >> $input.samples
# awk -F"\t" '{print $1"\t"$3"|"$4}' $input.samples > $input.samples.sim
# perl $myperl <(gunzip -c TCGA_phenotype_denseDataOnlyDownload.tsv.gz) <(head_line EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena | awk '{print $2}') 0 0 | awk -F"\t" '{print $1"\t"$5"|"$4}' > EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.phenotype
# Rscript /home1/04935/shaojf/myTools/UTHealthMut2Loop/TCGA.housekeeping/HKG.TSG.v1.5.R

##### use tcga_RSEM_gene_tpm to generate the hkg and tsg list, as they also use gencode version
##### the only difference to GTEx hkg & tsg is nullcount.sum ==0 ====> nullcount.sum < 10 
##### tissues <- names(samplecount[samplecount > 20]), including normal, primary & metastic
cut -f 1,132 v1.5.TSG.tsv | grep -v "," | sort -k2 | grep -v "tissueflags" > v1.5.TSG.sim
cut -f 1,132 v1.5.TSG.tsv | grep "," | sort -k2 >> v1.5.TSG.sim 
cut -f 1,2 v1.5.HKG.median.pheatmap.tsv | tail -n +2 | awk '{print $2"\thkg"$1}' | sort -k2 > v1.5.HKG.sim 
head -1 v1.5.log2tpm.median.tsv | cut -f2- | awk '{print "Gene\tFlags\t"$0}' > pancancer.hkg.tsg.tsv
perl $myperl v1.5.log2tpm.median.tsv <(cat v1.5.HKG.sim v1.5.TSG.sim) 0 0 | cut -f 1-2,4- >> pancancer.hkg.tsg.tsv

############
# library(data.table)
# library(pheatmap)
# # args <- commandArgs(TRUE)
# args <- c("pancancer.hkg.tsg.tsv")
# input <- fread(args[1], sep = "\t", header = T)
# scores <- data.matrix(input[,-c(1:2)])
# rownames(scores) <- 1:nrow(scores) # as.matrix(input[,1])
# annosR <- input[,2]
# rownames(annosR) <- rownames(scores)
# annosR[Flags == "hkg1" | Flags == "hkg2" | Flags == "hkg3" | Flags == "hkg4" | Flags == "hkg5"] <- "HKG"
# annosR[grep(",",Flags)] <- "mixTSG"
# annosR[Flags != "HKG" & Flags != "mixTSG"] <- "singleTSG"
# colors <- colorRampPalette(c("blue", "white", "red"))(100)
# scores[!is.na(scores) & scores > 10] <- 10
# scores[!is.na(scores) & scores < 0] <- 0
# png(filename = paste(args[1], "pheatmap.png", sep = "."), width = 1500, height = 1200)
# myplot <- pheatmap(scores, scale = "none", annotation_row = annosR, 
# 	show_rownames = F, show_colnames = T, color = colors, 
# 	cluster_cols = F, cluster_rows = F)
# dev.off()
# tsg <- scores[annosR$Flags != "HKG",]
# png(filename = paste(args[1], "TSG.pheatmap.png", sep = "."), width = 1500, height = 1200)
# myplot <- pheatmap(tsg, scale = "none", 
# 	show_rownames = F, show_colnames = T, color = colors, 
# 	cluster_cols = F, cluster_rows = F)
# dev.off()
############

cat <(cut -f 1 -d"." pancancer.hkg.tsg.tsv | tail -n +2) <(cut -f 1 -d"." hkg.tsg.srtbyPCA.class | tail -n +2) | sort | uniq > hkg.tsg.id
echo "gene GTEx pancancer" | tr " " "\t" > hkg.tsg.cat
perl $myperl <(tail -n +2 hkg.tsg.srtbyPCA.class | sed 's/\./\t/' | cut -f 1,3) hkg.tsg.id 0 0 | cut -f 1,3 | perl $myperl <(tail -n +2 pancancer.hkg.tsg.tsv | sed 's/\./\t/' | cut -f 1,3) /dev/stdin 0 0 | cut -f 1-2,4 | sort -k2 >> hkg.tsg.cat
head -1 hkg.tsg.cat > hkg.tsg.cat.srt
# grep hkg hkg.tsg.cat >> hkg.tsg.cat.srt
# grep -v -e "hkg" hkg.tsg.cat | grep -v "," | tail -n +2 >> hkg.tsg.cat.srt
# perl $myperl hkg.tsg.cat.srt hkg.tsg.cat 0 0 | awk -F"\t" '$4=="/"' | cut -f 1-3 >> hkg.tsg.cat.srt
awk -F"\t" '$2~/hkg/' hkg.tsg.cat >> hkg.tsg.cat.srt
awk -F"\t" 'NR>1 && $2!~/hkg/ && $2!="mixTSG" && $2!="/"' hkg.tsg.cat >> hkg.tsg.cat.srt
awk -F"\t" '$2=="mixTSG"' hkg.tsg.cat >> hkg.tsg.cat.srt
awk -F"\t" '$2=="/" && $3~/hkg/' hkg.tsg.cat >> hkg.tsg.cat.srt
awk -F"\t" '$2=="/" && $3!~/hkg/ && $3!~/,/' hkg.tsg.cat >> hkg.tsg.cat.srt
awk -F"\t" '$2=="/" && $3!~/hkg/ && $3~/,/' hkg.tsg.cat >> hkg.tsg.cat.srt

head -1 GTEx.log2tpm.median.tsv > sim.GTEx.log2tpm.median.tsv
paste <(tail -n +2 GTEx.log2tpm.median.tsv | cut -f 1 -d".") <(tail -n +2 GTEx.log2tpm.median.tsv | cut -f 2-) >> sim.GTEx.log2tpm.median.tsv
head -1 v1.5.log2tpm.median.tsv > sim.v1.5.log2tpm.median.tsv
paste <(tail -n +2 v1.5.log2tpm.median.tsv | cut -f 1 -d".") <(tail -n +2 v1.5.log2tpm.median.tsv | cut -f 2-) >> sim.v1.5.log2tpm.median.tsv

#### vim manually, add Gene in title
perl $myperl sim.GTEx.log2tpm.median.tsv hkg.tsg.cat.srt 0 0 | cut -f 1-3,5- > hkg.tsg.cat.srt.GTEx
perl $myperl sim.v1.5.log2tpm.median.tsv hkg.tsg.cat.srt 0 0 | cut -f 1-3,5- > hkg.tsg.cat.srt.pancancer
paste hkg.tsg.cat.srt.GTEx <(cut -f 4- hkg.tsg.cat.srt.pancancer) > hkg.tsg.cat.srt.GTEx.pancancer
###########
# library(data.table)
# library(pheatmap)
# args <- c("hkg.tsg.cat.srt.GTEx.pancancer")
# input <- fread(args[1], sep = "\t", header = T, na.strings = "/")
# scores <- data.matrix(input[,-c(1:3)])
# rownames(scores) <- as.matrix(input[,1]) # as.matrix(input[,1]) # 1:nrow(scores)
# annosR <- input[,2:3]
# rownames(annosR) <- rownames(scores)
# annosR[GTEx == "hkg1" | GTEx == "hkg2" | GTEx == "hkg3" | GTEx == "hkg4", 1] <- "HKG"
# annosR[pancancer == "hkg1" | pancancer == "hkg2" | pancancer == "hkg3" | pancancer == "hkg4" | pancancer == "hkg5", 2] <- "HKG"
# annosR[grep(",",pancancer),2] <- "mixTSG"
# scores[is.na(scores) | scores < 0] <- 0

# gtex.mean <- apply(scores[,1:49], 1, mean)
# pancancer.mean <- apply(scores[,50:90], 1, mean)
# logfc <- pancancer.mean - gtex.mean
# tmp <- data.frame(scores, logfc, annosR)
# scores.diff <- tmp[abs(tmp$logfc) > 1 & !is.na(tmp$logfc),]
# scores.diff.hkg.1 <- scores.diff[scores.diff$GTEx == "HKG" & scores.diff$pancancer == "HKG" & !is.na(scores.diff$GTEx) & !is.na(scores.diff$pancancer) & scores.diff$logfc > 1,]
# scores.diff.hkg.2 <- scores.diff[scores.diff$GTEx == "HKG" & scores.diff$pancancer == "HKG" & !is.na(scores.diff$GTEx) & !is.na(scores.diff$pancancer) & scores.diff$logfc < -1,]
# scores.diff.hkg.3 <- scores.diff[scores.diff$GTEx == "HKG" & is.na(scores.diff$pancancer),]
# scores.diff.hkg.4 <- scores.diff[is.na(scores.diff$GTEx) & scores.diff$pancancer == "HKG",]
# scores.diff.hkg <- rbind.data.frame(scores.diff.hkg.1, scores.diff.hkg.2, scores.diff.hkg.3, scores.diff.hkg.4)[,1:90]
# scores.diff.hkg.anno <- rbind.data.frame(scores.diff.hkg.1, scores.diff.hkg.2, scores.diff.hkg.3, scores.diff.hkg.4)[,92:93]
# write.table(data.frame(rownames(scores.diff.hkg),scores.diff.hkg, scores.diff.hkg.anno), file = paste(args[1], "diff.hkg.tsv", sep = "."), sep = "\t", row.names = F, quote = F)

# # annosR[Flags != "HKG" & Flags != "mixTSG"] <- "singleTSG"
# colors <- colorRampPalette(c("blue", "white", "red"))(100)
# # scores[!is.na(scores) & scores > 10] <- 10
# # scores[is.na(scores)] <- 0
# scores[scores > 10] <- 10
# png(filename = paste(args[1], "pheatmap.png", sep = "."), width = 1500, height = 1200)
# myplot <- pheatmap(scores, scale = "none", annotation_row = annosR, 
# 	show_rownames = F, show_colnames = T, color = colors, 
# 	cluster_cols = F, cluster_rows = F)
# dev.off()
# pdf(file = paste(args[1], "pheatmap.pdf", sep = "."), width = 15, height = 12)
# myplot <- pheatmap(scores, scale = "none", annotation_row = annosR, 
# 	show_rownames = F, show_colnames = T, color = colors, 
# 	cluster_cols = F, cluster_rows = F)
# dev.off()

# scores.diff.hkg[scores.diff.hkg > 10] <- 10
# png(filename = paste(args[1], "diff.hkg.pheatmap.png", sep = "."), width = 1500, height = 1200)
# myplot <- pheatmap(scores.diff.hkg, scale = "none", annotation_row = scores.diff.hkg.anno, 
# 	show_rownames = T, show_colnames = T, color = colors, fontsize_row = 1,
# 	cluster_cols = F, cluster_rows = F)
# dev.off()
# pdf(file = paste(args[1], "diff.hkg.pheatmap.pdf", sep = "."), width = 15, height = 12)
# myplot <- pheatmap(scores.diff.hkg, scale = "none", annotation_row = scores.diff.hkg.anno, 
# 	show_rownames = T, show_colnames = T, color = colors, fontsize_row = 1,
# 	cluster_cols = F, cluster_rows = F)
# dev.off()
############

awk -F"\t" '$2~/hkg/ && $3=="/"' hkg.tsg.cat.srt | cut -f1 | grep -wf /dev/stdin hkg.tsg.srtbyPCA.class > hkg.lost.in.pancancer
# paste <(cut -f 1 -d"." gencode.v19.gene.anntotation.txt) <(cut -f 2-3 gencode.v19.gene.anntotation.txt) > sim.gencode.v19.gene.anntotation.txt
# paste <(cut -f 1 -d"." gencode.v23.gene.anntotation.txt) <(cut -f 2-3 gencode.v23.gene.anntotation.txt) > sim.gencode.v23.gene.anntotation.txt 
perl $myperl sim.gencode.v19.gene.anntotation.txt <(cut -f 1,92-93 hkg.tsg.cat.srt.GTEx.pancancer.diff.hkg.tsv) 0 0 | cut -f 1-3,5-6 | perl $myperl sim.gencode.v23.gene.anntotation.txt /dev/stdin 0 0 | cut -f 1-5,7-8 > hkg.tsg.cat.srt.GTEx.pancancer.diff.hkg.anntotation

# paste <(gunzip -c gencode.v23.annotation.gene.probeMap.gz | cut -f 1 -d".") <(gunzip -c gencode.v23.annotation.gene.probeMap.gz | cut -f 2) | tail -n +2 > sim.gencode.v23.annotation.gene.probeMap

head -1 sim.v1.5.log2tpm.median.tsv > srt.sim.v1.5.log2tpm.median.tsv
perl $myperl sim.v1.5.log2tpm.median.tsv <(gunzip -c gencode.v23.annotation.gene.probeMap.gz | cut -f 1 -d"." | tail -n +2) 0 0 | cut -f 1,3- >> srt.sim.v1.5.log2tpm.median.tsv
head -1 hkg.tsg.cat > srt.sim.hkg.tsg.annotation
perl $myperl hkg.tsg.cat <(gunzip -c gencode.v23.annotation.gene.probeMap.gz | cut -f 1 -d"." | tail -n +2) 0 0 | cut -f 1,3- >> srt.sim.hkg.tsg.annotation


