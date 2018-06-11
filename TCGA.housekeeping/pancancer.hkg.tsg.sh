# rclone sync /home1/04935/shaojf/stampede2/housekeeping_genes/both.pc.and.nc.genes/PANCAN.hkg.tsg/ mygoogle:hkg_tsg/both.pc.and.nc.genes/PANCAN.hkg.tsg/
# rclone sync /home1/04935/shaojf/stampede2/housekeeping_genes/ mygoogle:hkg_tsg/
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

head -1 hkg.tsg.cat > srt.sim.hkg.tsg.annotation
perl $myperl hkg.tsg.cat <(gunzip -c gencode.v23.annotation.gene.probeMap.gz | cut -f 1 -d"." | tail -n +2) 0 0 | cut -f 1,3- >> srt.sim.hkg.tsg.annotation

awk -F"\t" -vOFS="\t" '{if($3~/,/){$3="mixTSG"}if($3~/hkg/){$3="HKG"}if($2~/hkg/){$2="HKG"}print $1,$2,$3}' srt.sim.hkg.tsg.annotation > srt.sim.hkg.tsg.annotation.1
awk -F"\t" -vOFS="\t" '{if(NR>1 && $2!="HKG" && $2!="mixTSG" && $2!="/"){$2="singleTSG"}if(NR>1 && $3!="HKG" && $3!="mixTSG" && $3!="/"){$3="singleTSG"}print $1,$2,$3}' srt.sim.hkg.tsg.annotation.1 > srt.sim.hkg.tsg.annotation.2

myheatmap=/home1/04935/shaojf/myTools/BioinformaticsDaily/RVisualization/pheatmap.with.row.annotation.R

head -1 sim.v1.5.log2tpm.median.tsv > srt.sim.v1.5.log2tpm.median.tsv
perl $myperl sim.v1.5.log2tpm.median.tsv <(gunzip -c gencode.v23.annotation.gene.probeMap.gz | cut -f 1 -d"." | tail -n +2) 0 0 | cut -f 1,3- >> srt.sim.v1.5.log2tpm.median.tsv
Rscript $myheatmap srt.sim.v1.5.log2tpm.median.tsv srt.sim.hkg.tsg.annotation.2 N

head -1 sim.GTEx.log2tpm.median.tsv > srt.sim.GTEx.log2tpm.median.tsv
perl $myperl sim.GTEx.log2tpm.median.tsv <(gunzip -c gencode.v23.annotation.gene.probeMap.gz | cut -f 1 -d"." | tail -n +2) 0 0 | cut -f 1,3- >> srt.sim.GTEx.log2tpm.median.tsv
Rscript $myheatmap srt.sim.GTEx.log2tpm.median.tsv srt.sim.hkg.tsg.annotation.2 N

paste srt.sim.GTEx.log2tpm.median.tsv <(cut -f 2- srt.sim.v1.5.log2tpm.median.tsv) > srt.sim.GTEx.pancancer.tsv
Rscript $myheatmap srt.sim.GTEx.pancancer.tsv srt.sim.hkg.tsg.annotation.2 N

paste <(gunzip -c gencode.v23.annotation.gene.probeMap.gz | tail -n +2 | awk -vOFS="\t" '{print $3,$4,$5}') <(gunzip -c gencode.v23.annotation.gene.probeMap.gz | tail -n +2 | cut -f1 -d".") <(gunzip -c gencode.v23.annotation.gene.probeMap.gz | tail -n +2 | cut -f2,6) > sim.gencode.v23.annotation.gene.probeMap.bed
# bedtools intersect -wao -a sim.gencode.v23.annotation.gene.probeMap.bed -b hESC_domains_hg19.bed > sim.gencode.v23.annotation.gene.probeMap.hESC_domains_hg19
# awk '{a="/";if($7!="."){a="TAD";}print $4"\t"a}' sim.gencode.v23.annotation.gene.probeMap.hESC_domains_hg19 | uniq > sim.gencode.v23.annotation.gene.probeMap.hESC_domains_hg19.tad
# perl $myperl sim.gencode.v23.annotation.gene.probeMap.hESC_domains_hg19.tad srt.sim.hkg.tsg.annotation.2 0 0 | cut -f 1-3,5 > srt.sim.hkg.tsg.annotation.3
# sed -i '1s?/?hESC_domains?' srt.sim.hkg.tsg.annotation.3
# Rscript $myheatmap srt.sim.GTEx.pancancer.tsv srt.sim.hkg.tsg.annotation.3 N

awk -vOFS="\t" '{print $1,$2-5000,$2+5000"\n"$1,$3-5000,$3+5000}' hESC_domains_hg19.bed | bedtools sort -i | uniq > hESC_domains_hg19.boundary.bed
# awk -vOFS="\t" '{print $1,$2-10000,$2+10000"\n"$1,$3-10000,$3+10000}' hESC_domains_hg19.bed | bedtools sort -i | uniq > hESC_domains_hg19.boundary.bed
bedtools intersect -wao -a sim.gencode.v23.annotation.gene.probeMap.bed -b hESC_domains_hg19.boundary.bed > sim.gencode.v23.annotation.gene.probeMap.hESC_domains_hg19.boundary
awk '{a="/";if($7!="."){a="TADboundary";}print $4"\t"a}' sim.gencode.v23.annotation.gene.probeMap.hESC_domains_hg19.boundary | uniq > sim.gencode.v23.annotation.gene.probeMap.hESC_domains_hg19.tad.boundary
perl $myperl sim.gencode.v23.annotation.gene.probeMap.hESC_domains_hg19.tad.boundary srt.sim.hkg.tsg.annotation.2 0 0 | cut -f 1-3,5 > srt.sim.hkg.tsg.annotation.4
sed -i '1s?/?hESC_domains?' srt.sim.hkg.tsg.annotation.4
Rscript $myheatmap srt.sim.GTEx.pancancer.tsv srt.sim.hkg.tsg.annotation.4 N

cut -f 2-4 srt.sim.hkg.tsg.annotation.4 | sort | uniq -c 
  # 47997 /	/	/
  #     1 GTEx	pancancer	hESC_domains
  #   994 /	HKG	/
  #   544 HKG	/	/
  #  1754 HKG	HKG	/
  #   172 HKG	HKG	TADboundary
  #   133 /	HKG	TADboundary
  #    64 HKG	/	TADboundary
  #   755 /	mixTSG	/
  #   822 mixTSG	/	/
  #   255 mixTSG	mixTSG	/
  #    33 mixTSG	mixTSG	TADboundary
  #   277 mixTSG	singleTSG	/
  #    39 mixTSG	singleTSG	TADboundary
  #   124 /	mixTSG	TADboundary
  #   107 mixTSG	/	TADboundary
  #  1084 /	singleTSG	/
  #  1827 singleTSG	/	/
  #   221 singleTSG	mixTSG	/
  #    16 singleTSG	mixTSG	TADboundary
  #   411 singleTSG	singleTSG	/
  #    35 singleTSG	singleTSG	TADboundary
  #   110 /	singleTSG	TADboundary
  #   160 singleTSG	/	TADboundary
  #  2564 /	/	TADboundary

##### cluster
bedtools cluster -d 10000 -i sim.gencode.v23.annotation.gene.probeMap.bed > cluster.sim.gencode.v23
cut -f 7 cluster.sim.gencode.v23 | sort | uniq -c | sort -k1,1nr | awk '{print $2"\t"$1}' > cluster.sim.gencode.v23.stats
# perl $myperl sim.gencode.v23.annotation.gene.probeMap.bed <(cut -f 1 srt.sim.GTEx.pancancer.tsv.expressed.tsv | tail -n +2) 3 0 | cut -f 2- > srt.sim.expressed.bed
bedtools cluster -d 10000 -i srt.sim.expressed.bed > cluster.srt.sim.expressed
cut -f 7 cluster.srt.sim.expressed | sort | uniq -c | sort -k1,1nr | awk '{print $2"\t"$1}' > cluster.srt.sim.expressed.stats
# perl $myperl sim.gencode.v23.annotation.gene.probeMap.bed <(grep -w HKG srt.sim.hkg.tsg.annotation.4 | cut -f 1) 3 0 | cut -f 2- > srt.sim.HKG.bed
bedtools cluster -d 10000 -i srt.sim.HKG.bed > cluster.srt.sim.HKG
cut -f 7 cluster.srt.sim.HKG | sort | uniq -c | sort -k1,1nr | awk '{print $2"\t"$1}' > cluster.srt.sim.HKG.stats
perl $myperl cluster.sim.gencode.v23 srt.sim.hkg.tsg.annotation.4 3 0 > srt.sim.hkg.tsg.annotation.4.ingencodecluster
awk '$2!="/" || $3!="/"{print $NF}' srt.sim.hkg.tsg.annotation.4.ingencodecluster | sort | uniq -c | sort -k1,1nr | awk '{print $2"\t"$1}' > srt.sim.hkg.tsg.annotation.4.ingencodecluster.stats
awk '$2=="HKG" || $3=="HKG"{print $NF}' srt.sim.hkg.tsg.annotation.4.ingencodecluster | sort | uniq -c | sort -k1,1nr | awk '{print $2"\t"$1}' > srt.sim.hkg.tsg.annotation.4.ingencodecluster.HKG.stats
#### HKG for NOT clustered


#######
# rclone sync /home1/04935/shaojf/stampede2/housekeeping_genes/fantom.ep/ mygoogle:hkg_tsg/fantom.ep/
# cd /home1/04935/shaojf/stampede2/housekeeping_genes/fantom.ep
mynullR=/home1/04935/shaojf/myTools/BioinformaticsDaily/textProcess/output_rows_gt0.R
Rscript $mynullR human_permissive_enhancers_phase_1_and_2_expression_tpm_matrix.txt
awk '{print $1"\t"$1}' human_permissive_enhancers_phase_1_and_2_expression_tpm_matrix.txt.notnull.tsv | tail -n +2 | sed 's/:/\t/;s/-/\t/' > fantom.hk.enhancers.bed
# bedtools closest -D a -a <(awk -F"\t" -vOFS="\t" '{a=$4-1;b=$4;if($8=="-"){a=$5-1;b=$5}if(a<0){a=0;}print $3,a,b,$6,$2,$8}' hkg.tsg.srtbyPCA.transcript.bed | bedtools sort -i -) -b fantom.hk.enhancers.bed > hkg.tsg.srtbyPCA.transcript.fantom.hk.enhancers
# bedtools closest -D b -a fantom.hk.enhancers.bed -b <(awk -F"\t" -vOFS="\t" '{a=$4-1;b=$4;if($8=="-"){a=$5-1;b=$5}if(a<0){a=0;}print $3,a,b,$6,$2,$8}' hkg.tsg.srtbyPCA.transcript.bed | bedtools sort -i -) > fantom.hk.enhancers.hkg.tsg.srtbyPCA.transcript
bedtools intersect -wao -a <(awk -vOFS="\t" '{print $1,$2-10000,$3+10000,$4}' fantom.hk.enhancers.bed) -b <(awk -F"\t" -vOFS="\t" '{a=$4-1000;b=$4+1000;if($8=="-"){a=$5-1000;b=$5+1000}if(a<0){a=0;}print $3,a,b,$6,$2,$8}' hkg.tsg.srtbyPCA.transcript.bed | bedtools sort -i -) > fantom.hk.enhancers.hkg.tsg.srtbyPCA.transcript

####
grep -e "GM12878" -e "testis" -e "liver" -e "brain" -e "glio" fantom.enhancer.tissues.txt | cut -f 1 | tr "\n" "," | sed 's/,$//' | xargs -I mycol cut -f 1,mycol human_permissive_enhancers_phase_1_and_2_expression_tpm_matrix.txt > fantom.enh.gm.liver.testis
# cut -f1,22,33,47,70,71,86 srt.sim.GTEx.pancancer.tsv > GTEx.pancancer.gm.liver.testis
head_line srt.sim.GTEx.pancancer.tsv | grep -i -e "EBV" -e "testis" -e "liver" -e "brain" | cut -f 1 -d" " | tr "\n" "," | sed 's/,$//' | xargs -I mycol cut -f 1,mycol srt.sim.GTEx.pancancer.tsv > GTEx.pancancer.gm.liver.testis

cat <(cut -f 1-4 human_permissive_enhancers_phase_1_and_2.bed) <(cut -f 1-4 sim.gencode.v23.annotation.gene.probeMap.bed) | bedtools sort -i - > enh.gene.srt.bed
paste <(head -1 GTEx.pancancer.gm.liver.testis) <(grep -e "GM12878" -e "testis" -e "liver" -e "brain" -e "glio" fantom.enhancer.tissues.txt | cut -f 3 | tr "\n" "\t" | sed 's/\t$//') > enh.gene.srt.tsv
# perl $myperl GTEx.pancancer.gm.liver.testis <(cut -f 4 enh.gene.srt.bed) 0 0 | cut -f 1,3- | perl $myperl fantom.enh.gm.liver.testis /dev/stdin 0 0 | cut -f 1-5,7- >> enh.gene.srt.tsv
perl $myperl GTEx.pancancer.gm.liver.testis <(cut -f 4 enh.gene.srt.bed) 0 0 | cut -f 1,3- | perl $myperl fantom.enh.gm.liver.testis /dev/stdin 0 0 | cut -f 1-22,24- >> enh.gene.srt.tsv

cat <(cut -f 1-3 srt.sim.hkg.tsg.annotation.4 | awk '{print $0"\t/"}') <(cut -f 1 fantom.enh.gm.liver.testis | awk -vOFS="\t" '{print $1,"/","/","enhancer"}') > enh.gene.annos
sed -i '1s?/?fantom?' enh.gene.annos
head -1 enh.gene.annos > enh.gene.srt.annos
perl $myperl enh.gene.annos <(cut -f 4 enh.gene.srt.bed) 0 0 | cut -f 1,3- >> enh.gene.srt.annos


paste <(head -1 srt.sim.GTEx.pancancer.tsv) <(sed -n '171,173p;769,$p' fantom.enhancer.tissues.txt | cut -f 3 | tr "\n" "\t" | sed 's/\t$//') > all.enh.gene.srt.tsv
perl $myperl srt.sim.GTEx.pancancer.tsv <(cut -f 4 enh.gene.srt.bed) 0 0 | cut -f 1,3- | perl $myperl <(cut -f 1,171-173,769- human_permissive_enhancers_phase_1_and_2_expression_tpm_matrix.txt) /dev/stdin 0 0 | cut -f 1-91,93- >> all.enh.gene.srt.tsv
############
library(data.table)
library(pheatmap)
# args <- commandArgs(TRUE)
# args <- c("enh.gene.srt.tsv", "enh.gene.srt.annos")
args <- c("all.enh.gene.srt.tsv", "enh.gene.srt.annos")
input <- fread(args[1], sep = "\t", header = T, na.strings = "/")
class <- fread(args[2], sep = "\t", header = T, na.strings = "/")
scores <- data.matrix(input[,-1])
rownames(scores) <- 1:nrow(scores) # as.matrix(input[,1])
annosR <- class[,-1]
rownames(annosR) <- rownames(scores)
scores[is.na(scores) | scores < 0] <- 0
# enh.score <- log2(scores[,22:35] + 1)
enh.score <- log2(scores[,91:1153] + 1)
##
# gene.score <- scores[,1:21]
gene.score <- scores[,1:90]
gene.score[gene.score > 10] <- 10
enh.score[enh.score > 5] <- 5
##
scores <- data.frame(gene.score,enh.score)
range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}
scores.scale <- apply(scores, 2, range01)
# colors <- colorRampPalette(c("blue", "white", "red"))(100)
colors <- colorRampPalette(c("white", "blue"))(100)
# scores[scores > 10] <- 10
# scores.scale <- scores
# scores.scale[scores.scale > 7] <- 7
png(filename = paste(args[1], "pheatmap.png", sep = "."), width = 1500, height = 2000)
myplot <- pheatmap(scores.scale, scale = "none", annotation_row = annosR,
	show_rownames = F, show_colnames = T, color = colors, 
	cluster_cols = F, cluster_rows = F)
dev.off()
png(filename = paste(args[1], "pheatmap.gene.png", sep = "."), width = 1500, height = 2000)
myplot <- pheatmap(gene.score, scale = "none", annotation_row = annosR,
	show_rownames = F, show_colnames = T, color = colors, 
	cluster_cols = F, cluster_rows = F)
dev.off()
png(filename = paste(args[1], "pheatmap.enh.png", sep = "."), width = 1500, height = 2000)
myplot <- pheatmap(enh.score, scale = "none", annotation_row = annosR,
	show_rownames = F, show_colnames = T, color = colors, 
	cluster_cols = F, cluster_rows = F)
dev.off()
# png(filename = paste(args[1], "pheatmap.scale.png", sep = "."), width = 1500, height = 2000)
# myplot <- pheatmap(scores, scale = "column", annotation_row = annosR,
# 	show_rownames = F, show_colnames = T, color = colors, 
# 	cluster_cols = F, cluster_rows = F)
# dev.off()

scores.g <- scores.scale[(annosR[,1]=="HKG" & !is.na(annosR[,1])) | (annosR[,2]=="HKG"  & !is.na(annosR[,2])) | (annosR[,3]=="enhancer" & !is.na(annosR[,3])),]
annosR.g <- annosR[annosR$GTEx=="HKG" | annosR$pancancer=="HKG" | annosR$fantom=="enhancer",]
rownames(annosR.g) <- rownames(scores.g)
png(filename = paste(args[1], "pheatmap.HKG.png", sep = "."), width = 1500, height = 2000)
myplot <- pheatmap(scores.g, scale = "none", annotation_row = annosR.g,
	show_rownames = F, show_colnames = T, color = colors, 
	cluster_cols = F, cluster_rows = F)
dev.off()

scores.o <- scores.scale[!is.na(annosR[,1]) | !is.na(annosR[,2]) | (annosR[,3]=="enhancer" & !is.na(annosR[,3])),]
annosR.o <- annosR[GTEx!="NA" | pancancer!="NA" | fantom=="enhancer",]
rownames(annosR.o) <- rownames(scores.o)
png(filename = paste(args[1], "pheatmap.HKG.TSG.png", sep = "."), width = 1500, height = 2000)
myplot <- pheatmap(scores.o, scale = "none", annotation_row = annosR.o,
	show_rownames = F, show_colnames = T, color = colors, 
	cluster_cols = F, cluster_rows = F)
dev.off()

scores.a <- scores.scale[!is.na(annosR[,1]) | !is.na(annosR[,2]) | (annosR[,3]=="enhancer" & !is.na(annosR[,3])),c(14,16,21:24,34:35)]
annosR.a <- annosR[GTEx!="NA" | pancancer!="NA" | fantom=="enhancer",]
rownames(annosR.a) <- rownames(scores.a)
png(filename = paste(args[1], "pheatmap.HKG.TSG.1.png", sep = "."), width = 1000, height = 1000)
myplot <- pheatmap(scores.a, scale = "none", annotation_row = annosR.a,
	show_rownames = F, show_colnames = T, color = colors, 
	cluster_cols = F, cluster_rows = F)
dev.off()
############
