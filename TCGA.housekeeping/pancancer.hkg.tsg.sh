# rclone sync /home1/04935/shaojf/stampede2/housekeeping_genes/both.pc.and.nc.genes/PANCAN.hkg.tsg/ mygoogle:hkg_tsg/both.pc.and.nc.genes/PANCAN.hkg.tsg/
# find hkg and tsg in TCGA pancancer data 
# use cohort.TCGA.Pan-Cancer.PANCAN log2(RSEM + 1)
myperl=/home1/04935/shaojf/myTools/BioinformaticsDaily/textProcess/add_any_2files_together.pl
# input=tcga_RSEM_gene_tpm
# head -1 PANCAN_clinicalMatrix | cut -f 1,21-22,25 > $input.samples
# perl $myperl <(cut -f 1,21-22,25 PANCAN_clinicalMatrix) <(head_line $input | awk '{print $2}' | tail -n +2) 0 0 | cut -f 1,3- >> $input.samples
# awk -F"\t" '{print $1"\t"$3"|"$4}' $input.samples > $input.samples.sim
perl $myperl <(gunzip -c TCGA_phenotype_denseDataOnlyDownload.tsv.gz) <(head_line EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena | awk '{print $2}') 0 0 | awk -F"\t" '{print $1"\t"$5"|"$4}' > EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.phenotype
Rscript /home1/04935/shaojf/myTools/UTHealthMut2Loop/TCGA.housekeeping/HKG.TSG.v1.5.R
cut -f 1,150 v1.5.TSG.tsv | grep -v "," | sort -k2 | grep -v "tissueflags" > v1.5.TSG.sim
cut -f 1,150 v1.5.TSG.tsv | grep "," | sort -k2 >> v1.5.TSG.sim 
cut -f 1,2 v1.5.HKG.median.pheatmap.tsv | tail -n +2 | awk '{print $2"\thkg"$1}' | sort -k2 > v1.5.HKG.sim 
head -1 v1.5.log2tpm.median.tsv | cut -f2- | awk '{print "Gene\tFlags\t"$0}' > pancancer.hkg.tsg.tsv
perl $myperl v1.5.log2tpm.median.tsv <(cat v1.5.HKG.sim v1.5.TSG.sim) 0 0 | cut -f 1-2,4- >> pancancer.hkg.tsg.tsv

############
library(data.table)
library(pheatmap)
# args <- commandArgs(TRUE)
args <- c("pancancer.hkg.tsg.tsv")
input <- fread(args[1], sep = "\t", header = T)
scores <- data.matrix(input[,-c(1:2)])
rownames(scores) <- 1:nrow(scores) # as.matrix(input[,1])
annosR <- input[,2]
rownames(annosR) <- rownames(scores)
annosR[Flags == "hkg1" | Flags == "hkg2" | Flags == "hkg3" | Flags == "hkg4" | Flags == "hkg5"] <- "HKG"
annosR[grep(",",Flags)] <- "mixTSG"
annosR[Flags != "HKG" & Flags != "mixTSG"] <- "singleTSG"
colors <- colorRampPalette(c("blue", "white", "red"))(100)
scores[!is.na(scores) & scores > 15] <- 15
png(filename = paste(args[1], "pheatmap.png", sep = "."), width = 1500, height = 1200)
myplot <- pheatmap(scores, scale = "none", annotation_row = annosR, 
	show_rownames = F, show_colnames = T, color = colors, 
	cluster_cols = F, cluster_rows = F)
dev.off()
tsg <- scores[annosR$Flags != "HKG",]
png(filename = paste(args[1], "TSG.pheatmap.png", sep = "."), width = 1500, height = 1200)
myplot <- pheatmap(tsg, scale = "none", 
	show_rownames = F, show_colnames = T, color = colors, 
	cluster_cols = F, cluster_rows = F)
dev.off()
############
