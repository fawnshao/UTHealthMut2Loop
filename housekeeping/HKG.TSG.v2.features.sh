#!/bin/sh
# rclone sync /home1/04935/shaojf/stampede2/housekeeping_genes/v2.feature/ mygoogle:hkg_tsg/v2.feature/
# args <- c("v2", "GTEx_sample.tissue.txt", "0.15", "1.5", "2", "0.25")
# outputpre <- args[1]
# tau.threshold <- as.numeric(args[3])
# sd.threshold <- as.numeric(args[4])
# median.threshold <- as.numeric(args[5])
# tau.threshold2 <- as.numeric(args[6])

myperl=/home1/04935/shaojf/myTools/BioinformaticsDaily/textProcess/add_any_2files_together.pl
mymerge=/home1/04935/shaojf/myTools/BioinformaticsDaily/textProcess/merge.rows.by.specificID.for.single.file.pl

# cd /home1/04935/shaojf/stampede2/housekeeping_genes/v2.feature
echo "Gene Type" | sed 's/ /\t/' > v2.hkg.tsg.txt
cut -f 1-2 v2.HKG.tsv | tail -n +2 | sort -k2,2nr | cut -f 1 | awk '{print $0"\tHKG.1"}' >> v2.hkg.tsg.txt
cut -f 1-2 v2.HKG.2.tsv | tail -n +2 | sort -k2,2nr | cut -f 1 | awk '{print $0"\tHKG.2"}' >> v2.hkg.tsg.txt
cut -f 1,100 v2.TSG.tsv | tail -n +2 | grep -v "," | sort -k2 >> v2.hkg.tsg.txt
cut -f 1,100 v2.TSG.tsv | tail -n +2 | grep "," | sort -k2 >> v2.hkg.tsg.txt

perl $myperl <(sed 's/\t/|/' gencode.v19.gene.anntotation.txt) v2.hkg.tsg.txt 0 0 | cut -f 1-2,4 > v2.hkg.tsg.anntotation.txt
perl $myperl sim.gencode.vert.known.motifs.mat v2.hkg.tsg.anntotation.txt 0 0 | cut -f 1-3,5- > v2.hkg.tsg.vert.known.motifs.mat
# Rscript ~/myTools/UTHealthMut2Loop/housekeeping/heatmap.for.feature.selection.R v2.hkg.tsg.vert.known.motifs.mat
############################## chromHMM
# # split -l 500000 hg19.roadmap.25_imputed12marks.bed
# # DAXX TSS is not accurate
# for pre in HKG v2.hkg.2 Testis Spleen Liver Brain EBV protein_coding gencode.v19
# do
# 	for i in x??
# 	do
# 		bedtools intersect -wo -a <(awk -vOFS="\t" '{a=$2-1000;b=$2+1000;if($6=="-"){a=$3-1000;b=$3+1000}if(a<0){a=0;}print $1,a,b,$4,"1000",$6}' $pre.gene.bed) -b $i >> roadmap.chromHMM.$pre.txt
# 	done
# done
# pre=promoter.withcpg
# for i in x??
# do
# 	bedtools intersect -wo -a $pre.gene.bed -b $i >> roadmap.chromHMM.$pre.txt
# done

# for f in roadmap.chromHMM.*.txt
# do
# 	grep -v -i -e "cell" -e "blood" -e "cultured" $f > tissue.$f &
# done

# grep -v -i -e "cell" -e "blood" -e "cultured" EIDlegend.txt | wc -l
# 50

for pre in HKG v2.hkg.2 Testis Spleen Liver Brain EBV protein_coding gencode.v19
do
	bedtools intersect -wo -a <(awk -vOFS="\t" '{a=$2-1000;b=$2+1000;if($6=="-"){a=$3-1000;b=$3+1000}if(a<0){a=0;}print $1,a,b,$4,"1000",$6}' $pre.gene.bed) -b hg19.tissue.roadmap.activeenhancer.bed > roadmap.chromHMM.activeenhancer.$pre.txt
done
pre=promoter.withcpg
bedtools intersect -wo -a $pre.gene.bed -b hg19.tissue.roadmap.activeenhancer.bed > roadmap.chromHMM.activeenhancer.$pre.txt
for pre in HKG v2.hkg.2 Testis Spleen Liver Brain EBV protein_coding gencode.v19 promoter.withcpg
do
	perl $mymerge <(cut -f 4,11 roadmap.chromHMM.activeenhancer.$pre.txt) > roadmap.chromHMM.activeenhancer.$pre.sim.txt
done

for pre in HKG v2.hkg.2 Testis Spleen Liver Brain EBV protein_coding gencode.v19 promoter.withcpg
do
	echo -n $pre": "
	awk -F":" '{print $1}' roadmap.chromHMM.activeenhancer.$pre.sim.txt | awk '$2>20' | wc -l
done
# HKG promoters do not tend to be enhancers



cd test.top.hkg.tsg/
echo "GeneType Motif MotifCount GeneCount GeneWithMotif" | tr " " "\t" > homer.vert.motif.percentage.txt
for pre in HKG v2.hkg.2 Testis Spleen Liver Brain EBV protein_coding promoter.withcpg gencode.v19
do
	all=`tail -n +2 $pre.promoter.vert.motifs.txt | wc -l`
	for i in `seq 22 384`
	do
		motif=`head_line $pre.promoter.vert.motifs.txt | awk -vvar=$i '$1==var{print $2}' | awk -F"/" '{print $1"/"$2}'`
		mcount=`cut -f 1,$i $pre.promoter.vert.motifs.txt | tail -n +2 | awk -F"\t" '$2~/,/' | wc -l`
		perc=`echo $mcount/$all | bc -l`
		echo $pre $motif $mcount $all $perc | tr " " "\t" >> homer.vert.motif.percentage.txt
	done
done
perl ~/myTools/BioinformaticsDaily/textProcess/make_matrixwith3col_from_single_file.pl <(awk -F"\t" -vOFS="\t" '{print $2,$1,$5}' homer.vert.motif.percentage.txt | tail -n +2) > homer.vert.motif.percentage.mat
####################################################
# args <- c("homer.vert.motif.percentage.mat")
# library(pheatmap)
# library(data.table)
# input <- fread(args[1], sep = "\t", header = T, na.strings = "/")
# scores <- data.matrix(input[,-1])
# colors <- colorRampPalette(c("white", "red"))(10)
# rownames(scores) <- as.matrix(input[,1])[,1]

# pdf(file = paste(args[1], "pheatmap.pdf", sep = "."), width = 10, height = 16)
# myplot <- pheatmap(scores, scale = "none", fontsize_row = 4, 
# 	show_rownames = T, show_colnames = T, color = colors,
# 	cluster_cols = T, cluster_rows = T)
# dev.off()
# diff.scores <- scores[(scores[,3] > 0.2 | scores[,7] > 0.2) & (scores[,3]/(scores[,7]+1e-10) > 1.5 | scores[,3]/(scores[,7]+1e-10) < 0.67), ]
# pdf(file = paste(args[1], "diff.pheatmap.pdf", sep = "."), width = 10, height = 8)
# myplot <- pheatmap(diff.scores, scale = "none", fontsize_row = 8, 
# 	show_rownames = T, show_colnames = T, color = colors,
# 	cluster_cols = T, cluster_rows = T)
# dev.off()
# score.min <- apply(scores, 1, min)
# score.max <- apply(scores, 1, max)
# diff.scores.1 <- scores[score.max / (score.min + 1e-10) > 2 & score.max > 0.2, ]
# pdf(file = paste(args[1], "diff.1.pheatmap.pdf", sep = "."), width = 10, height = 8)
# myplot <- pheatmap(diff.scores.1, scale = "none", fontsize_row = 8, 
# 	show_rownames = T, show_colnames = T, color = colors,
# 	cluster_cols = T, cluster_rows = T)
# dev.off()
####################################################

awk -vOFS="\t" '{print $1,$2,$3"\n"$4,$5,$6}' FitHiChIP.GM12878_cell_line.SRR5831489.Interactions_FULL_ALLFeatures.bed | sort | uniq > FitHiChIP.GM12878_cell_line.anchors.bed
for pre in HKG Testis Spleen Liver
do
	bedtools intersect -wo -a <(awk -vOFS="\t" '{a=$2-5000;b=$2+5000;if($6=="-"){a=$3-5000;b=$3+5000}if(a<0){a=0;}print $1,a,b,$4,"1000",$6}' $pre.gene.bed) -b FitHiChIP.GM12878_cell_line.SRR5831489.Interactions_FULL_ALLFeatures.bed > $pre.FitHiChIP.GM12878.txt &
done




# awk '$4/($2+1) < 0.67 && $2 > 3' YY1del.salmon.tximport.abundance.tsv | cut -f 1,2,4 > batch1.tpm
# perl $myperl <(sed 's/|/\t/' v2.hkg.tsg.anntotation.txt) batch1.tpm 0 0 > YY1del.salmon.hkg.tsg
# perl $myperl <(cut -f 1,2,4 YY1del.salmon.tximport.abundance.tsv) <(sed 's/|/\t/' v2.hkg.tsg.anntotation.txt) 0 0 > hkg.tsg.YY1del.salmon.tsv
# cut -f 1-4,6-7 hkg.tsg.YY1del.salmon.tsv | sed 's/\t/|/' > hkg.tsg.YY1del.forheatmap.tsv
# awk '$15 < 0.01' YY1del.RUV.tsv | cut -f 1,11 > YY1del.RUV.deg
# perl $myperl <(sed 's/|/\t/' v2.hkg.tsg.anntotation.txt) YY1del.RUV.deg 0 0 > YY1del.RUV.deg.hkg.tsg
# perl $myperl <(cut -f 1,6-9 YY1del.RUV.tsv) <(sed 's/|/\t/' v2.hkg.tsg.anntotation.txt) 0 0 > hkg.tsg.YY1del.RUV.tsv
# cut -f 1-2,6- hkg.tsg.YY1del.RUV.tsv | sed 's/\t/|/' > hkg.tsg.YY1del.RUV.normcounts
# cut -f 1,6-9 YY1del.RUV.tsv > YY1del.RUV.normcounts
####################################################
# args <- c("hkg.tsg.YY1del.RUV.normcounts")
# args <- c("YY1del.RUV.normcounts")
# args <- c("hkg.tsg.YY1del.forheatmap.tsv")
# library(pheatmap)
# library(data.table)
# input <- fread(args[1], sep = "\t", header = T, na.strings = "/")
# scores <- log2(data.matrix(input[,4:5])+1)
# colors <- colorRampPalette(c("blue", "white", "red"))(10)
# rownames(scores) <- as.matrix(input[,1])[,1]
# annos_row <- as.data.frame(input[,2:3])
# rownames(annos_row) <- rownames(scores)
# scores[scores > 10] <- 10

# png(filename = paste(args[1], "pheatmap.png", sep = "."), width = 800, height = 1000)
# myplot <- pheatmap(scores, scale = "none", annotation_row = annos_row, 
# 	show_rownames = F, show_colnames = T, color = colors,
# 	cluster_cols = F, cluster_rows = F)
# dev.off()
####################################################

# cut -f 2 v2.top.hkg.tsg.txt | sort | uniq -c | sort -k1,1n
   #  160 Brain - Cerebellar Hemisphere,Brain - Cerebellum
   #  170 Cells - EBV-transformed lymphocytes
   #  216 Liver
   #  246 Spleen
   #  472 HKG
   # 1314 Testis

# bedtools intersect -wo -a <(awk -vOFS="\t" '{a=$2-1000;b=$2+1000;if($6=="-"){a=$3-1000;b=$3+1000}if(a<0){a=0;}print $1,a,b,$4,"1000",$6}' gencode.v19.gene.bed) -b hg19.cpgIslandExt.withCount.bed > gencode.v19.promoter.withcpg.txt
# cut -f 1-6 gencode.v19.promoter.withcpg.txt > promoter.withcpg.gene.bed
pre=promoter.withcpg
findMotifsGenome.pl $pre.gene.bed hg19 $pre.homer.motifs -size given -mknown /home1/04935/shaojf/myTools/HOMER/data/knownTFs/vertebrates/known.motifs 1> $pre.findMotifsGenome.txt 2>&1 &
annotatePeaks.pl $pre.gene.bed hg19 -m /home1/04935/shaojf/myTools/HOMER/data/knownTFs/vertebrates/known.motifs -size given -cpu 68 1> $pre.promoter.vert.motifs.txt 2> $pre.promoter.vert.motifs.log &

for pre in HKG Testis Spleen Liver
do
	awk -v var=$pre -F"\t" '$2==var{print $1}' v2.top.hkg.tsg.txt | perl $myperl gencode.v19.gene.bed /dev/stdin 3 0 | cut -f 2- > $pre.gene.bed
	findMotifsGenome.pl <(awk -vOFS="\t" '{a=$2-1000;b=$2+1000;if($6=="-"){a=$3-1000;b=$3+1000}if(a<0){a=0;}print $1,a,b,$4,"1000",$6}' $pre.gene.bed) hg19 $pre.homer.motifs -size given -mknown /home1/04935/shaojf/myTools/HOMER/data/knownTFs/vertebrates/known.motifs 1> $pre.findMotifsGenome.txt 2>&1 &
	annotatePeaks.pl <(awk -vOFS="\t" '{a=$2-1000;b=$2+1000;if($6=="-"){a=$3-1000;b=$3+1000}if(a<0){a=0;}print $1,a,b,$4,"1000",$6}' $pre.gene.bed) hg19 -m /home1/04935/shaojf/myTools/HOMER/data/knownTFs/vertebrates/known.motifs -size given -cpu 68 1> $pre.promoter.vert.motifs.txt 2> $pre.promoter.vert.motifs.log &
done
for pre in Brain EBV
do
	awk -v var=$pre -F"\t" '$2~var{print $1}' v2.top.hkg.tsg.txt | perl $myperl gencode.v19.gene.bed /dev/stdin 3 0 | cut -f 2- > $pre.gene.bed
	findMotifsGenome.pl <(awk -vOFS="\t" '{a=$2-1000;b=$2+1000;if($6=="-"){a=$3-1000;b=$3+1000}if(a<0){a=0;}print $1,a,b,$4,"1000",$6}' $pre.gene.bed) hg19 $pre.homer.motifs -size given -mknown /home1/04935/shaojf/myTools/HOMER/data/knownTFs/vertebrates/known.motifs 1> $pre.findMotifsGenome.txt 2>&1 &
	annotatePeaks.pl <(awk -vOFS="\t" '{a=$2-1000;b=$2+1000;if($6=="-"){a=$3-1000;b=$3+1000}if(a<0){a=0;}print $1,a,b,$4,"1000",$6}' $pre.gene.bed) hg19 -m /home1/04935/shaojf/myTools/HOMER/data/knownTFs/vertebrates/known.motifs -size given -cpu 68 1> $pre.promoter.vert.motifs.txt 2> $pre.promoter.vert.motifs.log &
done

# cut -f 1,103 protein_coding.promoter.vert.motifs.txt | awk -F"\t" '$2~/,/' | wc -l
# cut -f 1,20-21,103,354 HKG.promoter.vert.motifs.txt | more
# grep -i -e ets -e elk -e etv -e elf -e erg -e fli-1 -e ehf -e gabp -e spi GTRD.tf > ETS
while read tf
do
	grep $tf hg19.human_meta_clusters.bed > $tf.meta_clusters.bed &
done < ETS.keyword

while read tf
do
	findMotifsGenome.pl <(awk -F"\t" -vOFS="\t" '{a=$2-30;b=$3+30;if(a<0){a=0;}print $1,a,b,$4}' $tf.meta_clusters.bed) hg19 $tf.homer.motifs -size given -mknown /home1/04935/shaojf/myTools/HOMER/data/knownTFs/vertebrates/known.motifs 1> $tf.findMotifsGenome.txt 2>&1 &
	annotatePeaks.pl <(awk -F"\t" -vOFS="\t" '{a=$2-30;b=$3+30;if(a<0){a=0;}print $1,a,b,$4}' $tf.meta_clusters.bed) hg19 -m /home1/04935/shaojf/myTools/HOMER/data/knownTFs/vertebrates/known.motifs -size given -cpu 68 1> $tf.promoter.vert.motifs.txt 2> $tf.promoter.vert.motifs.log &
done < ETS.keyword

# for tf in EHF ERG ETV1 ETV6 GABP PU.1 Spi-B Fli-1
for tf in ERG PU.1 Fli-1
do
	findMotifsGenome.pl <(awk -F"\t" -vOFS="\t" '{a=$2-30;b=$3+30;if(a<0){a=0;}print $1,a,b,$4}' $tf.meta_clusters.bed) hg19 $tf.homer.motifs -size given -mknown /home1/04935/shaojf/myTools/HOMER/data/knownTFs/vertebrates/known.motifs 1> $tf.findMotifsGenome.txt 2>&1 &
	annotatePeaks.pl <(awk -F"\t" -vOFS="\t" '{a=$2-30;b=$3+30;if(a<0){a=0;}print $1,a,b,$4}' $tf.meta_clusters.bed) hg19 -m /home1/04935/shaojf/myTools/HOMER/data/knownTFs/vertebrates/known.motifs -size given -cpu 68 1> $tf.promoter.vert.motifs.txt 2> $tf.promoter.vert.motifs.log &
done
while read tf
do
	findMotifsGenome.pl <(awk -F"\t" -vOFS="\t" '{a=$2-30;b=$3+30;if(a<0){a=0;}print $1,a,b,$4}' $tf.meta_clusters.bed) hg19 $tf.homer.motifs -size given -mknown /home1/04935/shaojf/myTools/HOMER/data/knownTFs/vertebrates/known.motifs 1> $tf.findMotifsGenome.txt 2>&1 &
done < ETS.keyword

while read pre
do
	ets=`cut -f 1,103 $pre.promoter.vert.motifs.txt | tail -n +2 | awk -F"\t" '$2~/,/' | wc -l`
	all=`tail -n +2 $pre.promoter.vert.motifs.txt | wc -l`
	perc=`echo $ets/$all | bc -l`
	echo $pre $ets $all $perc
done < ETS.keyword
# c-Ets-1 10873 295133 .03684101743959503003
# c-Ets-2 93 1399 .06647605432451751250
# EHF 3906 16011 .24395727937043282743
# Elf-1 27830 174938 .15908493294767288982
# Elf-2 479 1148 .41724738675958188153
# Elf-3 17874 126997 .14074348212949912202
# Elf-4 1 9 .11111111111111111111
# Elf-5 280 936 .29914529914529914529
# Elk-1 2999 9297 .32257717543293535549
# Elk-4 2617 11336 .23085744530698659139
# Runtime error (func=(main), adr=3): Divide by zero
# ERG 0 0
# ETV1 8858 48215 .18371875972207819143
# ETV4 4631 19033 .24331424368202595492
# ETV5 524 6468 .08101422387136672850
# ETV6 4952 41867 .11827931306279408603
# ETV7 859 9820 .08747454175152749490
# Runtime error (func=(main), adr=3): Divide by zero
# Fli-1 0 0
# GABP 21788 142268 .15314758062248713695
# PU.1 43551 452126 .09632491827499413880
# Spi-B 5819 45219 .12868484486609611004


for tf in YY1
do
	grep $tf hg19.human_meta_clusters.bed > $tf.meta_clusters.bed
	findMotifsGenome.pl <(awk -F"\t" -vOFS="\t" '{a=$2-30;b=$3+30;if(a<0){a=0;}print $1,a,b,$4}' $tf.meta_clusters.bed) hg19 $tf.homer.motifs -size given -mknown /home1/04935/shaojf/myTools/HOMER/data/knownTFs/vertebrates/known.motifs 1> $tf.findMotifsGenome.txt 2>&1 &
	annotatePeaks.pl <(awk -F"\t" -vOFS="\t" '{a=$2-30;b=$3+30;if(a<0){a=0;}print $1,a,b,$4}' $tf.meta_clusters.bed) hg19 -m /home1/04935/shaojf/myTools/HOMER/data/knownTFs/vertebrates/known.motifs -size given -cpu 68 1> $tf.promoter.vert.motifs.txt 2> $tf.promoter.vert.motifs.log &
done
# YY1 13810 309694 .04459240411502967445

# cut -f 1 ../v2.HKG.2.tsv | tail -n +2 | awk '{print $0"\tHKG.2"}' > v2.hkg.2.txt
# perl $myperl gencode.v19.gene.bed v2.hkg.2.txt 3 0 | cut -f 3- > v2.hkg.2.gene.bed
for pre in HKG v2.hkg.2 Testis Spleen Liver Brain EBV protein_coding gencode.v19 promoter.withcpg
do
	ets=`cut -f 1,103 $pre.promoter.vert.motifs.txt | tail -n +2 | awk -F"\t" '$2~/,/' | wc -l`
	all=`tail -n +2 $pre.promoter.vert.motifs.txt | wc -l`
	perc=`echo $ets/$all | bc -l`
	echo $pre $ets $all $perc
done
# HKG 252 472 .53389830508474576271
# v2.hkg.2 1479 2933 .50426184793726559836
# Testis 383 1314 .29147640791476407914
# Spleen 47 246 .19105691056910569105
# Liver 45 216 .20833333333333333333
# Brain 68 335 .20298507462686567164
# EBV 73 228 .32017543859649122807
# protein_coding 6815 20345 .33497173752764807077
# gencode.v19 14054 57820 .24306468350051885160
# promoter.withcpg 7356 19457 .37806444981240684586

for pre in HKG v2.hkg.2 Testis Spleen Liver Brain EBV protein_coding gencode.v19 promoter.withcpg
do
	yy1=`cut -f 1,354 $pre.promoter.vert.motifs.txt | tail -n +2 | awk -F"\t" '$2~/,/' | wc -l`
	all=`tail -n +2 $pre.promoter.vert.motifs.txt | wc -l`
	perc=`echo $yy1/$all | bc -l`
	echo $pre $yy1 $all $perc
done
# HKG 106 472 .22457627118644067796
# v2.hkg.2 580 2933 .19774974428912376406
# Testis 119 1314 .09056316590563165905
# Spleen 7 246 .02845528455284552845
# Liver 5 216 .02314814814814814814
# Brain 11 335 .03283582089552238805
# EBV 9 228 .03947368421052631578
# protein_coding 1952 20345 .09594494961907102482
# gencode.v19 3996 57820 .06911103424420615703
# promoter.withcpg 2086 19457 .10721077247263195765

for pre in HKG v2.hkg.2 Testis Spleen Liver Brain EBV protein_coding gencode.v19 promoter.withcpg
do
	rand=`cut -f 1,269 $pre.promoter.vert.motifs.txt | tail -n +2 | awk -F"\t" '$2~/,/' | wc -l`
	all=`tail -n +2 $pre.promoter.vert.motifs.txt | wc -l`
	perc=`echo $rand/$all | bc -l`
	echo $pre $rand $all $perc
done
# 269 PRDM14
# HKG 107 472 .22669491525423728813
# Testis 303 1314 .23059360730593607305
# Spleen 67 246 .27235772357723577235
# Liver 44 216 .20370370370370370370
# Brain 79 335 .23582089552238805970
# EBV 67 228 .29385964912280701754
# protein_coding 4325 20345 .21258294421233718358
# gencode.v19 12211 57820 .21118989968868903493
# v2.hkg.2 630 2933 .21479713603818615751

# 60 c-myc
# HKG 185 472 .39194915254237288135
# Testis 446 1314 .33942161339421613394
# Spleen 58 246 .23577235772357723577
# Liver 56 216 .25925925925925925925
# Brain 119 335 .35522388059701492537
# EBV 79 228 .34649122807017543859
# protein_coding 7729 20345 .37989678053575817154
# gencode.v19 16509 57820 .28552404012452438602
# v2.hkg.2 1311 2933 .44698261166041595635



########## roadmap histone
f=H3K4me3.Brain.Liver.Spleen
computeMatrix scale-regions \
	-R HKG.gene.bed v2.hkg.2.gene.bed Brain.gene.bed Liver.gene.bed Spleen.gene.bed \
	-S bigwigs/*Brain*H3K4me3*.bigwig bigwigs/*Liver*H3K4me3*.bigwig bigwigs/*Spleen*H3K4me3*.bigwig \
	-p 48 \
	-b 3000 -a 3000 \
	--regionBodyLength 5000 \
	--skipZeros -o ${f}_scaled.gz \
	--outFileNameMatrix ${f}_scaled.tab \
	--outFileSortedRegions ${f}_genes.bed
labs1=`ls bigwigs/*Brain*H3K4me3*.bigwig| awk -F"/" '{print $NF}' | sed 's/.bigwig//' | tr "\n" " "`
labs2=`ls bigwigs/*Liver*H3K4me3*.bigwig | awk -F"/" '{print $NF}' | sed 's/.bigwig//' | tr "\n" " "`
labs3=`ls bigwigs/*Spleen*H3K4me3*.bigwig | awk -F"/" '{print $NF}' | sed 's/.bigwig//' | tr "\n" " "`
plotHeatmap -m ${f}_scaled.gz \
	-o ${f}_scaled.plotHeatmap.pdf \
	--colorMap Greens --whatToShow 'heatmap and colorbar' \
	--plotFileFormat pdf --samplesLabel $labs1 $labs2 $labs3 \
	--heatmapWidth 10 \
	--perGroup &		

plotHeatmap -m ${f}_scaled.gz \
	-o ${f}_scaled.plotHeatmap.1.pdf \
	--colorMap Greens --whatToShow 'heatmap and colorbar' \
	--plotFileFormat pdf --samplesLabel $labs1 $labs2 $labs3 \
	--heatmapWidth 10 --heatmapHeight 20 \
	--perGroup &

plotHeatmap -m ${f}_scaled.gz \
    -o ${f}_scaled.plotHeatmap.2.pdf \
    --colorMap Greens --whatToShow 'heatmap and colorbar' \
    --plotFileFormat pdf --samplesLabel $labs1 $labs2 $labs3 \
    --heatmapWidth 10 --heatmapHeight 20 \
	--regionsLabel HKG.tier1 HKG.tier2 Brain Liver Spleen &
wait

for f in Fetal_Brain BI.Brain UCSF-UBC.Brain Liver Spleen
do
	computeMatrix scale-regions \
		-R HKG.gene.bed v2.hkg.2.gene.bed Brain.gene.bed Liver.gene.bed Spleen.gene.bed \
		-S bigwigs/*${f}*H3K4me3*.bigwig \
		-p 48 \
		-b 3000 -a 3000 \
		--regionBodyLength 5000 \
		--skipZeros -o H3K4me3_${f}_scaled.gz \
		--outFileNameMatrix H3K4me3_${f}_scaled.tab \
		--outFileSortedRegions H3K4me3_${f}_genes.bed
	labs=`ls bigwigs/*${f}*H3K4me3*.bigwig| awk -F"/" '{print $NF}' | sed 's/.bigwig//' | tr "\n" " "`
	plotHeatmap -m H3K4me3_${f}_scaled.gz \
		-o H3K4me3_${f}_scaled.plotHeatmap.pdf \
		--colorMap Greens --whatToShow 'heatmap and colorbar' \
		--plotFileFormat pdf --samplesLabel $labs \
		--heatmapWidth 10 \
		--perGroup &		

	plotHeatmap -m H3K4me3_${f}_scaled.gz \
		-o H3K4me3_${f}_scaled.plotHeatmap.1.pdf \
		--colorMap Greens --whatToShow 'heatmap and colorbar' \
		--plotFileFormat pdf --samplesLabel $labs \
		--heatmapWidth 10 --heatmapHeight 20 \
		--perGroup &

	plotHeatmap -m H3K4me3_${f}_scaled.gz \
	    -o H3K4me3_${f}_scaled.plotHeatmap.2.pdf \
	    --colorMap Greens --whatToShow 'heatmap and colorbar' \
	    --plotFileFormat pdf --samplesLabel $labs \
	    --heatmapWidth 10 --heatmapHeight 20 \
		--regionsLabel HKG.tier1 HKG.tier2 Brain Liver Spleen &
done
wait

################################################ ENCODE ChIP
# awk -F"\t" '$2=="bigWig" && $3=="fold change over control" && $38=="hg19"' metadata.tsv | cut -f 1,7,13,25 > hg19.bigwig.tsv
# awk -F"\t" '$2=="bigWig" && $3=="fold change over control" && $38=="hg19" && $13=="EP300-human"' ../metadata.tsv | cut -f 37 | xargs -n 1 curl -O -L
# awk -F"\t" '$2=="bigWig" && $3=="fold change over control" && $38=="hg19" && $13=="YY1-human"' ../metadata.tsv | cut -f 37 | xargs -n 1 curl -O -L
# awk -F"\t" '$2=="bigWig" && $3=="fold change over control" && $38=="hg19" && $13=="SP1-human"' ../metadata.tsv | cut -f 37 | xargs -n 1 curl -O -L

################################################ HiChIP hichipper loops
# awk '{print $1":"$2"-"$3"\t"$4":"$5"-"$6"\t"$8}' SRR5831489.filt.intra.loop_counts.bedpe > SRR5831489.filt.intra.loop_counts.bedpe.washU.txt
# SRR5831489.intra.loop_counts.bedpe
for pre in HKG v2.hkg.2 Brain Liver Spleen EBV Testis
do
	# bedtools pairtobed -a <(awk -vOFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$8}' SRR5831489.intra.loop_counts.bedpe) -b <(awk -vOFS="\t" '{a=$2-5000;b=$2+5000;if($6=="-"){a=$3-5000;b=$3+5000}if(a<0){a=0;}print $1,a,b,$4,"1000",$6}' $pre.gene.bed) > $pre.hichipper.GM12878_cell_line.SRR5831489.pairtobed.txt
	bedtools pairtobed -a <(awk -vOFS="\t" '$8 >= 5{print $1,$2,$3,$4,$5,$6,$7,$8}' SRR5831489.filt.intra.loop_counts.bedpe) -b <(awk -vOFS="\t" '{a=$2-1000;b=$2+1000;if($6=="-"){a=$3-1000;b=$3+1000}if(a<0){a=0;}print $1,a,b,$4,"1000",$6}' $pre.gene.bed) > $pre.hichipper.GM12878_cell_line.SRR5831489.pairtobed.txt
done
for pre in HKG v2.hkg.2 Brain Liver Spleen EBV Testis
do
	loopped=`cut -f 12 $pre.hichipper.GM12878_cell_line.SRR5831489.pairtobed.txt | sort | uniq | wc -l`
	all=`cat $pre.gene.bed | wc -l`
	perc=`echo $loopped/$all | bc -l`
	echo $pre $all $loopped $perc
done
# HKG 472 430 .91101694915254237288
# v2.hkg.2 2933 2625 .89498806682577565632
# Brain 335 57 .17014925373134328358
# Liver 216 50 .23148148148148148148
# Spleen 246 85 .34552845528455284552
# EBV 228 192 .84210526315789473684
# Testis 1314 399 .30365296803652968036

# awk '{print $1":"$2"-"$3"\t"$4":"$5"-"$6"\t"$8}' SRR5831489.filt.intra.loop_counts.bedpe > SRR5831489.filt.intra.loop_counts.bedpe.washU.txt
for f in SRR5831489.filt.intra.loop_counts.bedpe
do
	awk -vOFS="\t" '$8 >= 5{print $1,$2,$3"\n"$4,$5,$6}' $f | bedtools sort -i - | uniq | awk -vOFS="\t" '{print $1,$2,$3,$1":"$2"-"$3}' | bedtools intersect -wo -a - -b <(awk -vOFS="\t" '{a=$2-1000;b=$2+1000;if($6=="-"){a=$3-1000;b=$3+1000}if(a<0){a=0;}print $1,a,b,$4,"1000",$6}' gencode.v19.gene.bed) > $f.gencode.txt
	perl $myperl <(cut -f 4,8 $f.gencode.txt) <(awk '$3 >= 5' $f.washU.txt) 0 0 | cut -f 1-3,5 | perl $myperl <(cut -f 4,8 $f.gencode.txt) /dev/stdin 0 1 | cut -f 1-4,6 > $f.anntotation.txt
	perl $myperl v2.hkg.tsg.anntotation.txt $f.anntotation.txt 0 3 | cut -f 1-5,7-8 | perl $myperl v2.hkg.tsg.anntotation.txt /dev/stdin 0 4 | cut -f 1-7,9-10 > $f.anntotation.hkg.tsg
	awk -F"\t" -vOFS="\t" '$4!="/" && $5!="/" {print $4,$5,$6,$7,$8,$9}' $f.anntotation.hkg.tsg > $f.pp.txt
	awk -F"\t" -vOFS="\t" '{if($4!="/" && $5=="/") {print $4,$2,$6,$7,$8,$9}}' $f.anntotation.hkg.tsg > $f.pe.txt
	awk -F"\t" -vOFS="\t" '{if($4=="/" && $5!="/") {print $5,$1,$8,$9,$6,$7}}' $f.anntotation.hkg.tsg >> $f.pe.txt
	# cut -f 1,3 $f.pe.txt | awk -F"\t" '$2!="/"' | sort | uniq | cut -f 2 | sort | uniq -c | sort -k 1,1nr
	echo "Class Total promoter-promoter.count promoter-enhancer.count promoter-promoter% promoter-enhancer%"
	for pre in HKG.1 HKG.2 Brain Liver Spleen EBV Testis
	do
		pcount=`awk -vvar=$pre -F"\t" '{if($3==var){print $1}if($4==var){print $2}}' $f.pp.txt | sort | uniq | wc -l`
		ecount=`cut -f 1,3 $f.pe.txt | awk -F"\t" '$2!="/"' | sort | uniq | grep $pre | wc -l`
		all=`grep $pre v2.hkg.tsg.anntotation.txt | wc -l`
		pperc=`echo $pcount/$all | bc -l`
		eperc=`echo $ecount/$all | bc -l`
		echo $pre $all $pcount $ecount $pperc $eperc
	done
	pcount=`awk -F"\t" '{print $1"\n"$2}' $f.pp.txt | sort | uniq | wc -l`
	ecount=`cut -f 1 $f.pe.txt | sort | uniq | wc -l`
	all=`cat gencode.v19.gene.bed | wc -l`
	pperc=`echo $pcount/$all | bc -l`
	eperc=`echo $ecount/$all | bc -l`
	echo gencode $all $pcount $ecount $pperc $eperc

	pre=protein_coding
	pcount=`awk -vvar=$pre -F"\t" '{if($4==var){print $1}if($6==var){print $2}}' $f.pp.txt | sort | uniq | wc -l`
	ecount=`cut -f 1,4 $f.pe.txt | awk -F"\t" '$2!="/"' | sort | uniq | grep $pre | wc -l`
	all=`grep protein_coding gencode.v19.gene.anntotation.txt | wc -l`
	pperc=`echo $pcount/$all | bc -l`
	eperc=`echo $ecount/$all | bc -l`
	echo $pre $all $pcount $ecount $pperc $eperc

	pcount=`awk -F"\t" '{print $1"\n"$2}' $f.pp.txt | sort | uniq | perl $myperl <(awk '$22>2{print $1}' v2.log2tpm.median.tsv) /dev/stdin 0 0 | awk '$2!="/"' | wc -l`
	ecount=`cut -f 1 $f.pe.txt | sort | uniq | perl $myperl <(awk '$22>2{print $1}' v2.log2tpm.median.tsv) /dev/stdin 0 0 | awk '$2!="/"' | wc -l`
	all=`awk '$22>2{print $1}' v2.log2tpm.median.tsv | wc -l`
	pperc=`echo $pcount/$all | bc -l`
	eperc=`echo $ecount/$all | bc -l`
	echo EBV.expressed $all $pcount $ecount $pperc $eperc
done

########## compare loops between expressed and non-expressed genes
for f in SRR5831489.filt.intra.loop_counts.bedpe
do
	mkdir loops.$f
	awk -F"\t" '{print $1"\n"$2}' $f.pp.txt | sort | uniq | perl $myperl <(awk '$22>2{print $1}' v2.log2tpm.median.tsv) /dev/stdin 0 0 | awk '$2!="/"{print $1}' | perl $myperl gencode.v19.gene.bed /dev/stdin 3 0 | cut -f 2- > loops.$f/EBV.expressed.pp.gene.bed
	awk -F"\t" '{print $1"\n"$2}' $f.pp.txt | sort | uniq | perl $myperl <(awk '$22<1{print $1}' v2.log2tpm.median.tsv) /dev/stdin 0 0 | awk '$2!="/"{print $1}' | perl $myperl gencode.v19.gene.bed /dev/stdin 3 0 | cut -f 2- > loops.$f/EBV.nonexpressed.pp.gene.bed
	cut -f 1 $f.pe.txt | sort | uniq | perl $myperl <(awk '$22>2{print $1}' v2.log2tpm.median.tsv) /dev/stdin 0 0 | awk '$2!="/"{print $1}' | perl $myperl gencode.v19.gene.bed /dev/stdin 3 0 | cut -f 2- > loops.$f/EBV.expressed.pe.gene.bed
	cut -f 1 $f.pe.txt | sort | uniq | perl $myperl <(awk '$22<1{print $1}' v2.log2tpm.median.tsv) /dev/stdin 0 0 | awk '$2!="/"{print $1}' | perl $myperl gencode.v19.gene.bed /dev/stdin 3 0 | cut -f 2- > loops.$f/EBV.nonexpressed.pe.gene.bed
	cd loops.$f
	for pre in EBV.expressed.pp EBV.nonexpressed.pp EBV.expressed.pe EBV.nonexpressed.pe
	do
		findMotifsGenome.pl <(awk -vOFS="\t" '{a=$2-1000;b=$2+1000;if($6=="-"){a=$3-1000;b=$3+1000}if(a<0){a=0;}print $1,a,b,$4,"1000",$6}' $pre.gene.bed) hg19 $pre.homer.motifs -size given -mknown /home1/04935/shaojf/myTools/HOMER/data/knownTFs/vertebrates/known.motifs 1> $pre.findMotifsGenome.txt 2>&1 &
		annotatePeaks.pl <(awk -vOFS="\t" '{a=$2-1000;b=$2+1000;if($6=="-"){a=$3-1000;b=$3+1000}if(a<0){a=0;}print $1,a,b,$4,"1000",$6}' $pre.gene.bed) hg19 -m /home1/04935/shaojf/myTools/HOMER/data/knownTFs/vertebrates/known.motifs -size given 1> $pre.promoter.vert.motifs.txt 2> $pre.promoter.vert.motifs.log &
	done
	cd ../
done
########## compare loops between hkg, tsg, and other expressed genes
perl $myperl <(cut -f 1 v2.hkg.tsg.txt) gencode.v19.gene.bed 0 3 | grep "/" | cut -f 1-6 > other.gene.bed
perl $myperl <(awk '$22>2{print $1}' v2.log2tpm.median.tsv) other.gene.bed 0 3 | grep -v "/" | cut -f 1-6 > EBV.expressed.other.gene.bed
for f in SRR5831489.filt.intra.loop_counts.bedpe
do
	cd loops.$f/
	perl $myperl EBV.expressed.pp.gene.bed <(cut -f 4 HKG.gene.bed) 3 0 | grep -v "/" | cut -f 2- > EBV.expressed.HKG1.pp.gene.bed
	perl $myperl EBV.expressed.pe.gene.bed <(cut -f 4 HKG.gene.bed) 3 0 | grep -v "/" | cut -f 2- > EBV.expressed.HKG1.pe.gene.bed
	perl $myperl EBV.expressed.pp.gene.bed <(cut -f 4 v2.hkg.2.gene.bed) 3 0 | grep -v "/" | cut -f 2- > EBV.expressed.HKG2.pp.gene.bed
	perl $myperl EBV.expressed.pe.gene.bed <(cut -f 4 v2.hkg.2.gene.bed) 3 0 | grep -v "/" | cut -f 2- > EBV.expressed.HKG2.pe.gene.bed
	perl $myperl EBV.expressed.pp.gene.bed <(cut -f 4 other.gene.bed) 3 0 | grep -v "/" | cut -f 2- > EBV.expressed.other.pp.gene.bed
	perl $myperl EBV.expressed.pe.gene.bed <(cut -f 4 other.gene.bed) 3 0 | grep -v "/" | cut -f 2- > EBV.expressed.other.pe.gene.bed
	perl $myperl EBV.expressed.pp.gene.bed <(cut -f 4 EBV.gene.bed) 3 0 | grep -v "/" | cut -f 2- > EBV.expressed.tsg.pp.gene.bed
	perl $myperl EBV.expressed.pe.gene.bed <(cut -f 4 EBV.gene.bed) 3 0 | grep -v "/" | cut -f 2- > EBV.expressed.tsg.pe.gene.bed

	perl $myperl <(cat EBV.expressed.pp.gene.bed EBV.expressed.pe.gene.bed | sort | uniq) <(cut -f 4 HKG.gene.bed) 3 0 | grep "/" | cut -f 1 | perl $myperl gencode.v19.gene.bed /dev/stdin 3 0 | cut -f 2- | grep -wv chrM > EBV.expressed.HKG1.noloop.gene.bed
	perl $myperl <(cat EBV.expressed.pp.gene.bed EBV.expressed.pe.gene.bed | sort | uniq) <(cut -f 4 v2.hkg.2.gene.bed) 3 0 | grep "/" | cut -f 1 | perl $myperl gencode.v19.gene.bed /dev/stdin 3 0 | cut -f 2- | grep -wv chrM > EBV.expressed.HKG2.noloop.gene.bed
	perl $myperl <(cat EBV.expressed.pp.gene.bed EBV.expressed.pe.gene.bed | sort | uniq) <(cut -f 4 EBV.expressed.other.gene.bed) 3 0 | grep "/" | cut -f 1 | perl $myperl gencode.v19.gene.bed /dev/stdin 3 0 | cut -f 2- | grep -wv chrM > EBV.expressed.other.noloop.gene.bed
	perl $myperl <(cat EBV.expressed.pp.gene.bed EBV.expressed.pe.gene.bed | sort | uniq) <(cut -f 4 EBV.gene.bed) 3 0 | grep "/" | cut -f 1 | perl $myperl gencode.v19.gene.bed /dev/stdin 3 0 | cut -f 2- | grep -wv chrM > EBV.expressed.tsg.noloop.gene.bed

	for g in EBV.expressed.*.*.gene.bed
	do
		pre=`echo $g | sed 's/.gene.bed//'`
		findMotifsGenome.pl <(awk -vOFS="\t" '{a=$2-1000;b=$2+1000;if($6=="-"){a=$3-1000;b=$3+1000}if(a<0){a=0;}print $1,a,b,$4,"1000",$6}' $pre.gene.bed) hg19 $pre.homer.motifs -size given -mknown /home1/04935/shaojf/myTools/HOMER/data/knownTFs/vertebrates/known.motifs 1> $pre.findMotifsGenome.txt 2>&1 &
		annotatePeaks.pl <(awk -vOFS="\t" '{a=$2-1000;b=$2+1000;if($6=="-"){a=$3-1000;b=$3+1000}if(a<0){a=0;}print $1,a,b,$4,"1000",$6}' $pre.gene.bed) hg19 -m /home1/04935/shaojf/myTools/HOMER/data/knownTFs/vertebrates/known.motifs -size given 1> $pre.promoter.vert.motifs.txt 2> $pre.promoter.vert.motifs.log &
	done
	cd ../
done

for pre in `ls *promoter.vert.motifs.txt | sed 's/.promoter.vert.motifs.txt//'`
do
	ets=`cut -f 1,103 $pre.promoter.vert.motifs.txt | tail -n +2 | awk -F"\t" '$2~/,/' | wc -l`
	all=`tail -n +2 $pre.promoter.vert.motifs.txt | wc -l`
	perc=`echo $ets/$all | bc -l`
	echo $pre $ets $all $perc
done
for pre in `ls *promoter.vert.motifs.txt | sed 's/.promoter.vert.motifs.txt//'`
do
	yy1=`cut -f 1,354 $pre.promoter.vert.motifs.txt | tail -n +2 | awk -F"\t" '$2~/,/' | wc -l`
	all=`tail -n +2 $pre.promoter.vert.motifs.txt | wc -l`
	perc=`echo $yy1/$all | bc -l`
	echo $pre $yy1 $all $perc
done

echo "GeneType Motif MotifCount GeneCount GeneWithMotif" | tr " " "\t" > homer.vert.motif.percentage.txt
for pre in `ls *promoter.vert.motifs.txt | sed 's/.promoter.vert.motifs.txt//'`
do
	all=`tail -n +2 $pre.promoter.vert.motifs.txt | wc -l`
	for i in `seq 22 384`
	do
		motif=`head_line $pre.promoter.vert.motifs.txt | awk -vvar=$i '$1==var{print $2}' | awk -F"/" '{print $1"/"$2}'`
		mcount=`cut -f 1,$i $pre.promoter.vert.motifs.txt | tail -n +2 | awk -F"\t" '$2~/,/' | wc -l`
		perc=`echo $mcount/$all | bc -l`
		echo $pre $motif $mcount $all $perc | tr " " "\t" >> homer.vert.motif.percentage.txt
	done
done
perl ~/myTools/BioinformaticsDaily/textProcess/make_matrixwith3col_from_single_file.pl <(awk -F"\t" -vOFS="\t" '{print $2,$1,$5}' homer.vert.motif.percentage.txt | tail -n +2) > homer.vert.motif.percentage.mat
# ####################################################
# args <- c("homer.vert.motif.percentage.mat")
# library(pheatmap)
# library(data.table)
# input <- fread(args[1], sep = "\t", header = T, na.strings = "/")
# scores <- data.matrix(input[,-1])
# colors <- colorRampPalette(c("white", "red"))(10)
# rownames(scores) <- as.matrix(input[,1])[,1]

# pdf(file = paste(args[1], "pheatmap.pdf", sep = "."), width = 10, height = 16)
# myplot <- pheatmap(scores, scale = "none", fontsize_row = 4, 
# 	show_rownames = T, show_colnames = T, color = colors,
# 	cluster_cols = T, cluster_rows = T)
# dev.off()
# diff.scores <- scores[(scores[,3] > 0.2 | scores[,7] > 0.2) & (scores[,3]/(scores[,7]+1e-10) > 1.5 | scores[,3]/(scores[,7]+1e-10) < 0.67), ]
# pdf(file = paste(args[1], "diff.pheatmap.pdf", sep = "."), width = 10, height = 8)
# myplot <- pheatmap(diff.scores, scale = "none", fontsize_row = 8, 
# 	show_rownames = T, show_colnames = T, color = colors,
# 	cluster_cols = T, cluster_rows = T)
# dev.off()
# score.min <- apply(scores, 1, min)
# score.max <- apply(scores, 1, max)
# diff.scores.1 <- scores[score.max / (score.min + 1e-10) > 2 & score.max > 0.2, ]
# pdf(file = paste(args[1], "diff.1.pheatmap.pdf", sep = "."), width = 10, height = 8)
# myplot <- pheatmap(diff.scores.1, scale = "none", fontsize_row = 8, 
# 	show_rownames = T, show_colnames = T, color = colors,
# 	cluster_cols = T, cluster_rows = T)
# dev.off()
# ####################################################

