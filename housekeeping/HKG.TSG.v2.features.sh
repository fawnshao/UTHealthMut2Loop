#!/bin/sh
# rclone sync /home1/04935/shaojf/stampede2/housekeeping_genes/v2.feature/ mygoogle:hkg_tsg/v2.feature/
# args <- c("v2", "GTEx_sample.tissue.txt", "0.15", "1.5", "2", "0.25")
# outputpre <- args[1]
# tau.threshold <- as.numeric(args[3])
# sd.threshold <- as.numeric(args[4])
# median.threshold <- as.numeric(args[5])
# tau.threshold2 <- as.numeric(args[6])

myperl=/home1/04935/shaojf/myTools/BioinformaticsDaily/textProcess/add_any_2files_together.pl

# cd /home1/04935/shaojf/stampede2/housekeeping_genes/v2.feature
echo "Gene Type" | sed 's/ /\t/' > v2.hkg.tsg.txt
cut -f 1-2 v2.HKG.tsv | tail -n +2 | sort -k2,2nr | cut -f 1 | awk '{print $0"\tHKG.1"}' >> v2.hkg.tsg.txt
cut -f 1-2 v2.HKG.2.tsv | tail -n +2 | sort -k2,2nr | cut -f 1 | awk '{print $0"\tHKG.2"}' >> v2.hkg.tsg.txt
cut -f 1,100 v2.TSG.tsv | tail -n +2 | grep -v "," | sort -k2 >> v2.hkg.tsg.txt
cut -f 1,100 v2.TSG.tsv | tail -n +2 | grep "," | sort -k2 >> v2.hkg.tsg.txt

perl $myperl <(sed 's/\t/|/' gencode.v19.gene.anntotation.txt) v2.hkg.tsg.txt 0 0 | cut -f 1-2,4 > v2.hkg.tsg.anntotation.txt
perl $myperl sim.gencode.vert.known.motifs.mat v2.hkg.tsg.anntotation.txt 0 0 | cut -f 1-3,5- > v2.hkg.tsg.vert.known.motifs.mat
# Rscript ~/myTools/UTHealthMut2Loop/housekeeping/heatmap.for.feature.selection.R v2.hkg.tsg.vert.known.motifs.mat
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


awk -F"\t" '$2=="bigWig" && $3=="fold change over control" && $38=="hg19"' metadata.tsv | cut -f 1,7,13,25 > hg19.bigwig.tsv
awk -F"\t" '$2=="bigWig" && $3=="fold change over control" && $38=="hg19" && $13=="EP300-human"' ../metadata.tsv | cut -f 37 | xargs -n 1 curl -O -L
awk -F"\t" '$2=="bigWig" && $3=="fold change over control" && $38=="hg19" && $13=="YY1-human"' ../metadata.tsv | cut -f 37 | xargs -n 1 curl -O -L
awk -F"\t" '$2=="bigWig" && $3=="fold change over control" && $38=="hg19" && $13=="SP1-human"' ../metadata.tsv | cut -f 37 | xargs -n 1 curl -O -L