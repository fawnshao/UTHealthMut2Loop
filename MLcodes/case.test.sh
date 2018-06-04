#!/bin/sh
myperl=/home1/04935/shaojf/myScripts/add_any_2files_together.pl
# rclone sync /home1/04935/shaojf/stampede2/housekeeping_genes/both.pc.and.nc.genes/GM12878/ mygoogle:hkg_tsg/both.pc.and.nc.genes/GM12878/

bedtools intersect -wo -a hg19.cpgIslandExt.withCount.bed -b <(awk -vOFS="\t" '{a=$2-300;b=$2+100;if($6=="-"){a=$3-100;b=$3+300}if(a<0){a=0;}print $1,a,b,$4,$5,$6}' gencode.v19.transcript.bed) > gencode.cpgIslandExt.txt

gene=hkg.tsg.srtbyPCA.transcript.bed
pre=hkg.tsg.srtbyPCA.transcript
######## whole gene body plus +/- 5k #######
for post in GM12878 K562 Naive Th17 Treg
do
	postfile=Mumbach_NG2017_HiChIP_GSE101498/${post}_H3K27ac_Loops.txt
	bedtools intersect -wao -a <(cut -f 3-  $gene | awk -vOFS="\t" '{a=$2-5000;b=$3+5000;if($6=="-"){a=$2-5000;b=$3+5000}if(a<0){a=0;}print $1,a,b,$4,$5,$6}') -b <(tail -n +2 $postfile | awk '{print "chr"$1"\t"$2"\t"$3"\t"$4":"$5"-"$6"\nchr"$4"\t"$5"\t"$6"\t"$1":"$2"-"$3}') > $pre.HiChIP-Loop.$post.overlap.txt
done

for pre in hkg.tsg.srtbyPCA.transcript
do
	cut -f 4 $pre.HiChIP-Loop.GM12878.overlap.txt | sort | uniq > a
	cmd="paste a"
	for post in GM12878 K562 Naive Th17 Treg
	do
		awk '$7!="."{print $4}' $pre.HiChIP-Loop.$post.overlap.txt | sort | uniq -c | awk '{print $2"\t"$1}' | perl $myperl /dev/stdin a 0 0 | cut -f 3 > a.$post
		cmd=$cmd" a.$post"
	done
	echo "gene GM12878 K562 Naive Th17 Treg" | tr ' ' "\t" > $pre.HiChIP-Loop.count
	$cmd >> $pre.HiChIP-Loop.count
	rm a a.*
done

perl $myperl $pre.HiChIP-Loop.count $pre.class 0 2 | cut -f 1-3,5- > $pre.HiChIP-Loop.count.mat
cut -f 1-2,4- $pre.HiChIP-Loop.count.mat | uniq > $pre.HiChIP-Loop.count.sim.mat

head -1 $pre.HiChIP-Loop.count.sim.mat > hkg.gm.liver.HiChIP-Loop.count.sim.mat
awk -F"\t" '$2~/hkg/' $pre.HiChIP-Loop.count.sim.mat >> hkg.gm.liver.HiChIP-Loop.count.sim.mat
awk -F"\t" '$2=="Cells - EBV-transformed lymphocytes"' $pre.HiChIP-Loop.count.sim.mat >> hkg.gm.liver.HiChIP-Loop.count.sim.mat
awk -F"\t" '$2=="Liver"' $pre.HiChIP-Loop.count.sim.mat >> hkg.gm.liver.HiChIP-Loop.count.sim.mat

grep -v "/" hkg.gm.liver.HiChIP-Loop.count.mat > hkg.gm.liver.withHiChIP-Loop.txt
head -1 hkg.gm.liver.HiChIP-Loop.count.mat > hkg.gm.liver.withoutHiChIP-Loop.txt
awk -F"\t" '$3=="/" && $4=="/" && $5=="/" && $6=="/" && $7=="/" ' hkg.gm.liver.HiChIP-Loop.count.mat >> hkg.gm.liver.withoutHiChIP-Loop.txt

perl $myperl $pre.HiChIP-Loop.count $pre.class 0 2 | cut -f 2-3,5- > appris.tss.$pre.HiChIP-Loop.count.mat
head -1 appris.tss.$pre.HiChIP-Loop.count.mat > appris.tss.hkg.gm.liver.HiChIP-Loop.count.mat
awk -F"\t" '$1~/hkg/' appris.tss.$pre.HiChIP-Loop.count.mat >> appris.tss.hkg.gm.liver.HiChIP-Loop.count.mat
awk -F"\t" '$1=="Cells - EBV-transformed lymphocytes"' appris.tss.$pre.HiChIP-Loop.count.mat >> appris.tss.hkg.gm.liver.HiChIP-Loop.count.mat
awk -F"\t" '$1=="Liver"' appris.tss.$pre.HiChIP-Loop.count.mat >> appris.tss.hkg.gm.liver.HiChIP-Loop.count.mat

grep -v "/" appris.tss.hkg.gm.liver.HiChIP-Loop.count.mat > appris.tss.hkg.gm.liver.withHiChIP-Loop.txt
head -1 appris.tss.hkg.gm.liver.HiChIP-Loop.count.mat > appris.tss.hkg.gm.liver.withoutHiChIP-Loop.txt
awk -F"\t" '$3=="/" && $4=="/" && $5=="/" && $6=="/" && $7=="/" ' appris.tss.hkg.gm.liver.HiChIP-Loop.count.mat >> appris.tss.hkg.gm.liver.withoutHiChIP-Loop.txt


#######################
# homer
perl $myperl hkg.tsg.srtbyPCA.transcript.bed appris.tss.hkg.gm.liver.withHiChIP-Loop.txt 5 1 | cut -f 10- | tail -n +2 > appris.tss.hkg.gm.liver.withHiChIP-Loop.bed
perl $myperl hkg.tsg.srtbyPCA.transcript.bed appris.tss.hkg.gm.liver.withoutHiChIP-Loop.txt 5 1 | cut -f 10- | tail -n +2 > appris.tss.hkg.gm.liver.withoutHiChIP-Loop.bed

perl $myperl hkg.tsg.srtbyPCA.transcript.bed <(awk '$1~/^hkg/' appris.tss.hkg.gm.liver.withHiChIP-Loop.txt) 5 1 | cut -f 10- > hkg.withHiChIP-Loop.bed
perl $myperl hkg.tsg.srtbyPCA.transcript.bed <(awk '$1~/^hkg/' appris.tss.hkg.gm.liver.withoutHiChIP-Loop.txt) 5 1 | cut -f 10- > hkg.withoutHiChIP-Loop.bed
perl $myperl hkg.tsg.srtbyPCA.transcript.bed <(awk '$0~/\// && $1~/^hkg/ && ($3!="/" || $4!="/" || $5!="/" || $6!="/" || $7!="/")' appris.tss.hkg.gm.liver.HiChIP-Loop.count.mat) 5 1 | cut -f 10- > hkg.someHiChIP-Loop.bed

perl $myperl hkg.tsg.srtbyPCA.transcript.bed <(awk -F"\t" '$1~/EBV/ && $3!="/"' appris.tss.hkg.gm.liver.HiChIP-Loop.count.mat) 5 1 | cut -f 10- > gm.withHiChIP-Loop.bed
perl $myperl hkg.tsg.srtbyPCA.transcript.bed <(awk -F"\t" '$1~/EBV/ && $3=="/"' appris.tss.hkg.gm.liver.HiChIP-Loop.count.mat) 5 1 | cut -f 10- > gm.withoutHiChIP-Loop.bed

perl $myperl hkg.tsg.srtbyPCA.transcript.bed <(awk -F"\t" '$1~/Liver/ && $3!="/"' appris.tss.hkg.gm.liver.HiChIP-Loop.count.mat) 5 1 | cut -f 10- > liver.withHiChIP-Loop.bed
perl $myperl hkg.tsg.srtbyPCA.transcript.bed <(awk -F"\t" '$1~/Liver/ && $3=="/"' appris.tss.hkg.gm.liver.HiChIP-Loop.count.mat) 5 1 | cut -f 10- > liver.withoutHiChIP-Loop.bed

for gene in *with*HiChIP-Loop.bed
do
	pre=`echo $gene | sed 's/.bed//'`
	findMotifsGenome.pl <(awk -vOFS="\t" '{a=$2-1000;b=$2+1000;if($6=="-"){a=$3-1000;b=$3+1000}if(a<0){a=0;}print $1,a,b,$4,"1000",$6}' $gene) hg19 $pre.homer.motifs -size given -mknown /home1/04935/shaojf/myTools/HOMER/data/knownTFs/all.motifs 1> $pre.findMotifsGenome.txt 2>&1 &
	annotatePeaks.pl <(awk -vOFS="\t" '{a=$2-1000;b=$2+1000;if($6=="-"){a=$3-1000;b=$3+1000}if(a<0){a=0;}print $1,a,b,$4,"1000",$6}' $gene) hg19 -m /home1/04935/shaojf/myTools/HOMER/data/knownTFs/vertebrates/known.motifs -size given -cpu 60 1> $pre.promoter.vert.motifs.txt 2> $pre.promoter.vert.motifs.log &
done
gene=hkg.someHiChIP-Loop.bed
pre=`echo $gene | sed 's/.bed//'`
findMotifsGenome.pl <(awk -vOFS="\t" '{a=$2-1000;b=$2+1000;if($6=="-"){a=$3-1000;b=$3+1000}if(a<0){a=0;}print $1,a,b,$4,"1000",$6}' $gene) hg19 $pre.homer.motifs -size given -mknown /home1/04935/shaojf/myTools/HOMER/data/knownTFs/all.motifs -p 68 1> $pre.findMotifsGenome.txt 2>&1 &
annotatePeaks.pl <(awk -vOFS="\t" '{a=$2-1000;b=$2+1000;if($6=="-"){a=$3-1000;b=$3+1000}if(a<0){a=0;}print $1,a,b,$4,"1000",$6}' $gene) hg19 -m /home1/04935/shaojf/myTools/HOMER/data/knownTFs/vertebrates/known.motifs -size given -cpu 60 1> $pre.promoter.vert.motifs.txt 2> $pre.promoter.vert.motifs.log &

#####################################################################
###### find the pattern for hkg enhancers
pre=hkg.tsg.srtbyPCA.transcript
for post in GM12878 K562 Naive Th17 Treg
do
	# awk -F"\t" '$7!="."{print $4"\t"$10}' $pre.HiChIP-Loop.$post.overlap.txt | sed 's/|/\t/' | cut -f1,3 | sort | uniq > $pre.$post.enhancers
	awk -F"\t" '$7!="."{print $4"\t"$10}' $pre.HiChIP-Loop.$post.overlap.txt | sort | uniq > $pre.$post.enhancers
done
echo "gene%enhancer cell" | tr " " "\t" > $pre.enhancers.stats
for post in GM12878 K562 Naive Th17 Treg
do
	awk -F"\t" '$7!="."{print $4"\t"$10}' $pre.HiChIP-Loop.$post.overlap.txt | sort | uniq | awk -vvar=$post '{print $1"%"$2"\t"var}' >> $pre.enhancers.stats
done
perl /home1/04935/shaojf/myTools/BioinformaticsDaily/textProcess/make_binarymatrix_from_single_file.pl <(tail -n +2 $pre.enhancers.stats) > $pre.enhancers.mat
paste <(head -1 appris.tss.hkg.gm.liver.withHiChIP-Loop.txt | cut -f 1-2) <(head -1 $pre.enhancers.mat) > appris.tss.hkg.gm.liver.withHiChIP-Loop.enhancers.mat
perl $myperl <(tail -n +2 $pre.enhancers.mat | sed 's/%/\t/') <(cut -f 1-2 appris.tss.hkg.gm.liver.withHiChIP-Loop.txt | tail -n +2) 0 1 | cut -f 1-2,4- >> appris.tss.hkg.gm.liver.withHiChIP-Loop.enhancers.mat
awk -F"\t" '$4==1 && $5==1 && $6==1 && $7==1 && $8==1' appris.tss.hkg.gm.liver.withHiChIP-Loop.enhancers.mat > appris.tss.hkg.gm.liver.withHiChIP-Loop.enhancers.commonenhancers.mat
awk -F"\t" '$4==1 && $5==1 && $6==1 && $7==1 && $8==1 {print $3"\t"$2"%"$3}' appris.tss.hkg.gm.liver.withHiChIP-Loop.enhancers.mat | sed 's/:/\t/;s/-/\t/' > hkg.commonenhancers.bed
awk -F"\t" '$1~/hkg/ && $4+$5+$6+$7+$8 > 1 && $4+$5+$6+$7+$8 < 5 {print $3"\t"$2"%"$3}' appris.tss.hkg.gm.liver.withHiChIP-Loop.enhancers.mat | sed 's/:/\t/;s/-/\t/' > hkg.partialcommonenhancers.bed
for f in hkg.*commonenhancers.bed
do
	pre=`echo $f | sed 's/.bed//'`
	findMotifsGenome.pl <(cut -f 1-3 $f | sort | uniq | awk -vOFS="\t" '{print "chr"$1,$2,$3,NR}') hg19 $pre.homer.motifs -size given -mknown /home1/04935/shaojf/myTools/HOMER/data/knownTFs/vertebrates/known.motifs -p 60 1> $pre.findMotifsGenome.txt 2>&1 &
	annotatePeaks.pl <(cut -f 1-3 $f | sort | uniq | awk -vOFS="\t" '{print "chr"$1,$2,$3,NR}') hg19 -m /home1/04935/shaojf/myTools/HOMER/data/knownTFs/vertebrates/known.motifs -size given -cpu 60 1> $pre.vert.motifs.txt 2> $pre.vert.motifs.log &
done
################################################################################################
# library(pheatmap)
# library(data.table)
# args <- c("appris.tss.hkg.gm.liver.withHiChIP-Loop.enhancers.mat")
# input <- fread(args[1], sep = "\t", header = T, na.strings = "/")
# datax <- data.matrix(input[,4:8])
# rownames(datax) <- paste(as.matrix(input[,2]), as.matrix(input[,3]), sep = "%")
# colors <- colorRampPalette(c("white", "blue"))(100)
# annos <- input[,1]
# annos[Type=="Cells - EBV-transformed lymphocytes", 1] <- "EBV"
# annos[Type=="hkg1", 1] <- "HKG"
# annos[Type=="hkg2", 1] <- "HKG"
# annos[Type=="hkg3", 1] <- "HKG"
# annos[Type=="hkg4", 1] <- "HKG"
# rownames(annos) <- rownames(datax)
# pdf(file = paste(args[1], "pheatmap.pdf", sep = "."), width = 6, height = 6)
# # datax[order(rowSums(datax), decreasing = T),]
# myplot <- pheatmap(datax[order(datax[,1],datax[,2],datax[,3],datax[,4],datax[,5],decreasing = T),], scale = "none", annotation_row = annos,
# 	show_rownames = F, show_colnames = T, color = colors, fontsize_row = 3,
# 	cluster_cols = F, cluster_rows = F)
# dev.off()
################################################################################################

###for promoters with common enhancers
myperl=/home1/04935/shaojf/myScripts/add_any_2files_together.pl
for f in hkg.*commonenhancers.bed
do
	pre=`echo $f | sed 's/.bed//'`
	findMotifsGenome.pl <(sed 's/%/\t/' $f | cut -f 4 | perl $myperl hkg.tsg.srtbyPCA.transcript.bed /dev/stdin 5 0 | cut -f 4- | awk -vOFS="\t" '{a=$2-1000;b=$2+1000;if($6=="-"){a=$3-1000;b=$3+1000}if(a<0){a=0;}print $1,a,b,$4,"1000",$6}' | sort | uniq) hg19 $pre.promoter.homer.motifs -size given -mknown /home1/04935/shaojf/myTools/HOMER/data/knownTFs/vertebrates/known.motifs -p 60 1> $pre.promoter.findMotifsGenome.txt 2>&1 &
	annotatePeaks.pl <(sed 's/%/\t/' $f | cut -f 4 | perl $myperl hkg.tsg.srtbyPCA.transcript.bed /dev/stdin 5 0 | cut -f 4- | awk -vOFS="\t" '{a=$2-1000;b=$2+1000;if($6=="-"){a=$3-1000;b=$3+1000}if(a<0){a=0;}print $1,a,b,$4,"1000",$6}' | sort | uniq) hg19 -m /home1/04935/shaojf/myTools/HOMER/data/knownTFs/vertebrates/known.motifs -size given -cpu 60 1> $pre.promoter.vert.motifs.txt 2> $pre.promoter.vert.motifs.log &
done
##############for gm12878
perl $myperl hkg.tsg.srtbyPCA.transcript.GM12878.enhancers <(awk -F"\t" '$2=="Cells - EBV-transformed lymphocytes"' hkg.tsg.srtbyPCA.transcript.bed | cut -f 3-) 0 3 | awk -F"\t" '$7!="/"' | cut -f 1-6 | sort | uniq > gm.withenhancer.transcript.bed
perl $myperl hkg.tsg.srtbyPCA.transcript.GM12878.enhancers <(awk -F"\t" '$2=="Cells - EBV-transformed lymphocytes"' hkg.tsg.srtbyPCA.transcript.bed | cut -f 3-) 0 3 | awk -F"\t" '$7=="/"' | cut -f 1-6  | sort | uniq > gm.withoutenhancer.transcript.bed
perl $myperl hkg.tsg.srtbyPCA.transcript.GM12878.enhancers <(awk -F"\t" '$2=="Cells - EBV-transformed lymphocytes"' hkg.tsg.srtbyPCA.transcript.bed | cut -f 3-) 0 3 | awk -F"\t" '$7!="/"{print $8}' | sort | uniq | sed 's/:/\t/;s/-/\t/' | awk '{print $0"\t"NR}' > gm.withenhancer.enhancer.bed
perl $myperl hkg.tsg.srtbyPCA.transcript.GM12878.enhancers <(awk -F"\t" '$2=="Liver"' hkg.tsg.srtbyPCA.transcript.bed | cut -f 3-) 0 3 | awk -F"\t" '$7!="/"' | cut -f 1-6 | sort | uniq > liver.withenhancer.transcript.bed
perl $myperl hkg.tsg.srtbyPCA.transcript.GM12878.enhancers <(awk -F"\t" '$2=="Liver"' hkg.tsg.srtbyPCA.transcript.bed | cut -f 3-) 0 3 | awk -F"\t" '$7=="/"' | cut -f 1-6  | sort | uniq > liver.withoutenhancer.transcript.bed
perl $myperl hkg.tsg.srtbyPCA.transcript.GM12878.enhancers <(awk -F"\t" '$2=="Liver"' hkg.tsg.srtbyPCA.transcript.bed | cut -f 3-) 0 3 | awk -F"\t" '$7!="/"{print $8}' | sort | uniq | sed 's/:/\t/;s/-/\t/' | awk '{print $0"\t"NR}' > liver.withenhancer.enhancer.bed

for f in gm.with*enhancer.transcript.bed liver.with*enhancer.transcript.bed
do
	pre=`echo $f | sed 's/.bed//'`
	findMotifsGenome.pl <(awk -vOFS="\t" '{a=$2-1000;b=$2+1000;if($6=="-"){a=$3-1000;b=$3+1000}if(a<0){a=0;}print $1,a,b,$4,"1000",$6}' $f) hg19 $pre.promoter.homer.motifs -size given -mknown /home1/04935/shaojf/myTools/HOMER/data/knownTFs/vertebrates/known.motifs -p 60 1> $pre.promoter.findMotifsGenome.txt 2>&1 &
	annotatePeaks.pl <(awk -vOFS="\t" '{a=$2-1000;b=$2+1000;if($6=="-"){a=$3-1000;b=$3+1000}if(a<0){a=0;}print $1,a,b,$4,"1000",$6}' $f) hg19 -m /home1/04935/shaojf/myTools/HOMER/data/knownTFs/vertebrates/known.motifs -size given -cpu 60 1> $pre.promoter.vert.motifs.txt 2> $pre.promoter.vert.motifs.log &
done
### functional enhancer? maybe we should use genes not expressed but with enhancer looping
### most liver specific genes are not expressed in GM cells
for f in gm.withenhancer.enhancer.bed liver.withenhancer.enhancer.bed
do
	pre=`echo $f | sed 's/.bed//'`
	findMotifsGenome.pl $f hg19 $pre.homer.motifs -size given -mknown /home1/04935/shaojf/myTools/HOMER/data/knownTFs/vertebrates/known.motifs -p 68 1> $pre.findMotifsGenome.txt 2>&1 &
	annotatePeaks.pl $f hg19 -m /home1/04935/shaojf/myTools/HOMER/data/knownTFs/vertebrates/known.motifs -size given -cpu 68 1> $pre.vert.motifs.txt 2> $pre.promoter.vert.motifs.log &
done
### test all enhancers linked to promoters
# cd looped.ep.motifs/
myperl=/home1/04935/shaojf/myScripts/add_any_2files_together.pl
tail -n +2 hkg.tsg.srtbyPCA.transcript.enhancers.mat | cut -f 1 | cut -f 2 -d"%" | sort | uniq | sed 's/:/\t/;s/-/\t/' | awk '{print "chr"$0"\t"NR}' > all.looped.enhancers.bed
tail -n +2 hkg.tsg.srtbyPCA.transcript.enhancers.mat | cut -f 1 | cut -f 1 -d"%" | sort | uniq | perl $myperl hkg.tsg.srtbyPCA.transcript.bed /dev/stdin 5 0 | cut -f 4- > all.looped.promoters.bed
tail -n +2 hkg.tsg.srtbyPCA.transcript.enhancers.mat | cut -f 1 | cut -f 1 -d"%" | sort | uniq | perl $myperl /dev/stdin hkg.tsg.srtbyPCA.transcript.bed 0 5 | grep "/" | cut -f 3-9 > all.notlooped.promoters.bed
for f in all.looped.enhancers.bed
do
	pre=`echo $f | sed 's/.bed//'`
	findMotifsGenome.pl $f hg19 $pre.homer.motifs -size given -mknown /home1/04935/shaojf/myTools/HOMER/data/knownTFs/vertebrates/known.motifs -p 68 1> $pre.findMotifsGenome.txt 2>&1 &
	annotatePeaks.pl $f hg19 -m /home1/04935/shaojf/myTools/HOMER/data/knownTFs/vertebrates/known.motifs -size given -cpu 68 1> $pre.vert.motifs.txt 2> $pre.promoter.vert.motifs.log &
done
for f in all.*looped.promoters.bed
do
	pre=`echo $f | sed 's/.bed//'`
	findMotifsGenome.pl <(awk -vOFS="\t" '{a=$2-1000;b=$2+1000;if($6=="-"){a=$3-1000;b=$3+1000}if(a<0){a=0;}print $1,a,b,$4,"1000",$6}' $f) hg19 $pre.homer.motifs -size given -mknown /home1/04935/shaojf/myTools/HOMER/data/knownTFs/vertebrates/known.motifs -p 68 1> $pre.findMotifsGenome.txt 2>&1 &
	annotatePeaks.pl <(awk -vOFS="\t" '{a=$2-1000;b=$2+1000;if($6=="-"){a=$3-1000;b=$3+1000}if(a<0){a=0;}print $1,a,b,$4,"1000",$6}' $f) hg19 -m /home1/04935/shaojf/myTools/HOMER/data/knownTFs/vertebrates/known.motifs -size given -cpu 68 1> $pre.vert.motifs.txt 2> $pre.vert.motifs.log &
done

perl $myperl <(tail -n +2 hkg.tsg.srtbyPCA.transcript.enhancers.mat | cut -f 1 | sed 's/%/\t/') <(awk -F"\t" '$2~/hkg/' hkg.tsg.srtbyPCA.transcript.bed) 0 5 | grep -v "/" | cut -f 10 | sort | uniq | sed 's/:/\t/;s/-/\t/' | awk '{print "chr"$0"\t"NR}' > all.looped2hkg.enhancers.bed
perl $myperl <(tail -n +2 hkg.tsg.srtbyPCA.transcript.enhancers.mat | cut -f 1 | sed 's/%/\t/') <(awk -F"\t" '$2!~/hkg/' hkg.tsg.srtbyPCA.transcript.bed) 0 5 | grep -v "/" | cut -f 10 | sort | uniq | sed 's/:/\t/;s/-/\t/' | awk '{print "chr"$0"\t"NR}' > all.looped2tsg.enhancers.bed
for f in all.looped2??g.enhancers.bed
do
	pre=`echo $f | sed 's/.bed//'`
	findMotifsGenome.pl $f hg19 $pre.homer.motifs -size given -mknown /home1/04935/shaojf/myTools/HOMER/data/knownTFs/vertebrates/known.motifs -p 68 1> $pre.findMotifsGenome.txt 2>&1 &
	annotatePeaks.pl $f hg19 -m /home1/04935/shaojf/myTools/HOMER/data/knownTFs/vertebrates/known.motifs -size given -cpu 68 1> $pre.vert.motifs.txt 2> $pre.vert.motifs.log &
done

#######combine all the known vert motifs
# cd motifs.compare/
# rclone sync /home1/04935/shaojf/stampede2/housekeeping_genes/both.pc.and.nc.genes/GM12878/motifs.compare/ mygoogle:hkg_tsg/both.pc.and.nc.genes/GM12878/motifs.compare/
grep -e "/Promoter/" *.homer.motifs/knownResults.txt | cut -f1,7
grep -e "CG-repeat" *.homer.motifs/knownResults.txt | cut -f1,7
for f in *.promoter*.homer.motifs/knownResults.txt hkg.withoutHiChIP-Loop.homer.motifs/knownResults.txt
do
	pre=`echo $f | cut -f 1-2 -d"."`
	awk -F"\t" -vvar=$pre -vOFS="\t" '{print var,$1,$5,$7}' $f >> promoter.motifs
done

grep "/Promoter/" promoter.motifs | cut -f 1-2,4 | sed 's/%//' > promoter.Promoter.motifs
Rscript /home1/04935/shaojf/myTools/BioinformaticsDaily/RVisualization/barplot4threecolsdata.R promoter.Promoter.motifs
grep "/Promoter/" promoter.motifs | grep -e "CRE" -e "ETS" -e "Sp1" -e "YY1" | cut -f 1-2,4 | sed 's/%//' > promoter.simPromoter.motifs
Rscript /home1/04935/shaojf/myTools/BioinformaticsDaily/RVisualization/barplot4threecolsdata.R promoter.simPromoter.motifs

# sed 's/%//' enhancer.motifs | awk '$4>50 && $3<0.01'
grep -e "ELF1" -e "Elk4" -e "Elk1" -e "ETV1" -e "ETS(" promoter.motifs | cut -f 1-2,4 | sed 's/%//' | grep -v "Jaspar" | grep -v "withenhancer" > promoter.ETS.motifs
grep -e "NRF" -e "Klf9" -e "ZNF143" -e "E2A(bHLH),near_PU.1" promoter.motifs | cut -f 1-2,4 | sed 's/%//' | grep -v "Jaspar" | grep -v "withenhancer" > promoter.others.motifs
grep -e "YY1" -e "CTCF(" -e "TATA-Box" -e "Sp1(" promoter.motifs | cut -f 1-2,4 | sed 's/%//' | grep -v "Jaspar" | grep -v "withenhancer" > promoter.others.1.motifs
Rscript /home1/04935/shaojf/myTools/BioinformaticsDaily/RVisualization/barplot4threecolsdata.R promoter.ETS.motifs
Rscript /home1/04935/shaojf/myTools/BioinformaticsDaily/RVisualization/barplot4threecolsdata.R promoter.others.motifs
Rscript /home1/04935/shaojf/myTools/BioinformaticsDaily/RVisualization/barplot4threecolsdata.R promoter.others.1.motifs


for f in *enhancer*.homer.motifs/knownResults.txt
do
	pre=`echo $f | cut -f 1-2 -d"."`
	awk -F"\t" -vvar=$pre -vOFS="\t" '{print var,$1,$5,$7}' $f >> enhancer.motifs
done

grep "/Promoter/" enhancer.motifs | cut -f 1-2,4 | sed 's/%//' > enhancer.Promoter.motifs
Rscript /home1/04935/shaojf/myTools/BioinformaticsDaily/RVisualization/barplot4threecolsdata.R enhancer.Promoter.motifs

grep -e "ELF1" -e "Elk4" -e "Elk1" -e "ETV1" -e "ETS(" enhancer.motifs | cut -f 1-2,4 | sed 's/%//' | grep -v "Jaspar" | grep -v "withenhancer"  > enhancer.ETS.motifs
grep -e "NRF" -e "Klf9" -e "ZNF143" -e "E2A(bHLH),near_PU.1" enhancer.motifs | cut -f 1-2,4 | sed 's/%//' | grep -v "Jaspar" | grep -v "withenhancer" > enhancer.others.motifs
grep -e "YY1" -e "CTCF(" -e "TATA-Box" -e "Sp1(" enhancer.motifs | cut -f 1-2,4 | sed 's/%//' | grep -v "Jaspar" | grep -v "withenhancer" > enhancer.others.1.motifs
Rscript /home1/04935/shaojf/myTools/BioinformaticsDaily/RVisualization/barplot4threecolsdata.R enhancer.ETS.motifs
Rscript /home1/04935/shaojf/myTools/BioinformaticsDaily/RVisualization/barplot4threecolsdata.R enhancer.others.motifs
Rscript /home1/04935/shaojf/myTools/BioinformaticsDaily/RVisualization/barplot4threecolsdata.R enhancer.others.1.motifs


cat <(awk -v OFS="\t" '{print $2,"E."$1,$4}' enhancer.motifs) <(awk -v OFS="\t" '{print $2,"P."$1,$4}' promoter.motifs) | grep -v "Motif Name" | sed 's/%//' | perl /home1/04935/shaojf/myTools/BioinformaticsDaily/textProcess/make_matrixwith3col_from_single_file.pl /dev/stdin > EP.motif.mat
head -1 EP.motif.mat > sim.EP.motif.mat
grep -e "ETS" -e "NRF" -e "bHLH" -e "Zf" -e "CTCF" -e "YY1" EP.motif.mat >> sim.EP.motif.mat
################################################################################################
# library(pheatmap)
# library(data.table)
# # args <- c("EP.motif.mat")
# args <- c("sim.EP.motif.mat")
# input <- fread(args[1], sep = "\t", header = T, na.strings = "/")
# datax <- data.matrix(input[,-1])
# rownames(datax) <- as.matrix(input[,1])[,1]
# colors <- colorRampPalette(c("white", "blue"))(100)

# nullcount <- apply(datax, 1, function(x){length(x[x==0])})
# pdf(file = paste(args[1], "pheatmap.pdf", sep = "."), width = 10, height = 16)
# myplot <- pheatmap(datax[nullcount == 0, ], scale = "none", 
# 	# show_rownames = T, show_colnames = T, color = colors, fontsize_row = 3,
# 	show_rownames = T, show_colnames = T, color = colors, fontsize_row = 5,
# 	cluster_cols = T, cluster_rows = T)
# dev.off()
# pdf(file = paste(args[1], "manually.pheatmap.pdf", sep = "."), width = 10, height = 16)
# myplot <- pheatmap(datax[nullcount == 0, ], scale = "row", 
# 	show_rownames = T, show_colnames = T, color = colors, fontsize_row = 5,
# 	cluster_cols = F, cluster_rows = T)
# dev.off()
# pdf(file = paste(args[1], "manually.pheatmap.1.pdf", sep = "."), width = 10, height = 6)
# myplot <- pheatmap(datax[nullcount == 0 & (datax[,13]-datax[,12]>10 | datax[,13]-datax[,11]>10), ], scale = "none", 
# 	show_rownames = T, show_colnames = T, color = colors, fontsize_row = 10,
# 	cluster_cols = F, cluster_rows = T)
# dev.off()
# pdf(file = paste(args[1], "manually.pheatmap.2.pdf", sep = "."), width = 4, height = 6)
# myplot <- pheatmap(datax[nullcount == 0 & (abs(datax[,13]-datax[,12])>5 | abs(datax[,13]-datax[,11])>5), 8:14], scale = "none", 
# 	show_rownames = T, show_colnames = T, color = colors, fontsize_row = 4,
# 	cluster_cols = T, cluster_rows = T)
# dev.off()
# pdf(file = paste(args[1], "scaled.pheatmap.pdf", sep = "."), width = 10, height = 16)
# myplot <- pheatmap(datax[nullcount == 0, ], scale = "row", 
# 	show_rownames = T, show_colnames = T, color = colors, fontsize_row = 3,
# 	cluster_cols = T, cluster_rows = T)
# dev.off()
################################################################################################






for gene in *HiChIP-Loop.bed
do
	class=`echo $gene | sed 's/.bed//'`
	for post in GM12878 K562 Naive Th17 Treg
	do
		perl $myperl $pre.$post.enhancers <(cut -f 4 $gene) 0 0 | cut -f 1,3 | grep -v "/" > $class.$post.enhancers
	done
done
for gene in *HiChIP-Loop.bed
do
	class=`echo $gene | sed 's/.bed//'`
	for post in GM12878 K562 Naive Th17 Treg
	do
		findMotifsGenome.pl <(sed 's/:/\t/;s/-/\t/' $class.$post.enhancers | awk -vOFS="\t" '{print "chr"$2,$3,$4,$1}') hg19 $class.$post.enhancers.homer.motifs -size given -mknown /home1/04935/shaojf/myTools/HOMER/data/knownTFs/vertebrates/known.motifs -p 60 1> $class.$post.enhancers.findMotifsGenome.txt 2>&1 &
		annotatePeaks.pl <(sed 's/:/\t/;s/-/\t/' $class.$post.enhancers | awk -vOFS="\t" '{print "chr"$2,$3,$4,$1}') hg19 -m /home1/04935/shaojf/myTools/HOMER/data/knownTFs/vertebrates/known.motifs -size given -cpu 60 1> $class.$post.enhancers.vert.motifs.txt 2> $class.$post.enhancers.vert.motifs.log &
	done
done
#####################################################################


#########
# dotplot for motif compare
head_line hkg.withoutHiChIP-Loop.promoter.vert.motifs.txt | awk '{print $2}' | tail -n +22 > known.vert.motif.list
echo "motifs qvalue PercentageWithMotif class" | tr " " "\t" > known.vert.motif.tsv
for f in `ls *HiChIP-Loop.homer.motifs/knownResults.txt`
do
	pre=`echo $f | sed 's?.homer.motifs/knownResults.txt??'`
	awk -F"\t" -v OFS="\t" -v var=$pre '$5<0.01{print $1,$5,$7,var}' $f | grep -wf known.vert.motif.list | sed 's/%//' >> known.vert.motif.tsv
done
awk '$4!~/appris/' known.vert.motif.tsv | sed 's?/Homer??' > selected.known.vert.motif.tsv

for f in `ls *HiChIP-Loop.homer.motifs/knownResults.txt | grep -v "appris"`
do
	awk -F"\t" '$5<0.01{print $1}' $f | grep -wf known.vert.motif.list | head -10 >> top.known.vert.motif.tmp
done

echo "motifs qvalue PercentageWithMotif class" | tr " " "\t" > top.known.vert.motif.tsv
for f in `ls *HiChIP-Loop.homer.motifs/knownResults.txt | grep -v "appris"`
do
	pre=`echo $f | sed 's?.homer.motifs/knownResults.txt??'`
	awk -F"\t" -v OFS="\t" -v var=$pre '{print $1,$5,$7,var}' $f | grep -wf <(sort top.known.vert.motif.tmp | uniq) | sed 's/%//;s?/Homer??' >> top.known.vert.motif.tsv
done
# head_line hkg.gm.liver.promoter.vert.motifs.mat | grep -wf <(cut -f1-2 -d"/" top.known.vert.motif.tmp) | cut -f1 -d" " | tr "\n" ","
cut -f 1,38,54,55,56,57,69,72,73,83,89,95,108,121,122,134,136,149,208,210,255,268,291,294,299,334 hkg.gm.liver.promoter.vert.motifs.mat > hkg.gm.liver.promoter.vert.motifs.top.mat
################################################################################################
# library(data.table)
# library(ggplot2)
# # args <- c("known.vert.motif.tsv")
# # args <- c("selected.known.vert.motif.tsv")
# args <- c("top.known.vert.motif.tsv")
# input <- fread(args[1], sep = "\t", header = T, na.strings = "/")
# # png(filename = paste(args[1], "enrichment.plot.png", sep = "."), width = 1000, height = 1600)
# # pdf(file = paste(args[1], "enrichment.plot.pdf", sep = "."), width = 10, height = 20)
# pdf(file = paste(args[1], "enrichment.plot.pdf", sep = "."), width = 10, height = 10)
# myplot <- ggplot(data = input, aes(x = class, y = motifs)) + 
# 		geom_point(aes(colour = qvalue, size = PercentageWithMotif)) +
# 		# scale_colour_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0.005, space = "Lab", na.value = "grey50", guide = "colourbar") + 
# 		scale_colour_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0.5, space = "Lab", na.value = "grey50", guide = "colourbar") + 
# 		labs(title = args[1], caption = date()) + 
# 		theme(axis.text.x = element_text(angle = 60, hjust = 1))
# print(myplot)
# dev.off()
# pdf(file = paste(args[1], "sim.enrichment.plot.pdf", sep = "."), width = 10, height = 15)
# myplot <- ggplot(data = input[class != "hkg.someHiChIP-Loop"], aes(x = class, y = motifs)) + 
# 		geom_point(aes(colour = qvalue, size = PercentageWithMotif)) +
# 		scale_colour_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0.005, space = "Lab", na.value = "grey50", guide = "colourbar") + 
# 		labs(title = args[1], caption = date()) + 
# 		theme(axis.text.x = element_text(angle = 60, hjust = 1))
# print(myplot)
# dev.off()
################################################################################################

cut -f 1-13,18-20,23,25-26 hkg.gm.liver.promoter.vert.motifs.top.mat > sim.hkg.gm.liver.promoter.vert.motifs.top.mat
################################################################################################
# library(pheatmap)
# library(data.table)
# # args <- c("hkg.gm.liver.promoter.vert.motifs.top.mat")
# args <- c("sim.hkg.gm.liver.promoter.vert.motifs.top.mat")
# input <- fread(args[1], sep = "\t", header = T, na.strings = "/")
# datax <- data.matrix(input[,-1])
# rownames(datax) <- as.matrix(input[,1])[,1]
# datax[datax > 5] <- 5

# annos <- data.frame(class = c(rep("hkg+",302), rep("hkg+-",2821), rep("hkg-",1058), 
# 	rep("gm+",117), rep("gm-", 96), rep("liver+",53), rep("liver-", 265)))
# rownames(annos) <- rownames(datax)
# colors <- colorRampPalette(c("white", "blue"))(10)

# pdf(file = paste(args[1], "pheatmap.pdf", sep = "."), width = 10, height = 8)
# myplot <- pheatmap(datax, scale = "none", annotation_row = annos,
# 	show_rownames = F, show_colnames = T, color = colors, fontsize_col = 8,
# 	cluster_cols = T, cluster_rows = F)
# dev.off()
# pdf(file = paste(args[1], "pheatmap.1.pdf", sep = "."), width = 10, height = 8)
# myplot <- pheatmap(datax, scale = "none", annotation_row = annos,
# 	show_rownames = F, show_colnames = T, color = colors, fontsize_col = 8,
# 	cluster_cols = T, cluster_rows = T)
# dev.off()
# datay <- datax[1:4181,]
# annos.sim <- data.frame(class = c(rep("hkg+",302), rep("hkg+-",2821), rep("hkg-",1058)))
# rownames(annos.sim) <- rownames(datay)
# pdf(file = paste(args[1], "pheatmap.2.pdf", sep = "."), width = 10, height = 8)
# myplot <- pheatmap(datay, scale = "none", annotation_row = annos.sim,
# 	show_rownames = F, show_colnames = T, color = colors, fontsize_col = 8,
# 	cluster_cols = T, cluster_rows = T)
# dev.off()
# datay <- datax[c(1:302,3124:4181),]
# annos.sim <- data.frame(class = c(rep("hkg+",302), rep("hkg-",1058)))
# rownames(annos.sim) <- rownames(datay)
# pdf(file = paste(args[1], "pheatmap.3.pdf", sep = "."), width = 10, height = 8)
# myplot <- pheatmap(datay, scale = "none", annotation_row = annos.sim,
# 	show_rownames = F, show_colnames = T, color = colors, fontsize_col = 8,
# 	cluster_cols = T, cluster_rows = T)
# dev.off()
# datay <- datax[c(1:302,3124:4181),c(9,11,15:18)]
# annos.sim <- data.frame(class = c(rep("hkg+",302), rep("hkg-",1058)))
# rownames(annos.sim) <- rownames(datay)
# pdf(file = paste(args[1], "pheatmap.4.pdf", sep = "."), width = 6, height = 8)
# myplot <- pheatmap(datay, scale = "none", annotation_row = annos.sim,
# 	show_rownames = F, show_colnames = T, color = colors, fontsize_col = 8,
# 	cluster_cols = T, cluster_rows = T)
# dev.off()
# pdf(file = paste(args[1], "pheatmap.5.pdf", sep = "."), width = 6, height = 8)
# myplot <- pheatmap(datay, scale = "none", annotation_row = annos.sim,
# 	show_rownames = F, show_colnames = T, color = colors, fontsize_col = 8,
# 	cluster_cols = T, cluster_rows = F)
# dev.off()
################################################################################################


# for gene in gm.with*HiChIP-Loop.bed
# do
# 	pre=`echo $gene | sed 's/.bed//'`
# 	findMotifsGenome.pl <(awk -vOFS="\t" '{a=$2-1000;b=$2+1000;if($6=="-"){a=$3-1000;b=$3+1000}if(a<0){a=0;}print $1,a,b,$4,"1000",$6}' $gene) hg19 $pre.homer.motifs -size given -mknown /home1/04935/shaojf/myTools/HOMER/data/knownTFs/all.motifs -p 30 > $pre.findMotifsGenome.txt &
# done
for f in *HiChIP-Loop.promoter.vert.motifs.txt
do
	pre=`echo $f | sed 's/.txt//'`
	perl ~/myTools/UTHealthMut2Loop/relatedScripts/make.matrix.from.homer.motif.file.pl <(cut -f 1,22- $f | sed '1s?/Homer Distance From Peak(sequence,strand,conservation)??g') > $pre.mat
done
f=hkg.someHiChIP-Loop.promoter.vert.motifs.txt
pre=`echo $f | sed 's/.txt//'`
perl ~/myTools/UTHealthMut2Loop/relatedScripts/make.matrix.from.homer.motif.file.pl <(cut -f 1,22- $f | sed '1s?/Homer Distance From Peak(sequence,strand,conservation)??g') > $pre.mat

cat hkg.withHiChIP-Loop.promoter.vert.motifs.mat <(tail -n +2 hkg.someHiChIP-Loop.promoter.vert.motifs.mat) <(tail -n +2 hkg.withoutHiChIP-Loop.promoter.vert.motifs.mat) <(tail -n +2 gm.withHiChIP-Loop.promoter.vert.motifs.mat) <(tail -n +2 gm.withoutHiChIP-Loop.promoter.vert.motifs.mat) <(tail -n +2 liver.withHiChIP-Loop.promoter.vert.motifs.mat) <(tail -n +2 liver.withoutHiChIP-Loop.promoter.vert.motifs.mat) > hkg.gm.liver.promoter.vert.motifs.mat
################################################################################################
# library(pheatmap)
# library(data.table)
# args <- c("hkg.gm.liver.promoter.vert.motifs.mat")
# input <- fread(args[1], sep = "\t", header = T, na.strings = "/")
# datax <- data.matrix(input[,-1])
# rownames(datax) <- as.matrix(input[,1])[,1]
# datax[datax > 10] <- 10

# annos <- data.frame(class = c(rep("hkg+",302), rep("hkg+-",2821), rep("hkg-",1058), 
# 	rep("gm+",117), rep("gm-", 96), rep("liver+",53), rep("liver-", 265)))
# rownames(annos) <- rownames(datax)
# colors <- colorRampPalette(c("white", "blue"))(10)

# # png(filename = paste(pre, "pheatmap.png", sep = "."), width = 1000, height = 800)
# pdf(file = paste(args[1], "pheatmap.pdf", sep = "."), width = 20, height = 12)
# myplot <- pheatmap(datax, scale = "none", annotation_row = annos,
# 	show_rownames = F, show_colnames = T, color = colors, fontsize_col = 2,
# 	cluster_cols = T, cluster_rows = F)
# dev.off()
# # colsums <- apply(datax[1:1360,], 2, sum)
# nullcount <- apply(datax, 2, function(x){length(x[x==0])})
# pdf(file = paste(args[1], "selected.pheatmap.pdf", sep = "."), width = 20, height = 12)
# myplot <- pheatmap(datax[, nullcount < 3000], scale = "none", annotation_row = annos,
# 	show_rownames = F, show_colnames = T, color = colors, fontsize_col = 6,
# 	cluster_cols = T, cluster_rows = F)
# dev.off()
# pdf(file = paste(args[1], "selected.pheatmap.1.pdf", sep = "."), width = 20, height = 12)
# myplot <- pheatmap(datax[, nullcount < 3000], scale = "none", annotation_row = annos,
# 	show_rownames = F, show_colnames = T, color = colors, fontsize_col = 6,
# 	cluster_cols = T, cluster_rows = T)
# dev.off()
# annos.sim <- data.frame(class = c(rep("hkg+",302), rep("hkg-",1058)))
# rownames(annos.sim) <- rownames(datax)[c(1:302,3124:4181)]
# pdf(file = paste(args[1], "onlyhkg.selected.pheatmap.pdf", sep = "."), width = 20, height = 12)
# myplot <- pheatmap(datax[c(1:302,3124:4181), nullcount < 3000], scale = "none", annotation_row = annos.sim,
# 	show_rownames = F, show_colnames = T, color = colors, fontsize_col = 6,
# 	cluster_cols = T, cluster_rows = T)
# dev.off()
# datay <- datax[c(1:302,3124:4181), ]
# nullcount.y <- apply(datay, 2, function(x){length(x[x==0])})
# pdf(file = paste(args[1], "onlyhkg.selected.pheatmap.1.pdf", sep = "."), width = 20, height = 12)
# myplot <- pheatmap(datay[, nullcount.y < 1000], scale = "none", annotation_row = annos.sim,
# 	show_rownames = F, show_colnames = T, color = colors, fontsize_col = 6,
# 	cluster_cols = T, cluster_rows = F)
# dev.off()
################################################################################################
# perl ~/myTools/UTHealthMut2Loop/relatedScripts/motif.dis.from.homer.pl <(perl $myperl gencode.wg.promoter.motifs.sim hkg.tsg.srtbyPCA.class 0 0 | awk '$2~/hkg/ || NR==1' | cut -f 3- | sed 's/),/);/g') > hkg.promotermotif.dis.txt
#######################


############
# YY1
# cd /home1/04935/shaojf/stampede2/housekeeping_genes/both.pc.and.nc.genes/YY1_ENCODE
# rclone sync /home1/04935/shaojf/stampede2/housekeeping_genes/both.pc.and.nc.genes/YY1_ENCODE/ mygoogle:hkg_tsg/both.pc.and.nc.genes/YY1_ENCODE
gene=hkg.tsg.srtbyPCA.transcript.bed
echo "gene cell accession counts" | tr " " "\t" > hkg.tsg.srtbyPCA.transcript.YY1.stats
for f in hg19.*.bed.gz
do
	cell=`echo $f | cut -f2 -d"."`
	acc=`echo $f | cut -f4 -d"."`
	bedtools intersect -wo -a <(cut -f 3-  $gene | awk -vOFS="\t" '{a=$2-1000;b=$2+1000;if($6=="-"){a=$3-1000;b=$3+1000}if(a<0){a=0;}print $1,a,b,$4,"1000",$6}') -b <(gunzip -c $f | cut -f 1-3) > hkg.tsg.srtbyPCA.transcript.$cell.$acc.overlap.txt
	cut -f 4 hkg.tsg.srtbyPCA.transcript.$cell.$acc.overlap.txt | sort | uniq -c | awk -vOFS="\t" -v c=$cell -v a=$acc '{print $2,c,a,$1}' >> hkg.tsg.srtbyPCA.transcript.YY1.stats
done
perl make.peaks.matrix.pl <(awk -vOFS="\t" '{print $1,$2"."$3,$4}' hkg.tsg.srtbyPCA.transcript.YY1.stats | tail -n +2) > hkg.tsg.srtbyPCA.transcript.YY1.mat
head -1 hkg.tsg.srtbyPCA.transcript.YY1.mat | awk '{print $0"\tYY1.motif\tETS(ETS)/Promoter"}' > hkg.gm.liver.YY1.mat
perl $myperl hkg.tsg.srtbyPCA.transcript.YY1.mat <(cut -f1 hkg.gm.liver.promoter.vert.motifs.mat | tail -n +2) 0 0 | cut -f 1,3- | sed 's?/?0?g' > a
paste a <(cut -f 334 hkg.gm.liver.promoter.vert.motifs.mat | tail -n +2) <(cut -f 83 hkg.gm.liver.promoter.vert.motifs.mat | tail -n +2) >> hkg.gm.liver.YY1.mat
################################################################################################
# library(pheatmap)
# library(data.table)
# args <- c("hkg.gm.liver.YY1.mat")
# input <- fread(args[1], sep = "\t", header = T, na.strings = "/")
# datax <- data.matrix(input[,-1])
# rownames(datax) <- as.matrix(input[,1])[,1]
# datax[datax > 10] <- 10

# annos <- data.frame(class = c(rep("hkg+",302), rep("hkg+-",2821), rep("hkg-",1058), 
# 	rep("gm+",117), rep("gm-", 96), rep("liver+",53), rep("liver-", 265)))
# rownames(annos) <- rownames(datax)
# colors <- colorRampPalette(c("white", "blue"))(10)

# pdf(file = paste(args[1], "pheatmap.pdf", sep = "."), width = 10, height = 8)
# myplot <- pheatmap(datax, scale = "none", annotation_row = annos,
# 	show_rownames = F, show_colnames = T, color = colors, fontsize_col = 6,
# 	cluster_cols = T, cluster_rows = F)
# dev.off()
# nullcount <- apply(datax, 2, function(x){length(x[x==0])})
# pdf(file = paste(args[1], "selected.pheatmap.pdf", sep = "."), width = 10, height = 8)
# myplot <- pheatmap(datax[, nullcount < 4000], scale = "none", annotation_row = annos,
# 	show_rownames = F, show_colnames = T, color = colors, fontsize_col = 6,
# 	cluster_cols = T, cluster_rows = F)
# dev.off()
################################################################################################
head -1 hkg.tsg.srtbyPCA.transcript.YY1.mat > hkg.tsg.srtbyPCA.transcript.YY1.peaks.mat
perl $myperl hkg.tsg.srtbyPCA.transcript.YY1.mat <(cut -f6 hkg.tsg.srtbyPCA.transcript.bed) 0 0 | cut -f 1,3- | sed 's?/?0?g' >> hkg.tsg.srtbyPCA.transcript.YY1.peaks.mat
################################################################################################
# library(pheatmap)
# library(data.table)
# args <- c("hkg.tsg.srtbyPCA.transcript.YY1.peaks.mat", "hkg.tsg.srtbyPCA.transcript.bed")
# input <- fread(args[1], sep = "\t", header = T, na.strings = "/")
# datax <- data.matrix(input[,-1])
# rownames(datax) <- as.matrix(input[,1])[,1]
# datax[datax > 10] <- 10

# annos <- fread(args[2], sep = "\t", header = F, na.strings = "/")[,2]
# rownames(annos) <- rownames(datax)
# colors <- colorRampPalette(c("white", "blue"))(10)

# pdf(file = paste(args[1], "pheatmap.pdf", sep = "."), width = 10, height = 12)
# myplot <- pheatmap(datax, scale = "none", annotation_row = annos,
# 	show_rownames = F, show_colnames = T, color = colors, fontsize_col = 6,
# 	cluster_cols = T, cluster_rows = F)
# dev.off()
################################################################################################

############
# NRF1
# cd /home1/04935/shaojf/stampede2/housekeeping_genes/both.pc.and.nc.genes/NRF1_ENCODE
# rclone sync /home1/04935/shaojf/stampede2/housekeeping_genes/both.pc.and.nc.genes/NRF1_ENCODE/ mygoogle:hkg_tsg/both.pc.and.nc.genes/NRF1_ENCODE
gene=hkg.tsg.srtbyPCA.transcript.bed
echo "gene cell accession counts" | tr " " "\t" > hkg.tsg.srtbyPCA.transcript.NRF1.stats
for f in hg19.*.bed.gz
do
	cell=`echo $f | cut -f2 -d"."`
	acc=`echo $f | cut -f4 -d"."`
	bedtools intersect -wo -a <(cut -f 3-  $gene | awk -vOFS="\t" '{a=$2-1000;b=$2+1000;if($6=="-"){a=$3-1000;b=$3+1000}if(a<0){a=0;}print $1,a,b,$4,"1000",$6}') -b <(gunzip -c $f | cut -f 1-3) > hkg.tsg.srtbyPCA.transcript.$cell.$acc.overlap.txt
	cut -f 4 hkg.tsg.srtbyPCA.transcript.$cell.$acc.overlap.txt | sort | uniq -c | awk -vOFS="\t" -v c=$cell -v a=$acc '{print $2,c,a,$1}' >> hkg.tsg.srtbyPCA.transcript.NRF1.stats
done
perl make.peaks.matrix.pl <(awk -vOFS="\t" '{print $1,$2"."$3,$4}' hkg.tsg.srtbyPCA.transcript.NRF1.stats | tail -n +2) > hkg.tsg.srtbyPCA.transcript.NRF1.mat
head -1 hkg.tsg.srtbyPCA.transcript.NRF1.mat | awk '{print $0"\tNRF1(NRF)\tNRF(NRF)/Promoter\tYY1.motif\tETS(ETS)/Promoter"}' > hkg.gm.liver.NRF1.mat
perl $myperl hkg.tsg.srtbyPCA.transcript.NRF1.mat <(cut -f1 hkg.gm.liver.promoter.vert.motifs.mat | tail -n +2) 0 0 | cut -f 1,3- | sed 's?/?0?g' > a
paste a <(cut -f 208 hkg.gm.liver.promoter.vert.motifs.mat | tail -n +2) <(cut -f 210 hkg.gm.liver.promoter.vert.motifs.mat | tail -n +2) <(cut -f 334 hkg.gm.liver.promoter.vert.motifs.mat | tail -n +2) <(cut -f 83 hkg.gm.liver.promoter.vert.motifs.mat | tail -n +2) >> hkg.gm.liver.NRF1.mat
################################################################################################
library(pheatmap)
library(data.table)
args <- c("hkg.gm.liver.NRF1.mat")
input <- fread(args[1], sep = "\t", header = T, na.strings = "/")
datax <- data.matrix(input[,-1])
rownames(datax) <- as.matrix(input[,1])[,1]
datax[datax > 1] <- 1

annos <- data.frame(class = c(rep("hkg+",302), rep("hkg+-",2821), rep("hkg-",1058), 
	rep("gm+",117), rep("gm-", 96), rep("liver+",53), rep("liver-", 265)))
rownames(annos) <- rownames(datax)
colors <- colorRampPalette(c("white", "blue"))(10)

pdf(file = paste(args[1], "pheatmap.pdf", sep = "."), width = 10, height = 8)
myplot <- pheatmap(datax, scale = "none", annotation_row = annos,
	show_rownames = F, show_colnames = T, color = colors, fontsize_col = 6,
	cluster_cols = T, cluster_rows = F)
dev.off()
nullcount <- apply(datax, 2, function(x){length(x[x==0])})
pdf(file = paste(args[1], "selected.pheatmap.pdf", sep = "."), width = 10, height = 8)
myplot <- pheatmap(datax[, nullcount < 4000], scale = "none", annotation_row = annos,
	show_rownames = F, show_colnames = T, color = colors, fontsize_col = 6,
	cluster_cols = F, cluster_rows = F)
dev.off()
datay <- datax[1:4181,16:19]
annos.sim <- data.frame(class = c(rep("hkg+",302), rep("hkg+-",2821), rep("hkg-",1058)))
rownames(annos.sim) <- rownames(datay)
pdf(file = paste(args[1], "selected.motifs.pheatmap.pdf", sep = "."), width = 4, height = 6)
myplot <- pheatmap(datay, scale = "none", annotation_row = annos.sim,
	show_rownames = F, show_colnames = T, color = colors, fontsize_col = 6,
	cluster_cols = T, cluster_rows = T)
dev.off()
################################################################################################
head -1 hkg.tsg.srtbyPCA.transcript.NRF1.mat > hkg.tsg.srtbyPCA.transcript.NRF1.peaks.mat
perl $myperl hkg.tsg.srtbyPCA.transcript.NRF1.mat <(cut -f6 hkg.tsg.srtbyPCA.transcript.bed) 0 0 | cut -f 1,3- | sed 's?/?0?g' >> hkg.tsg.srtbyPCA.transcript.NRF1.peaks.mat
################################################################################################
# library(pheatmap)
# library(data.table)
# args <- c("hkg.tsg.srtbyPCA.transcript.NRF1.peaks.mat", "hkg.tsg.srtbyPCA.transcript.bed")
# input <- fread(args[1], sep = "\t", header = T, na.strings = "/")
# datax <- data.matrix(input[,-1])
# rownames(datax) <- as.matrix(input[,1])[,1]
# datax[datax > 5] <- 5

# annos <- fread(args[2], sep = "\t", header = F, na.strings = "/")[,2]
# rownames(annos) <- rownames(datax)
# colors <- colorRampPalette(c("white", "blue"))(10)

# pdf(file = paste(args[1], "pheatmap.pdf", sep = "."), width = 10, height = 12)
# myplot <- pheatmap(datax, scale = "none", annotation_row = annos,
# 	show_rownames = F, show_colnames = T, color = colors, fontsize_col = 6,
# 	cluster_cols = T, cluster_rows = F)
# dev.off()
################################################################################################










### paste tpm value into one matrix 
cut -f 1-2 A172__transcript.quantifications.ENCFF185SEV.tsv | awk '{print $2"|"$1}' > encode.transcript.tsv
ls *transcript*.tsv | grep -v "dendritic.cell" | while read files
do
	echo $files > a
	cut -f 6 $files | tail -n +2 >> a
	paste encode.transcript.tsv a > b
	mv b encode.transcript.tsv
	rm a
done

### prepare the tss for each transcript
# gunzip -c gencode.v19.annotation.gtf.gz | awk '$3=="transcript"' | awk -vOFS="\t" '{print $1,$4,$5,$10"|"$12,".",$7}' | sed 's/"//g;s/;//g' > gencode.v19.transcript.bed
# login2.stampede2(1052)$ gunzip -c gencode.v19.annotation.gtf.gz | awk '$3=="transcript"' | grep "appris_principal" | grep "level 1" | wc -l
# 464
# login2.stampede2(1053)$ gunzip -c gencode.v19.annotation.gtf.gz | awk '$3=="transcript"' | grep "appris_principal" | grep -e "level 1" -e "level 2"| wc -l
# 21269
# login2.stampede2(1021)$ gunzip -c gencode.v19.annotation.gtf.gz | awk '$3=="transcript"' | grep -w -e "appris_principal" -e "appris_candidate" | wc -l
# 28036
# gunzip -c gencode.v19.annotation.gtf.gz | awk '$3=="transcript"' | grep "appris_principal" | grep -e "level 1" -e "level 2" | awk -vOFS="\t" '{print $1,$4,$5,$10"|"$12,".",$7}' | sed 's/"//g;s/;//g' > gencode.v19.transcript.appris_principal.level12.bed
# ln -s gencode.v19.transcript.appris_principal.level12.bed gencode.v19.transcript.bed
gunzip -c gencode.v19.annotation.gtf.gz | awk '$3=="transcript"' | grep -w -e "appris_principal" -e "appris_candidate" | awk -vOFS="\t" '{print $1,$4,$5,$10"|"$12,".",$7}' | sed 's/"//g;s/;//g' > gencode.v19.transcript.appris.bed
ln -s gencode.v19.transcript.appris.bed gencode.v19.transcript.bed
perl $myperl <(sed 's/|/\t/' gencode.v19.transcript.bed) <(sed 's/|/\t/' hkg.tsg.srtbyPCA.class) 3 0 | awk -F"\t" -vOFS="\t" '{print $1"|"$2,$3,$4,$5,$6,$7"|"$8,$9,$10}' | tail -n +2 | grep "/" > hkg.tsg.srtbyPCA.no.transcript
perl $myperl <(sed 's/|/\t/' gencode.v19.transcript.bed) <(sed 's/|/\t/' hkg.tsg.srtbyPCA.class) 3 0 | awk -F"\t" -vOFS="\t" '{print $1"|"$2,$3,$4,$5,$6,$7"|"$8,$9,$10}' | tail -n +2 | grep -v "/" > hkg.tsg.srtbyPCA.transcript.bed
cut -f 1-2,6 hkg.tsg.srtbyPCA.transcript.bed > hkg.tsg.srtbyPCA.transcript.class

### find the overlapped features for transcript
bedtools intersect -wo -a hg19.cpgIslandExt.withCount.bed -b <(awk -vOFS="\t" '{a=$2-300;b=$2+100;if($6=="-"){a=$3-100;b=$3+300}if(a<0){a=0;}print $1,a,b,$4,$5,$6}' gencode.v19.transcript.bed) > gencode.cpgIslandExt.txt

### for ENCODE ChIP-seq
ls ../hg19.released.files/*peaks.* | while read files
do
	exp=`echo $files | awk -F"/" '{print $3}' | sed 's/.bed.gz//'`
	bedtools intersect -wo -a <(gunzip -c $files | awk -vOFS="\t" -v name=$exp '{print $1,$2,$3,name}') -b <(awk -vOFS="\t" '{a=$2-300;b=$2+100;if($6=="-"){a=$3-100;b=$3+300}if(a<0){a=0;}print $1,a,b,$4,$5,$6}' gencode.v19.transcript.bed) > gencode.peaks.$exp &
	sleep 0.1s
done

# cat gencode.peaks.* | cut -f 4,8 | sort | uniq -c | awk -vOFS="\t" '{print $2,$3,$1}' | perl /home1/04935/shaojf/stampede2/myTools/UTHealthMut2Loop/relatedScripts/make.occurence_count.matrix.from.GTRD.pl /dev/stdin > ../formated.files/encode.peaks.tsv
# perl ../make.peaks.matrix.pl <(cat gencode.peaks.*) > ../formated.files/encode.peaks.tsv 

for f in `ls gencode.peaks.*`
do
	# cut -f 4,8 $f | sort | uniq -c | awk -vOFS="\t" '{print $2,$3,$1}' >> encode.peaks.txt
	cut -f 4,8 $f | sort | uniq -c | awk -vOFS="\t" '{print $2,$3,$1}' > tmp.$f &
	sleep 10s
done
cat tmp.gencode.peaks.* | perl /home1/04935/shaojf/stampede2/myTools/UTHealthMut2Loop/relatedScripts/make.occurence_count.matrix.from.GTRD.pl /dev/stdin > ../formated.files/encode.peaks.tsv


### ChIA-PET
# awk -F"\t" -v OFS="\t" '$2!="bam" && $2!~"bigBed" && $2!="fastq" && $2!="bigWig" && $39=="hg19" && $42=="released"{print $38,$7"_"$13"_"$3}' metadata.tsv | awk -F"/" '{print $7}' | sed 's/ /./g' | awk '{print $0"."$1}' |  while read line
# do
# 	old=`echo $line | awk '{print $1}'`
# 	new=`echo $line | awk '{print $2}'`
# 	if [ -e $old ]
# 	then
# 		mv $old $new
# 	fi
# done
ls *.bed.gz | while read files
do
	exp=`echo $files | sed 's/.bed.gz//'`
	bedtools intersect -wo -a <(gunzip -c $files | awk -vOFS="\t" -v name=$exp '{print $1,$2,$3,name}') -b <(awk -vOFS="\t" '{a=$2-300;b=$2+100;if($6=="-"){a=$3-100;b=$3+300}if(a<0){a=0;}print $1,a,b,$4,$5,$6}' gencode.v19.transcript.bed) > gencode.loop.$exp &
	sleep 0.1s
done
perl ../make.peaks.matrix.pl <(cat gencode.loop.*) > ../formated.files/encode.loop.tsv 



##### finalize the matrix
perl $myperl encode.transcript.tsv hkg.tsg.srtbyPCA.transcript.class 0 2 | cut -f 4- > hkg.tsg.srtbyPCA.transcript.tpm
perl $myperl encode.peaks.tsv hkg.tsg.srtbyPCA.transcript.class 0 2 | cut -f 4- > hkg.tsg.srtbyPCA.transcript.peaks
perl $myperl encode.loop.tsv hkg.tsg.srtbyPCA.transcript.class 0 2 | cut -f 4- > hkg.tsg.srtbyPCA.transcript.loop

