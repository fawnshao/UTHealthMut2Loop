#!/bin/sh
# rclone sync /home1/04935/shaojf/stampede2/housekeeping_genes/v2/ mygoogle:hkg_tsg/v2/
# cut -f 1,100 v2.TSG.tsv | awk '$2=="Liver"{print $1}' | cut -f 2 -d "|" > v2.TSG.Liver.singleTSG.gene.txt
# cut -f 1,100 v2.TSG.tsv | awk -F"\t" '$2=="Cells - EBV-transformed lymphocytes"{print $1}' | cut -f 2 -d "|" > v2.TSG.GM.singleTSG.gene.txt
# cut -f 1,100 v2.TSG.tsv | awk -F"\t" '$2=="Testis"{print $1}' | cut -f 2 -d "|" > v2.TSG.Testis.singleTSG.gene.txt
# cut -f 1,100 v2.TSG.tsv | grep "Brain" | awk -F"\t" '{print $1}' | cut -f 2 -d "|" > v2.TSG.Brain.singleTSG.gene.txt

# cat <(cut -f 1 v2.HKG.tsv | tail -n +2) <(cut -f 1 v2.HKG.2.tsv | tail -n +2) | perl ~/myTools/BioinformaticsDaily/textProcess/add_any_2files_together.pl hkg.tsg.srtbyPCA.class /dev/stdin 0 0 > v2.HKG.cmp
# cat <(cut -f 1 v2.HKG.tsv | tail -n +2) <(cut -f 1 v2.HKG.2.tsv | tail -n +2) | perl ~/myTools/BioinformaticsDaily/textProcess/add_any_2files_together.pl /dev/stdin hkg.tsg.srtbyPCA.class 0 0 > v2.HKG.cmp.2

# more v2.HKG.cmp.2 | grep hkg | grep / | cut -f 1 | grep -wf /dev/stdin v2.allstats.tsv > lost.in.v2.HKG
# more v2.HKG.cmp | grep / | cut -f 1 | grep -wf /dev/stdin v2.allstats.tsv > gained.in.v2.HKG 

myperl=/home1/04935/shaojf/myTools/BioinformaticsDaily/textProcess/add_any_2files_together.pl

# cd /home1/04935/shaojf/stampede2/housekeeping_genes/v2/test.top.hkg.tsg
echo "Gene Type" | sed 's/ /\t/' > v2.top.hkg.tsg.txt
cut -f 1 v2.HKG.tsv | tail -n +2 | awk '{print $0"\tHKG"}' >> v2.top.hkg.tsg.txt
cut -f 1,100 v2.TSG.tsv | tail -n +2 >> v2.top.hkg.tsg.txt

perl $myperl <(sed 's/\t/|/' gencode.v19.gene.anntotation.txt) v2.top.hkg.tsg.txt 0 0 | cut -f 1-2,4 > v2.top.hkg.tsg.anntotation.txt
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