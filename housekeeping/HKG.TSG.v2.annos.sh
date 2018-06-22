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

for tf in EHF ERG ETV1 ETV6 GABP PU.1 Spi-B
do
	findMotifsGenome.pl <(awk -F"\t" -vOFS="\t" '{a=$2-30;b=$3+30;if(a<0){a=0;}print $1,a,b,$4}' $tf.meta_clusters.bed) hg19 $tf.homer.motifs -size given -mknown /home1/04935/shaojf/myTools/HOMER/data/knownTFs/vertebrates/known.motifs 1> $tf.findMotifsGenome.txt 2>&1 &
	annotatePeaks.pl <(awk -F"\t" -vOFS="\t" '{a=$2-30;b=$3+30;if(a<0){a=0;}print $1,a,b,$4}' $tf.meta_clusters.bed) hg19 -m /home1/04935/shaojf/myTools/HOMER/data/knownTFs/vertebrates/known.motifs -size given -cpu 68 1> $tf.promoter.vert.motifs.txt 2> $tf.promoter.vert.motifs.log &
done

# cut -f 1 ../v2.HKG.2.tsv | tail -n +2 | awk '{print $0"\tHKG.2"}' > v2.hkg.2.txt
# perl $myperl gencode.v19.gene.bed v2.hkg.2.txt 3 0 | cut -f 3- > v2.hkg.2.gene.bed
for pre in HKG Testis Spleen Liver Brain EBV protein_coding gencode.v19 v2.hkg.2
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
