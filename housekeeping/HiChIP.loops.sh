#!/bin/sh
# rclone sync /home1/04935/shaojf/stampede2/housekeeping_genes/both.pc.and.nc.genes/HiChIP.loops mygoogle:hkg_tsg/both.pc.and.nc.genes/HiChIP.loops
# cd /home1/04935/shaojf/stampede2/housekeeping_genes/both.pc.and.nc.genes/HiChIP.loops
myperl=/home1/04935/shaojf/myTools/BioinformaticsDaily/textProcess/add_any_2files_together.pl
matperl=/home1/04935/shaojf/myTools/BioinformaticsDaily/textProcess/make_matrixwith3col_from_single_file.pl 
bedperl=/home1/04935/shaojf/myTools/UTHealthMut2Loop/housekeeping/get.HiChIP.loops.annotation.from.bedtools.intersect.pl
annoperl=/home1/04935/shaojf/myTools/UTHealthMut2Loop/housekeeping/add.annotation.to.a.combination.pl
for f in *_H3K27ac_Loops.txt
do
	pre=`echo $f | cut -f 1 -d"_"`
	tail -n +2 $f | awk -v var=$pre '{print $1":"$2"-"$3"%"$4":"$5"-"$6"\t"var"\t1"}' >> tmp.txt
done
perl $matperl tmp.txt > HiChIP.loops.mat
rm tmp.txt
# tail -n +2 HiChIP.loops.mat | awk -vOFS="\t" '{print $1,$1}' | sed 's/%/\t/' | awk '{print $1"\t"$3"\n"$2"\t"$3}' | uniq | sed 's/:/\t/;s/-/\t/' | awk '{print "chr"$0}' | bedtools sort -i - > HiChIP.loops.bed
tail -n +2 HiChIP.loops.mat | awk '{print $1}' | sed 's/%/\t/' | awk '{print $1"\t"$1"\n"$2"\t"$2}' | uniq | sed 's/:/\t/;s/-/\t/' | awk '{print "chr"$0}' | bedtools sort -i - > HiChIP.loops.bed
# bedtools intersect -wao -a <(awk -vOFS="\t" '{print $1,$2-5000,$3+5000,$4}' HiChIP.loops.bed) -b <(awk -F"\t" -vOFS="\t" '{a=$2-5000;b=$2+5000;if($6=="-"){a=$3-5000;b=$3+5000}if(a<0){a=0;}print $1,a,b,$4,$5,$6}' gencode.v19.transcript.bed | bedtools sort -i -) > HiChIP.loops.transcript
# bedtools intersect -wao -a HiChIP.loops.bed -b <(awk -F"\t" -vOFS="\t" '{a=$2-5000;b=$2+5000;if($6=="-"){a=$3-5000;b=$3+5000}if(a<0){a=0;}print $1,a,b,$4,$5,$6}' gencode.v19.transcript.bed | bedtools sort -i -) > HiChIP.loops.transcript
bedtools intersect -wao -a HiChIP.loops.bed -b <(awk -F"\t" -vOFS="\t" '{a=$2-1000;b=$2+1000;if($6=="-"){a=$3-1000;b=$3+1000}if(a<0){a=0;}print $1,a,b,$4,$5,$6}' gencode.v19.transcript.bed | bedtools sort -i -) > HiChIP.loops.transcript
perl $bedperl HiChIP.loops.transcript > HiChIP.loops.transcript.sim

head -1 HiChIP.loops.mat | awk '{print $0"\tGeneAnnotation"}' > HiChIP.loops.mat.annotation
perl $myperl HiChIP.loops.transcript.sim <(sed 's/%/\t/' HiChIP.loops.mat | tail -n +2) 0 0 | cut -f 1-7,9 | perl $myperl HiChIP.loops.transcript.sim /dev/stdin 0 1 | cut -f 1-8,10 | sed 's/\t/%/' >> HiChIP.loops.mat.annotation
# head -1 HiChIP.loops.mat.annotation > HiChIP.loops.mat.annotation.diffloop
# sed 's/%/\t/' HiChIP.loops.mat.annotation | tail -n +2 | awk -F"\t" '$1!=$2' >> HiChIP.loops.mat.annotation.diffloop

echo "ID cell.count GeneAnnotation1 GeneAnnotation2" | tr " " "\t" > HiChIP.loops.mat.annotation.sim
tail -n +2 HiChIP.loops.mat.annotation | awk -vOFS="\t" '{sum=0;for(i=2;i<=6;i++){sum+=$i}print $1,sum,$7,$8}' >> HiChIP.loops.mat.annotation.sim

head -1 HiChIP.loops.mat.annotation.sim > HiChIP.loops.mat.annotation.sim.hkg.tsg
perl $annoperl <(tail -n +2 hkg.tsg.srtbyPCA.class | sed 's/|/\t/' | awk -F"\t" -vOFS="\t" '{print $1,$2"|"$3}') <(tail -n +2 HiChIP.loops.mat.annotation.sim) >> HiChIP.loops.mat.annotation.sim.hkg.tsg

awk '$2>2' HiChIP.loops.mat.annotation.sim.hkg.tsg | grep hkg