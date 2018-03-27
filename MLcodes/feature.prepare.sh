#!/bin/sh
myperl=/home1/04935/shaojf/myScripts/add_any_2files_together.pl
pathology=/work/04935/shaojf/stampede2/refs/ProteinAtlas/pathology.tsv
myR=/home1/04935/shaojf/myTools/UTHealthMut2Loop/MLcodes/feature_generation_learning.R
myboxplot=/home1/04935/shaojf/myTools/UTHealthMut2Loop/MLcodes/features.boxplot.R
mybarplot=/home1/04935/shaojf/myTools/UTHealthMut2Loop/MLcodes/features.percentage.R
mystatR=/home1/04935/shaojf/myTools/UTHealthMut2Loop/MLcodes/output.percentage.R 

list=hkg.tsg.srtbyPCA.class
# rclone sync /home1/04935/shaojf/stampede2/housekeeping_genes/both.pc.and.nc.genes/ mygoogle:hkg_tsg/both.pc.and.nc.genes/
###### sequence features
ln -s /work/04935/shaojf/stampede2/refs/UCSC/hg19/hg19.cpgIslandExt.withCount.bed 
ln -s /work/04935/shaojf/stampede2/refs/UCSC/hg19/hg19.simpleRepeat.withFRepeatLengthandTimes.bed
ln -s /work/04935/shaojf/stampede2/refs/GENCODE/gencode.v19.gene.exoncount.txt 
ln -s /work/04935/shaojf/stampede2/refs/GENCODE/gencode.v19.gene.length.txt
bedtools intersect -wo -a hg19.cpgIslandExt.withCount.bed -b <(awk -vOFS="\t" '{a=$2-1000;b=$2+1000;if($6=="-"){a=$3-1000;b=$3+1000}if(a<0){a=0;}print $1,a,b,$4,$5,$6}' gencode.v19.gene.bed) > gencode.cpgIslandExt.txt
bedtools intersect -wo -a hg19.simpleRepeat.withFRepeatLengthandTimes.bed -b <(awk -vOFS="\t" '{a=$2-1000;b=$2+1000;if($6=="-"){a=$3-1000;b=$3+1000}if(a<0){a=0;}print $1,a,b,$4,$5,$6}' gencode.v19.gene.bed) > gencode.simpleRepeat.txt
perl ~/myTools/UTHealthMut2Loop/relatedScripts/sum_cpg_repeat_length.pl gencode.cpgIslandExt.txt > gencode.cpgIslandExt.length.txt
perl ~/myTools/UTHealthMut2Loop/relatedScripts/sum_cpg_repeat_length.pl gencode.simpleRepeat.txt > gencode.simpleRepeat.length.txt

awk '{print $4"\t"$3-$2+1}' gencode.v19.firstI.bed > gencode.v19.firstI.length
awk '{print $4"\t"$3-$2+1}' gencode.v19.firstE.bed | sed 's/:/\t/' | cut -f 1,3 > gencode.v19.firstE.length

### 5th: mean0 - average over bases with non-covered bases counting as zeroes
### 6th: mean - average over just covered bases
echo "Gene Length ExonCount ExonLength IntronLenth FirstExonLength FirstIntronLength FirstExonCons FirstIntronCons PromoterCons cpgIslandLength simpleRepeatLength CAGE.count" | tr " " "\t" > gencode.sequenceFeatures.txt
perl $myperl gencode.v19.gene.length.txt <(cut -f 4 gencode.v19.gene.bed) 0 0 | cut -f 1,3 | perl $myperl gencode.v19.gene.exoncount.txt /dev/stdin 0 0 | cut -f 1-2,4 | perl $myperl gencode.v19.exon.length.txt /dev/stdin 0 0 | cut -f 1-3,5 | perl $myperl gencode.v19.intron.length.txt /dev/stdin 0 0 | cut -f 1-4,6 | perl $myperl gencode.v19.firstE.length /dev/stdin 0 0 | cut -f 1-5,7 | perl $myperl gencode.v19.firstI.length /dev/stdin 0 0 | cut -f 1-6,8 | perl $myperl <(cut -f 1,5 gencode.100way.firstE.tab | sed 's/:/\t/' | cut -f 1,3) /dev/stdin 0 0 | cut -f 1-7,9 | perl $myperl <(cut -f 1,5 gencode.100way.firstI.tab) /dev/stdin 0 0 | cut -f 1-8,10 | perl $myperl <(cut -f 1,5 gencode.100way.promoter.tab) /dev/stdin 0 0 | cut -f 1-9,11 | perl $myperl gencode.cpgIslandExt.length.txt /dev/stdin 0 0 | cut -f 1-10,12 | perl $myperl gencode.simpleRepeat.length.txt /dev/stdin 0 0 | cut -f 1-11,13 | perl $myperl gencode.cage_peak_phase1and2combined.count.txt /dev/stdin 0 0 | cut -f 1-12,14 >> gencode.sequenceFeatures.txt

###### general information from sequence, no matter in which cell or tissue
head -1 gencode.sequenceFeatures.txt | cut -f 2- > a
perl $myperl gencode.sequenceFeatures.txt $list 0 0 | cut -f 4- | tail -n +2 | sed 's?/?0?g' >> a
paste $list a > $list.sequenceFeatures
rm a

head -1 gencode.vert.known.motifs.mat | cut -f2- | tr "\t" "\n" | awk -F"/" '{print "Homer."$1"/"$2}' | tr "\n" "\t" | sed 's/\t$/\n/' > a
perl $myperl gencode.vert.known.motifs.mat $list 0 0 | cut -f 4- | tail -n +2 | sed 's?/?0?g' >> a
paste $list a > $list.Homer

head -1 gencode.vert.known.motifs.mat | cut -f2- | tr "\t" "\n" | awk -F"/" '{print "Homer."$1"/"$2}' | tr "\n" "\t" | sed 's/\t$/\n/' > a
perl $myperl gencode.vert.known.motifs.mat $list 0 0 | cut -f 4- | tail -n +2 | sed 's?/?0?g' >> a
paste $list.sequenceFeatures a > $list.sequenceFeatures.Homer
rm a

###### GTRD peaks and roadmap histone marks/DNA methylation, with cell or tissue
head -1 gencode.GTRD.bin.mat | cut -f 2- | tr "\t" "\n" | awk '{print "GTRD."$0}' | tr "\n" "\t" | sed 's/\t$/\n/' > a
perl $myperl gencode.GTRD.bin.mat $list 0 0 | cut -f 4- | tail -n +2 | sed 's?/?0?g' >> a
paste $list a > $list.GTRD
rm a

head -1 gencode.histone.bin.mat | cut -f 2- | tr "\t" "\n" | awk '{print "roadmap."$0}' | tr "\n" "\t" | sed 's/\t$/\n/' > a
perl $myperl gencode.histone.bin.mat $list 0 0 | cut -f 4- | tail -n +2 | sed 's?/?0?g' >> a
paste $list a > $list.roadmap.histone
rm a

head -1 gencode.DNase.bin.mat | cut -f 2- | tr "\t" "\n" | awk '{print "DHS."$0}' | tr "\n" "\t" | sed 's/\t$/\n/' > a
perl $myperl gencode.DNase.bin.mat $list 0 0 | cut -f 4- | tail -n +2 | sed 's?/?0?g' >> a
paste $list a > $list.roadmap.DNase
rm a

head -1 gencode.meth.mean.mat | cut -f 2- | tr "\t" "\n" | awk '{print "meth."$1}' | tr "\n" "\t" | sed 's/\t$/\n/' > a
perl $myperl gencode.meth.mean.mat $list 0 0 | cut -f 4- | tail -n +2 | sed 's?/?0?g' >> a
paste $list a > $list.roadmap.meth
rm a

###### loops, with cell information
head -1 gencode.HiC-Loop.count | cut -f 3- | tr "\t" "\n" | awk '{print "HiC."$1}' | tr "\n" "\t" | sed 's/\t$/\n/' > a
perl $myperl gencode.HiC-Loop.count $list 0 0 | cut -f 5- | tail -n +2 | sed 's?/?0?g' >> a
paste $list a > $list.HiC
rm a

head -1 gencode.HiChIP-Loop.count | cut -f 2- | tr "\t" "\n" | awk '{print "HiChIP."$1}' | tr "\n" "\t" | sed 's/\t$/\n/' > a
perl $myperl gencode.HiChIP-Loop.count $list 0 0 | cut -f 4- | tail -n +2 | sed 's?/?0?g' >> a
paste $list a > $list.HiChIP
rm a

head -1 gencode.Cell_Javierre_17cells.count | cut -f 2- | tr "\t" "\n" | awk '{print "PCHiC."$1}' | tr "\n" "\t" | sed 's/\t$/\n/' > a
perl $myperl gencode.Cell_Javierre_17cells.count $list 0 0 | cut -f 4- | tail -n +2 | sed 's?/?0?g' >> a
paste $list a > $list.PCHiC
rm a

head -1 gencode.Cell_Javierre_17cells.oe.count | cut -f 2- | tr "\t" "\n" | awk '{print "PCHiC.oe."$1}' | tr "\n" "\t" | sed 's/\t$/\n/' > a
perl $myperl gencode.Cell_Javierre_17cells.oe.count $list 0 0 | cut -f 4- | tail -n +2 | sed 's?/?0?g' >> a
paste $list.PCHiC a > $list.PCHiC.oe
rm a

paste hkg.tsg.srtbyPCA.class <(cut -f 3- hkg.tsg.srtbyPCA.class.sequenceFeatures.Homer) <(cut -f 3- hkg.tsg.srtbyPCA.class.GTRD) <(cut -f 3- hkg.tsg.srtbyPCA.class.roadmap.meth) <(cut -f 3- hkg.tsg.srtbyPCA.class.roadmap.histone) <(cut -f 3- hkg.tsg.srtbyPCA.class.roadmap.DNase) <(cut -f 3- hkg.tsg.srtbyPCA.class.HiC) <(cut -f 3- hkg.tsg.srtbyPCA.class.HiChIP) <(cut -f 3- hkg.tsg.srtbyPCA.class.PCHiC) <(cut -f 3- hkg.tsg.srtbyPCA.class.PCHiC.oe) > hkg.tsg.srtbyPCA.all

# Rscript $myR $list.HiC.HiChIP.PCHiC.oe &
# Rscript $myR $list.GTRD.roadmap.meth &
# Rscript $myR $list.sequenceFeatures.cage.phastCons.Homer &
for f in hkg.tsg.srtbyPCA.class.sequenceFeatures
do
	Rscript $myboxplot $f &
	Rscript $mybarplot $f &
done
