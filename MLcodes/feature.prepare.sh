#!/bin/sh
myperl=/home1/04935/shaojf/myScripts/add_any_2files_together.pl
pathology=/work/04935/shaojf/stampede2/refs/ProteinAtlas/pathology.tsv
myR=/home1/04935/shaojf/myTools/UTHealthMut2Loop/MLcodes/feature_generation_learning.R
myboxplot=/home1/04935/shaojf/myTools/UTHealthMut2Loop/MLcodes/features.boxplot.R
mybarplot=/home1/04935/shaojf/myTools/UTHealthMut2Loop/MLcodes/features.percentage.R
mystatR=/home1/04935/shaojf/myTools/UTHealthMut2Loop/MLcodes/output.percentage.R 

list=hkg.tsg.srtbyPCA.class
# rclone sync -L /home1/04935/shaojf/stampede2/housekeeping_genes/both.pc.and.nc.genes/ mygoogle:hkg_tsg/both.pc.and.nc.genes/
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

paste hkg.tsg.srtbyPCA.class <(cut -f 3- hkg.tsg.srtbyPCA.class.sequenceFeatures.Homer) <(cut -f 3- hkg.tsg.srtbyPCA.class.GTRD) <(cut -f 3- hkg.tsg.srtbyPCA.class.roadmap.meth) <(cut -f 3- hkg.tsg.srtbyPCA.class.roadmap.histone) <(cut -f 3- hkg.tsg.srtbyPCA.class.roadmap.DNase) <(cut -f 3- hkg.tsg.srtbyPCA.class.HiC) <(cut -f 3- hkg.tsg.srtbyPCA.class.HiChIP) <(cut -f 3- hkg.tsg.srtbyPCA.class.PCHiC.oe) > hkg.tsg.srtbyPCA.all

#########################
# for additional test
head -1 gencode.v19.enhanceratlas.mat  | cut -f 2- | tr "\t" "\n" | awk '{print "enhanceratlas."$1}' | tr "\n" "\t" | sed 's/\t$/\n/' > a
perl $myperl gencode.v19.enhanceratlas.mat $list 0 0 | cut -f 4- | tail -n +2 | sed 's?/?0?g' >> a
paste $list a > $list.enhanceratlas
rm a

paste hkg.tsg.srtbyPCA.class <(cut -f 3- hkg.tsg.srtbyPCA.class.sequenceFeatures) <(cut -f 3- hkg.tsg.srtbyPCA.class.HiC) <(cut -f 3- hkg.tsg.srtbyPCA.class.HiChIP) <(cut -f 3- hkg.tsg.srtbyPCA.class.PCHiC.oe) <(cut -f 3- hkg.tsg.srtbyPCA.class.enhanceratlas) > hkg.tsg.srtbyPCA.class.part

head -1 gencode.hESC.IMR90.count | cut -f 2- | tr "\t" "\n" | awk '{print "TAD."$1}' | tr "\n" "\t" | sed 's/\t$/\n/' > a
perl $myperl gencode.hESC.IMR90.count $list 0 0 | cut -f 4- | tail -n +2 | sed 's?/?0?g' >> a
paste $list a > $list.TAD
rm a

paste hkg.tsg.srtbyPCA.class <(cut -f 3- hkg.tsg.srtbyPCA.class.sequenceFeatures) <(cut -f 3- hkg.tsg.srtbyPCA.class.HiC) <(cut -f 3- hkg.tsg.srtbyPCA.class.HiChIP) <(cut -f 3- hkg.tsg.srtbyPCA.class.PCHiC.oe) <(cut -f 3- hkg.tsg.srtbyPCA.class.enhanceratlas) <(cut -f 3- $list.TAD) > $list.part2

#########################
perl ~/myScripts/add_any_2files_together.pl gencode.neighbors.txt hkg.tsg.srtbyPCA.class 3 0 > hkg.tsg.srtbyPCA.class.neighbors
cut -f 1-2,8,12,14,15 hkg.tsg.srtbyPCA.class.neighbors > hkg.tsg.srtbyPCA.class.neighbors.sim

perl ~/myScripts/add_any_2files_together.pl gencode.tss.neighbors.txt hkg.tsg.srtbyPCA.class 3 0 > hkg.tsg.srtbyPCA.class.tss.neighbors
cut -f 1-2,8,12,14,15 hkg.tsg.srtbyPCA.class.tss.neighbors > hkg.tsg.srtbyPCA.class.tss.neighbors.sim

#########################
# add enhanceratlas
paste hkg.tsg.srtbyPCA.all <(cut -f 3- hkg.tsg.srtbyPCA.class.enhanceratlas) > hkg.tsg.srtbyPCA.all2
cat hkg.tsg.srtbyPCA.all.dnase.FeatureSelection.union_feat.tsv hkg.tsg.srtbyPCA.all.pchic.FeatureSelection.union_feat.tsv hkg.tsg.srtbyPCA.all.hic.FeatureSelection.union_feat.tsv hkg.tsg.srtbyPCA.all.seqfeature.FeatureSelection.union_feat.tsv hkg.tsg.srtbyPCA.all.homermotif.FeatureSelection.union_feat.tsv hkg.tsg.srtbyPCA.all.meth.FeatureSelection.union_feat.tsv hkg.tsg.srtbyPCA.all.gtrdpeaks.FeatureSelection.union_feat.tsv hkg.tsg.srtbyPCA.all.histone.FeatureSelection.union_feat.tsv hkg.tsg.srtbyPCA.class.enhanceratlas.enhanceratlas.FeatureSelection.union_feat.tsv | awk '$2>0.8' > important.features.tsv
cat hkg.tsg.srtbyPCA.all.dnase.sim.FeatureSelection.union_feat.tsv hkg.tsg.srtbyPCA.all.pchic.sim.FeatureSelection.union_feat.tsv hkg.tsg.srtbyPCA.all.hic.sim.FeatureSelection.union_feat.tsv hkg.tsg.srtbyPCA.all.seqfeature01.sim.FeatureSelection.union_feat.tsv hkg.tsg.srtbyPCA.all.homermotif.sim.FeatureSelection.union_feat.tsv hkg.tsg.srtbyPCA.all.meth.sim.FeatureSelection.union_feat.tsv hkg.tsg.srtbyPCA.all.gtrdpeaks.sim.FeatureSelection.union_feat.tsv hkg.tsg.srtbyPCA.all.histone.sim.FeatureSelection.union_feat.tsv hkg.tsg.srtbyPCA.class.enhanceratlas.enhanceratlas.sim.FeatureSelection.union_feat.tsv | awk '$2>0.8' > sim.important.features.tsv

awk '$5 > 0.5 && $5 - $7 > 0.3' hkg.tsg.srtbyPCA.percentage.tsv | grep -v "Homer" > bigdiff.features.tsv
awk '$5 - $7 > 0.1' hkg.tsg.srtbyPCA.percentage.tsv | grep "Homer" >> bigdiff.features.tsv
awk '$7 > 0.5 && $7 - $5 > 0.3' hkg.tsg.srtbyPCA.percentage.tsv >> bigdiff.features.tsv



###### try motif again
perl $myperl gencode.v19.gene.bed <(awk '$2=="hkg1"' $list) 3 0 | cut -f 3- > hkg1.gene.bed
perl $myperl gencode.v19.gene.bed <(awk '$2=="hkg2"' $list) 3 0 | cut -f 3- > hkg2.gene.bed
perl $myperl gencode.v19.gene.bed <(awk '$2=="hkg3"' $list) 3 0 | cut -f 3- > hkg3.gene.bed
perl $myperl gencode.v19.gene.bed <(awk '$2=="hkg4"' $list) 3 0 | cut -f 3- > hkg4.gene.bed
perl $myperl gencode.v19.gene.bed <(awk -F"\t" '$2!~/hkg/ && $2!="mixTSG"' $list | tail -n +2) 3 0 | cut -f 3- > singleTSG.gene.bed
perl $myperl gencode.v19.gene.bed <(awk '$2=="mixTSG"' $list) 3 0 | cut -f 3- > mixTSG.gene.bed
cat hkg1.gene.bed hkg2.gene.bed hkg3.gene.bed hkg4.gene.bed > hkg.gene.bed

for gene in hkg1.gene.bed  hkg2.gene.bed  hkg3.gene.bed  hkg4.gene.bed  hkg.gene.bed  mixTSG.gene.bed  singleTSG.gene.bed
do
	pre=`echo $gene | awk -F"." '{print $1}'`
	findMotifsGenome.pl <(awk -vOFS="\t" '{a=$2-1000;b=$2+1000;if($6=="-"){a=$3-1000;b=$3+1000}if(a<0){a=0;}print $1,a,b,$4,"1000",$6}' $gene) hg19 $pre.homer.motifs -size given -p 68 -mknown /home1/04935/shaojf/myTools/HOMER/data/knownTFs/all.motifs > $pre.findMotifsGenome.txt &
done

# cat *.homer.motifs/knownResults.txt | cut -f 1-2,5,7 | more
cat *.homer.motifs/knownResults.txt | awk -F"\t" '$5<0.001 && $7-$9>5{print $1"\t"$2}' | sort | uniq | grep -v "Motif Name" > known.motif.id
echo "motif seq hkg1.q hkg1.ratio hkg2.q hkg2.ratio hkg3.q hkg3.ratio hkg4.q hkg4.ratio hkg.q hkg.ratio singleTSG.q singleTSG.ratio mixTSG.q mixTSG.ratio" | tr " " "\t" > known.motif.mat
perl $myperl <(awk -F"\t" '$5<0.001 && $7-$9>5' hkg1.homer.motifs/knownResults.txt | cut -f 1-2,5,7) known.motif.id 0 0 | cut -f 1-2,5- | perl $myperl <(awk -F"\t" '$5<0.001 && $7-$9>5' hkg2.homer.motifs/knownResults.txt | cut -f 1-2,5,7) /dev/stdin 0 0 | cut -f 1-4,7- | perl $myperl <(awk -F"\t" '$5<0.001 && $7-$9>5' hkg3.homer.motifs/knownResults.txt | cut -f 1-2,5,7) /dev/stdin 0 0 | cut -f 1-6,9- | perl $myperl <(awk -F"\t" '$5<0.001 && $7-$9>5' hkg4.homer.motifs/knownResults.txt | cut -f 1-2,5,7) /dev/stdin 0 0 | cut -f 1-8,11- | perl $myperl <(awk -F"\t" '$5<0.001 && $7-$9>5' hkg.homer.motifs/knownResults.txt | cut -f 1-2,5,7) /dev/stdin 0 0 | cut -f 1-10,13- | perl $myperl <(awk -F"\t" '$5<0.001 && $7-$9>5' singleTSG.homer.motifs/knownResults.txt | cut -f 1-2,5,7) /dev/stdin 0 0 | cut -f 1-12,15- | perl $myperl <(awk -F"\t" '$5<0.001 && $7-$9>5' mixTSG.homer.motifs/knownResults.txt | cut -f 1-2,5,7) /dev/stdin 0 0 | cut -f 1-14,17- | sed 's/%//g' >> known.motif.mat

awk -F"\t" '$5<0.001 && $7-$9>10{print $1"\t"$2}' hkg.homer.motifs/knownResults.txt | grep -v "Motif Name" > sim.known.motif.id
echo "motif seq hkg.q hkg.ratio singleTSG.q singleTSG.ratio mixTSG.q mixTSG.ratio" | tr " " "\t" > sim.known.motif.mat
perl $myperl <(cut -f 1-2,5,7 hkg.homer.motifs/knownResults.txt) sim.known.motif.id 0 0 | cut -f 1-2,5- | perl $myperl <(cut -f 1-2,5,7 singleTSG.homer.motifs/knownResults.txt) /dev/stdin 0 0 | cut -f 1-4,7- | perl $myperl <(cut -f 1-2,5,7 mixTSG.homer.motifs/knownResults.txt) /dev/stdin 0 0 | cut -f 1-6,9- | sed 's/%//g' >> sim.known.motif.mat

######################
cut -f 1-5,87,89,90,91,92,93,96,100,101,102,103,104,105,106,107,108,109,128,275,276,277,313,314,311,354 gencode.vert.known.motifs.txt > gencode.vert.known.ETS.Sp1.YY1.txt
perl ~/myTools/UTHealthMut2Loop/relatedScripts/count.motif.bin.from.homer.pl <(sed 's/),/);/g' gencode.vert.known.ETS.Sp1.YY1.txt) > gencode.vert.known.ETS.Sp1.YY1.bin.txt
perl ~/myTools/UTHealthMut2Loop/relatedScripts/count.motif.bin.from.homer.pl <(perl $myperl gencode.vert.known.ETS.Sp1.YY1.txt hkg.tsg.srtbyPCA.class 0 0 | awk '$2~/hkg/ || NR==1' | cut -f 3- | sed 's/),/);/g') > hkg.vert.known.ETS.Sp1.YY1.bin.txt
perl ~/myTools/UTHealthMut2Loop/relatedScripts/count.motif.bin.from.homer.pl <(perl $myperl gencode.vert.known.ETS.Sp1.YY1.txt hkg.tsg.srtbyPCA.class 0 0 | awk '$2~/mixTSG/ || NR==1' | cut -f 3- | sed 's/),/);/g') > mixTSG.vert.known.ETS.Sp1.YY1.bin.txt
perl ~/myTools/UTHealthMut2Loop/relatedScripts/count.motif.bin.from.homer.pl <(perl $myperl gencode.vert.known.ETS.Sp1.YY1.txt hkg.tsg.srtbyPCA.class 0 0 | awk '($2!~/singleTSG/ && $2!~/mixTSG/) || NR==1' | cut -f 3- | sed 's/),/);/g') > singleTSG.vert.known.ETS.Sp1.YY1.bin.txt
