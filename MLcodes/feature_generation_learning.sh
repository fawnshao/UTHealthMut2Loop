#!/bin/sh
myperl=/home1/04935/shaojf/myScripts/add_any_2files_together.pl
pathology=/work/04935/shaojf/stampede2/refs/ProteinAtlas/pathology.tsv
myR=/home1/04935/shaojf/myTools/UTHealthMut2Loop/MLcodes/feature_generation_learning.R
myboxplot=/home1/04935/shaojf/myTools/UTHealthMut2Loop/MLcodes/features.boxplot.R
mybarplot=/home1/04935/shaojf/myTools/UTHealthMut2Loop/MLcodes/features.percentage.R
mystatR=/home1/04935/shaojf/myTools/UTHealthMut2Loop/MLcodes/output.percentage.R 

list=total.type.srt
# rclone sync formated_features/ mygoogle:hkg_tsg/formated_features/
##### input files, prepared elsewhere
ln -s ../motifs/total.type
head -1 total.type > total.type.srt
tail -n +2 total.type | sort -k2 >> total.type.srt 

ln -s ../conservation/promoter.1k/GTEx.100way.tab

ln -s ../motifs/GTEx.vert.known.motifs.mat
ln -s ../v1.5/GTEx.cage_peak_phase1and2combined.count.txt
ln -s ../v1.5/GTEx.sequenceFeatures.txt
ln -s ../v1.5/gencode.Cell_Javierre_17cells.count 
ln -s ../v1.5/gencode.Cell_Javierre_17cells.oe.count 
ln -s ../histone.marks/GTEx.roadmap.count.mat
ln -s ../dnamethylation/GTEx.meth.mean.mat 
ln -s ../motifs/GTEx.GTRD.mat

cut -f 1-2,4- ../motifs/total.type.HiC.HiChIP > total.type.HiC.HiChIP

###### general information from sequence, no matter in which cell or tissue
head -1 GTEx.sequenceFeatures.txt | cut -f 2- > a
perl $myperl GTEx.sequenceFeatures.txt $list 0 0 | cut -f 4- | tail -n +2 | sed 's?/?0?g' >> a
paste $list a > $list.sequenceFeatures
rm a

echo CAGE.count > a
perl $myperl GTEx.cage_peak_phase1and2combined.count.txt $list 0 0 | cut -f 4 | tail -n +2 | sed 's?/?0?g' >> a
paste $list.sequenceFeatures a > $list.sequenceFeatures.cage
rm a

### 5th: mean0 - average over bases with non-covered bases counting as zeroes
### 6th: mean - average over just covered bases
echo promoter.mean.100way.phastCons > a
perl $myperl <(cut -f 1,5 GTEx.100way.tab) $list 0 0 | cut -f 4 | tail -n +2 | sed 's?/?0?g' >> a
paste $list.sequenceFeatures.cage a > $list.sequenceFeatures.cage.phastCons
rm a

head -1 GTEx.vert.known.motifs.mat | cut -f2- | tr "\t" "\n" | awk -F"/" '{print "Homer."$1"/"$2}' | tr "\n" "\t" | sed 's/\t$/\n/' > a
perl $myperl GTEx.vert.known.motifs.mat $list 0 0 | cut -f 4- | tail -n +2 | sed 's?/?0?g' >> a
paste $list.sequenceFeatures.cage.phastCons a > $list.sequenceFeatures.cage.phastCons.Homer
rm a


###### GTRD peaks and roadmap histone marks/DNA methylation, with cell or tissue
head -1 GTEx.GTRD.mat | cut -f2- | tr "\t" "\n" | awk '{print "GTRD."$0}' | tr "\n" "\t" | sed 's/\t$/\n/' > a
perl $myperl GTEx.GTRD.mat $list 0 0 | cut -f 4- | tail -n +2 | sed 's?/?0?g' >> a
paste $list a > $list.GTRD
rm a

head -1 GTEx.roadmap.count.mat | cut -f2- | tr "\t" "\n" | awk '{print "roadmap."$0}' | tr "\n" "\t" | sed 's/\t$/\n/' > a
perl $myperl GTEx.roadmap.count.mat $list 0 0 | cut -f 4- | tail -n +2 | sed 's?/?0?g' >> a
paste $list.GTRD a > $list.GTRD.roadmap
rm a

head -1 GTEx.meth.mean.mat | cut -f 2- | tr "\t" "\n" | awk '{print "meth."$1}' | tr "\n" "\t" | sed 's/\t$/\n/' > a
perl $myperl GTEx.meth.mean.mat $list 0 0 | cut -f 4- | tail -n +2 | sed 's?/?0?g' >> a
paste $list.GTRD.roadmap a > $list.GTRD.roadmap.meth
rm a

###### loops, with cell information
head -1 total.type.HiC.HiChIP | cut -f 3- > a
perl $myperl total.type.HiC.HiChIP $list 0 0 | cut -f 5- | tail -n +2 | sed 's?/?0?g' >> a
paste $list a > $list.HiC.HiChIP
rm a

head -1 gencode.Cell_Javierre_17cells.count | cut -f 2- | tr "\t" "\n" | awk '{print "PCHiC."$1}' | tr "\n" "\t" | sed 's/\t$/\n/' > a
perl $myperl gencode.Cell_Javierre_17cells.count $list 0 0 | cut -f 4- | tail -n +2 | sed 's?/?0?g' >> a
paste $list.HiC.HiChIP a > $list.HiC.HiChIP.PCHiC
rm a

head -1 gencode.Cell_Javierre_17cells.oe.count | cut -f 2- | tr "\t" "\n" | awk '{print "PCHiC.oe."$1}' | tr "\n" "\t" | sed 's/\t$/\n/' > a
perl $myperl gencode.Cell_Javierre_17cells.oe.count $list 0 0 | cut -f 4- | tail -n +2 | sed 's?/?0?g' >> a
paste $list.HiC.HiChIP.PCHiC a > $list.HiC.HiChIP.PCHiC.oe
rm a

Rscript $myR $list.HiC.HiChIP.PCHiC.oe &
Rscript $myR $list.GTRD.roadmap.meth &
Rscript $myR $list.sequenceFeatures.cage.phastCons.Homer &


##########################
ln -s ../GM_Liver_as_control/pc.hkg.subsettsg.genes
subs=pc.hkg.subsettsg.genes
for post in HiC.HiChIP.PCHiC.oe GTRD.roadmap.meth sequenceFeatures.cage.phastCons.Homer
do
	head -1 $list.$post > $subs.$post
	perl $myperl $list.$post $subs 0 0 | cut -f 1-2,5- >> $subs.$post
done

Rscript $myR $subs.HiC.HiChIP.PCHiC.oe &
Rscript $myR $subs.GTRD.roadmap.meth &
Rscript $myR $subs.sequenceFeatures.cage.phastCons.Homer &

for post in HiC.HiChIP.PCHiC.oe GTRD.roadmap.meth sequenceFeatures.cage.phastCons.Homer
do
	Rscript $myboxplot pc.hkg.subsettsg.genes.$post &
	Rscript $myboxplot total.type.srt.$post &
done

for post in HiC.HiChIP.PCHiC.oe GTRD.roadmap.meth sequenceFeatures.cage.phastCons.Homer
do
	Rscript $mybarplot pc.hkg.subsettsg.genes.$post &
	Rscript $mybarplot total.type.srt.$post &
done
for post in HiC.HiChIP.PCHiC.oe GTRD.roadmap.meth sequenceFeatures.cage.phastCons.Homer
do
	Rscript $mystatR total.type.srt.$post &
done

# for post in sequenceFeatures.cage.phastCons.Homer GTRD.roadmap.meth
# do
# 	# sed -i '1s?/?.?g' pc.hkg.subsettsg.genes.$post
# 	# sed -i '1s?/?.?g' total.type.srt.$post
# 	Rscript $myboxplot pc.hkg.subsettsg.genes.$post &
# 	Rscript $myboxplot total.type.srt.$post &
# done




### use the following order to concatenate all the fetures together
list=total.type
perl $myperl <(cut -f 4 housekeepinggene.byGTEx.gene.bed) <(cut -f 4 GTEx.proteincoding.noMT.bed) 0 0 | perl $myperl <(cut -f 4 tissuespecificgene.byGTEx.gene.bed) /dev/stdin 0 0 | tail -n +2 | awk -vOFS="\t" 'BEGIN{print "Gene","type"}{a="other";if($2!="/"){a="hkg";}if($3!="/"){a="tsg";}print $1,a}' > $list
echo "HiC.CH12-LX HiC.GM12878_primary+replicate HiC.HeLa HiC.HMEC HiC.HUVEC HiC.IMR90 HiC.K562 HiC.KBM7 HiC.NHEK" | tr " " "\t" > a
perl $myperl <(cut -f 1-2,4,6- gencode.HiC-Loop.count) $list 0 0 | tail -n +2 | cut -f 4- | sed 's?/?0?g' >> a
paste $list a > $list.HiC
rm a
echo "HiChIP.GM12878 HiChIP.K562 HiChIP.Naive HiChIP.Th17 HiChIP.Treg" | tr " " "\t" > a
perl $myperl gencode.HiChIP-Loop.count $list 0 0 | tail -n +2 | cut -f 4- | sed 's?/?0?g' >> a
paste $list.HiC a > $list.HiC.HiChIP
rm a
head -1 GTEx.GTRD.mat | cut -f2- | tr "\t" "\n" | awk '{print "GTRD.maxcelllinecount."$0}' | tr "\n" "\t" | sed 's/\t$/\n/' > a
perl $myperl GTEx.GTRD.mat $list 0 0 | tail -n +2 | cut -f 4- | sed 's?/?0?g' >> a
paste $list.HiC.HiChIP a > $list.HiC.HiChIP.GTRD
rm a
head -1 GTEx.known.motifs.mat | cut -f2- | tr "\t" "\n" | awk -F "/" '{print "homer.motifs."$1"."$2}' | tr "\n" "\t" | sed 's/\t$/\n/' > a
perl $myperl GTEx.known.motifs.mat $list 0 0 | tail -n +2 | cut -f 4- | sed 's?/?0?g' >> a
paste $list.HiC.HiChIP.GTRD a > $list.HiC.HiChIP.GTRD.homer
rm a
# GTEx.cage_peak_phase1and2combined.count.txt
echo CAGE.count > a
perl $myperl GTEx.cage_peak_phase1and2combined.count.txt $list 0 0 | cut -f 4 | tail -n +2 | sed 's?/?0?g' >> a
paste $list.HiC.HiChIP.GTRD.homer a > $list.HiC.HiChIP.GTRD.homer.cage
rm a
head -1 GTEx.sequenceFeatures.txt | cut -f 2- > a
perl $myperl GTEx.sequenceFeatures.txt $list 0 0 | cut -f 4- | tail -n +2 | sed 's?/?0?g' >> a
paste $list.HiC.HiChIP.GTRD.homer.cage a > $list.HiC.HiChIP.GTRD.homer.cage.sequenceFeatures
rm a
head -1 GTEx.roadmap.count.mat | cut -f 2- > a
perl $myperl GTEx.roadmap.count.mat $list 0 0 | cut -f 4- | tail -n +2 | sed 's?/?0?g' >> a
paste $list.HiC.HiChIP.GTRD.homer.cage.sequenceFeatures a > $list.HiC.HiChIP.GTRD.homer.cage.sequenceFeatures.roadmap
# paste $list a > $list.roadmap
rm a
