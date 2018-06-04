#!/bin/sh
myperl=/home1/04935/shaojf/myScripts/add_any_2files_together.pl
# rclone sync /home1/04935/shaojf/stampede2/housekeeping_genes/both.pc.and.nc.genes/HKG.pmut/ mygoogle:hkg_tsg/both.pc.and.nc.genes/HKG.pmut/
awk '$2~/hkg/{print $1}' hkg.tsg.srtbyPCA.class | cut -f2 -d"|" > hkg.name

# awk '$1>20{print $5}' hkg.pMUT |cut -f 1 -d"|"
awk -v OFS="\t" '$1>5{print $2,$3,$4,$6":"$1"@"$5}' hkg.pMUT > hkg.pMUT.gt5.bed
awk -v OFS="\t" '$1==1{print $2,$3,$4,$6":"$1"@"$5}' hkg.pMUT | shuf | head -n 984 > hkg.pMUT.rand1.bed
awk -v OFS="\t" '$1==1{print $2,$3,$4,$6":"$1"@"$5}' all.cancer.pMUT.counts | shuf | head -n 984 > rand1.pMUT.bed
awk -v OFS="\t" '$1>5{print $2,$3,$4,$6":"$1"@"$5}' all.cancer.pMUT.counts | shuf | head -n 984 > randgt5.pMUT.bed
for f in hkg.pMUT.gt5.bed rand1.pMUT.bed hkg.pMUT.rand1.bed randgt5.pMUT.bed
do
	pre=`echo $f | sed 's/.bed//'`
	findMotifsGenome.pl <(awk -vOFS="\t" '{print "chr"$1,$2-20,$3+20,$4}' $f) hg19 $pre.homer.motifs -size given -mknown /home1/04935/shaojf/myTools/HOMER/data/knownTFs/vertebrates/known.motifs -p 68 1> $pre.findMotifsGenome.txt 2>&1 &
done

for f in hkg.pMUT.gt5.bed rand1.pMUT.bed hkg.pMUT.rand1.bed randgt5.pMUT.bed
do
	pre=`echo $f | sed 's/.bed//'`
	annotatePeaks.pl <(awk -vOFS="\t" '{print "chr"$1,$2-20,$3+20,$4}' $f) hg19 -m /home1/04935/shaojf/myTools/HOMER/data/knownTFs/vertebrates/known.motifs -size given -cpu 68 1> $pre.vert.motifs.txt 2> $pre.vert.motifs.log &
done

for f in *homer.motifs/knownResults.txt
do
	pre=`echo $f | sed 's?.homer.motifs/knownResults.txt??'`
	awk -F"\t" -vvar=$pre -vOFS="\t" '{print var,$1,$5,$7}' $f >> pMUT.motifs
done

grep -e "/Promoter/" -e "ETS" -e "NRF" -e "Klf9" pMUT.motifs | cut -f 1-2,4 | sed 's/%//' > pMUT.Promoter.motifs
Rscript /home1/04935/shaojf/myTools/BioinformaticsDaily/RVisualization/barplot4threecolsdata.R pMUT.Promoter.motifs

more ../appris.tss/hkg.tsg.srtbyPCA.transcript.class | grep "RPL13A"
more ../appris.tss/hkg.tsg.srtbyPCA.transcript.HiChIP-Loop.count | grep "ENSG00000142541.12|ENST00000391857.4"