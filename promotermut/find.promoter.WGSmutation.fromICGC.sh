#!/bin/sh
# rclone sync /home1/04935/shaojf/stampede2/mutations/promotermut/ mygoogle:NCmutation/promotermut/
tssBED=hg19.refGene.tss.uniq.srt.bed
promoterLEN=1000
for mutationTSV in simple_somatic_mutation.open.*.tsv.gz
do
	outpre=`echo $mutationTSV | sed 's/simple_somatic_mutation.open.//;s/.tsv.gz//;'`
	gunzip -c $mutationTSV | awk -F"\t" -vOFS="\t" '{print $9,$10-1,$11,$4"|"$3"|"$13"|"$14"|"$16"|"$17"|"$34}' | tail -n +2 | uniq | bedtools sort -i - > $mutationTSV.srt.bed &
done
wait
# for mutations in simple_somatic_mutation.open.*.tsv.gz.srt.bed
for mutations in `ls simple_somatic_mutation.open.*.tsv.gz.srt.bed | grep -v "MELA-AU"`
do
	outpre=`echo $mutations | sed 's/simple_somatic_mutation.open.//;s/.tsv.gz.srt.bed//;'`
	bedtools closest -D b -a <(uniq $mutations) -b $tssBED | awk -F "\t" -v var=$promoterLEN -v OFS="\t" '$(NF-6)!="." && $NF > -1 * var && $NF < var {print $1,$2,$3,$4,$8,$11}' > ${outpre}.pMUT &
done
wait
split -l 1000000 simple_somatic_mutation.open.MELA-AU.tsv.gz.srt.bed 
for mutations in x??
do
	outpre=$mutations
	bedtools closest -D b -a <(uniq $mutations) -b $tssBED | awk -F "\t" -v var=$promoterLEN -v OFS="\t" '$(NF-6)!="." && $NF > -1 * var && $NF < var {print $1,$2,$3,$4,$8,$11}' > ${outpre}.pMUT &
done
cat xa?.pMUT > MELA-AU.pMUT
rm xa*

cat *.pMUT | cut -f 1-3,5-6 | sort | uniq -c | sort -k1,1nr > all.cancer.pMUT.counts
# awk -v OFS="\t" '$1 >= 5 && $6 < 100 && $6 > -200{print $2,$3,$4,$5,$6,$1}' all.cancer.pMUT.counts | grep -v "MIR" | bedtools sort -i > all.cancer.recurrentpmut.bed
awk -v OFS="\t" '$1>10 && $5!~/MIR/ {print $2,$3,$4,$6"@"$5"#"$1}' all.cancer.pMUT.counts > all.cancer.toprecurrentpmut.bed
awk -v OFS="\t" '$1>5 && $5!~/MIR/ {print $2,$3,$4,$6"@"$5"#"$1}' all.cancer.pMUT.counts > all.cancer.highrecurrentpmut.bed

awk -v OFS="\t" '$1>1 && $5!~/MIR/ {print $2,$3,$4,$6"@"$5"#"$1}' all.cancer.pMUT.counts > all.cancer.recurrentpmut.bed
awk -v OFS="\t" '$1==1 && $5!~/MIR/ {print $2,$3,$4,$6"@"$5"#"$1}' all.cancer.pMUT.counts > all.cancer.nonrecurrentpmut.bed
for gene in all.cancer.*recurrentpmut.bed
do
	pre=`echo $gene | sed 's/.bed//'`
	findMotifsGenome.pl <(awk -vOFS="\t" '{print "chr"$1,$2-20,$3+20,$4}' $gene) hg19 $pre.homer.motifs -size given -mknown /home1/04935/shaojf/myTools/HOMER/data/knownTFs/vertebrates/known.motifs -p 60 1> $pre.findMotifsGenome.txt 2>&1 &
	annotatePeaks.pl <(awk -vOFS="\t" '{print "chr"$1,$2-20,$3+20,$4}' $gene) hg19 -m /home1/04935/shaojf/myTools/HOMER/data/knownTFs/vertebrates/known.motifs -size given -cpu 60 1> $pre.promoter.vert.motifs.txt 2> $pre.promoter.vert.motifs.log &
done

for gene in all.cancer.*recurrentpmut.bed
do
	pre=`echo $gene | sed 's/.bed//'`
	grep -wf Cosmic.CancerGeneCensus.oncogene.gene $gene > $pre.oncogene.bed &
	grep -vwf Cosmic.CancerGeneCensus.oncogene.gene $gene > $pre.other.bed &
done

for gene in all.cancer.*recurrentpmut.bed
do
	pre=`echo $gene | sed 's/.bed//'`
	bedtools closest -D b -a <(bedtools sort -i $gene) -b hg19.refGene.tss.Cosmic.CancerGeneCensus.oncogene.gene.bed > $pre.2oncogene.txt &
done

for gene in all.cancer.*recurrentpmut.other.bed
do
	pre=`echo $gene | sed 's/.bed//'`
	bedtools closest -D b -a <(bedtools sort -i $gene) -b hg19.refGene.tss.Cosmic.CancerGeneCensus.oncogene.gene.bed > $pre.2oncogene.txt &
done

for f in all.cancer.*recurrentpmut.other.2oncogene.txt
do
	pre=`echo $f | sed 's/.2oncogene.txt//'`
	awk  -vOFS="\t" '$5!="." && $11 > -200000 && $11 < 200000{print $1,$2,$3,$4}' $f > $pre.onconeighbor.bed
	awk  -vOFS="\t" '$11 < -20000000 || $11 > 20000000{print $1,$2,$3,$4}' $f > $pre.oncodistal.bed
	findMotifsGenome.pl <(awk -vOFS="\t" '{print "chr"$1,$2-20,$3+20,$4}' $pre.onconeighbor.bed) hg19 $pre.onconeighbor.homer.motifs -size given -mknown /home1/04935/shaojf/myTools/HOMER/data/knownTFs/vertebrates/known.motifs -p 60 1> $pre.onconeighbor.findMotifsGenome.txt 2>&1 &
	findMotifsGenome.pl <(awk -vOFS="\t" '{print "chr"$1,$2-20,$3+20,$4}' $pre.oncodistal.bed) hg19 $pre.oncodistal.homer.motifs -size given -mknown /home1/04935/shaojf/myTools/HOMER/data/knownTFs/vertebrates/known.motifs -p 60 1> $pre.oncodistal.findMotifsGenome.txt 2>&1 &
done
# non recurrent pmut might due to limited sequenced samples, no big difference to 2 or 3.
# so compare recurrent pmut to oncogene and non-oncogene

































#######################################################################
tssBED=hg19.refGene.tss.uniq.srt.bed
promoterLEN=1000
for mutationTSV in simple_somatic_mutation.open.*.tsv.gz
do
	outpre=`echo $mutationTSV | sed 's/simple_somatic_mutation.open.//;s/.tsv.gz//;'`
	gunzip -c $mutationTSV | awk -F"\t" -vOFS="\t" '$34=="WGS"{print $9,$10-1,$11,$2,".","+"}' | uniq | bedtools sort -i - > $mutationTSV.WGS.srt.bed &
done

# ls -lS *.WGS.srt.bed | awk '$5>0' | grep -v "MELA-AU" | awk '{print $9}' | xargs -n 1 -I old mv old selected.cancers/
for mutations in simple_somatic_mutation.open.*.tsv.gz.WGS.srt.bed
do
	outpre=`echo $mutations | sed 's/simple_somatic_mutation.open.//;s/.tsv.gz.WGS.srt.bed//;'`
	bedtools closest -D b -a <(uniq $mutations) -b $tssBED | awk -F "\t" -v var=$promoterLEN -v OFS="\t" '$(NF-6)!="." && $NF > -1 * var && $NF < var {print $1,$2,$3,$4,$10,$13}' > ${outpre}.pMUT &
done

for pre in `cut -d " " -f2 mut2promoter/selected.cancers/all.cancer.wgs1M.cancertypes`
do
	awk -F"\t" -v var=$pre -vOFS="\t" '{$4=var"|"$4;print $0}' mut2promoter/selected.cancers/tss.1k.mut/${pre}.pMUT >> 13cancers.pMUT
done
for pre in `cut -d " " -f2 mut2promoter/selected.cancers/all.cancer.wgs1M.cancertypes`
do
	samplecount=`cut -f 4 mut2promoter/selected.cancers/simple_somatic_mutation.open.$pre.tsv.gz.WGS.srt.bed | sort | uniq | wc -l`
	echo $pre $samplecount
done
# BRCA-EU 569
# ESAD-UK 253
# LICA-FR 39
# LIRI-JP 258
# LMS-FR 67
# LUSC-KR 30
# MALY-DE 241
# MELA-AU 183
# OV-AU 93
# PACA-AU 252
# PACA-CA 209
# PBCA-DE 457
# SKCA-BR 80
for pre in `cut -d " " -f2 mut2promoter/selected.cancers/all.cancer.wgs1M.cancertypes`
do
	mutcount=`cat mut2promoter/selected.cancers/simple_somatic_mutation.open.$pre.tsv.gz.WGS.srt.bed | wc -l`
	echo $pre $mutcount
done
# BRCA-EU	3973947
# ESAD-UK	8104957
# LICA-FR	1246616
# LIRI-JP	3883294
# LMS-FR	3668778
# LUSC-KR	1717259
# MALY-DE	3274936
# MELA-AU	23536812
# OV-AU	1070420
# PACA-AU	1366211
# PACA-CA	1969896
# PBCA-DE	1143929
# SKCA-BR	3445643


####################################################################################
############################## get the mutation information ##############################
for pre in `cut -d " " -f2 mut2promoter/selected.cancers/all.cancer.wgs1M.cancertypes`
do
	rawmutation=mut2promoter/simple_somatic_mutation.open.$pre.tsv.gz
	gunzip -c $rawmutation | awk -F"\t" -vOFS="\t" '$34=="WGS"' | cut -f 2-3,9-11,14,16-17 | uniq > $pre.WGS.rawmut.txt
done

# grep -w -e "PAX5" -e "TERT" oncogene.200kb.nononcogene.neighbors.pMUT.sim > ../ICGC.WGS.mut2tss.vis/13cancers.2oncogene.nononcogene.neighbors.sim
paste <(cut -f 1-4 13cancers.2oncogene.nononcogene.neighbors.sim | sed 's/|/\t/') <(cut -f 5 13cancers.2oncogene.nononcogene.neighbors.sim | cut -f 1 -d"|") <(cut -f 6 13cancers.2oncogene.nononcogene.neighbors.sim | cut -f 1 -d"|") | sort | uniq > 13cancers.2oncogene.sim
grep -w -e "TERT" -e "PAX5" ../ICGC.WGS.mut2oncogene/13cancers.pMUT > 13cancers.2oncogene.pMUT
paste <(cut -f 1-4 13cancers.2oncogene.pMUT | sed 's/|/\t/') <(cut -f 5 13cancers.2oncogene.pMUT | cut -f 1 -d"|") | sort | uniq > 13cancers.2oncogene.pMUT.sim

for pre in `cut -f 4 13cancers.2oncogene.sim | sort | uniq`
do
	perl ~/myScripts/add_any_2files_together.pl <(awk -F"\t" -vOFS="\t" '{print $2":"$1":"$3":"$5,$0}' $pre.WGS.rawmut.txt) <(awk -F"\t" -vOFS="\t" '{print $4":"$5":"$1":"$3,$6,$7}' 13cancers.2oncogene.sim | grep -w $pre) 0 0 >> 13cancers.2oncogene.sim.tmp
done
for pre in `cut -f 4 13cancers.2oncogene.sim | sort | uniq`
do
	perl ~/myScripts/add_any_2files_together.pl <(awk -F"\t" -vOFS="\t" '{print $2":"$1":"$3":"$5,$0}' $pre.WGS.rawmut.txt) <(awk -F"\t" -vOFS="\t" '{print $4":"$5":"$1":"$3,$6,"Oncogene"}' 13cancers.2oncogene.pMUT.sim | grep -w $pre) 0 0 >> 13cancers.2oncogene.sim.tmp.1
done
echo "Hugo_Symbol Chromosome Start_Position End_Position Variant_Classification Variant_Type Tumor_Sample_Barcode Gene_Type Reference_Allele Tumor_Seq_Allele1 Tumor_Seq_Allele2 Cancer_Type" | tr " " "\t" > 13cancers.2oncogene.maf
cat 13cancers.2oncogene.sim.tmp 13cancers.2oncogene.sim.tmp.1 | cut -f 2-3,5- | awk -F"\t" -vOFS="\t" '{print $1,$5,$6,$7,$8,$8,$3,$2,$9,$9,$10,$4}' >> 13cancers.2oncogene.maf


####################################################################################
# rclone sync /work/04935/shaojf/stampede2/mutations/wenbo.paper.check/ICGC.WGS.mut2oncogene/ mygoogle:NCmutation/mutations/wenbo.paper.check/ICGC.WGS.mut2oncogene/
############ get the neighboring promoter mutations for random 279 nononcogene #############
# get_seeded_random()
# {
#   seed="$1";
#   openssl enc -aes-256-ctr -pass pass:"$seed" -nosalt </dev/zero 2>/dev/null;
# }
# shuf --random-source=<(get_seeded_random 123456) <(seq 1 10)
# shuf --random-source=<(echo 123456) <(seq 1 10)
# shuf --random-source=<(echo 123456) hg19.refGene.tss.uniq.srt.bed | grep "NM" | cut -f 4 | cut -f 1 -d "|" |  grep -wvf Cosmic.CancerGeneCensus.oncogene.gene | head -n 273 > seed123456.rand273.gene

# echo 123456 | shuf --random-source=/dev/stdin hg19.refGene.tss.uniq.srt.bed | grep "NM" | cut -f 4 | cut -f 1 -d "|" | grep -wvf Cosmic.CancerGeneCensus.oncogene.gene | head -n 278 | sort | uniq | wc -l
echo 123456 | shuf --random-source=/dev/stdin hg19.refGene.tss.uniq.srt.bed | grep "NM" | cut -f 4 | cut -f 1 -d "|" | grep -wvf Cosmic.CancerGeneCensus.oncogene.gene | head -n 278 | sort | uniq > seed123456.rand273.gene
grep -wf seed123456.rand273.gene hg19.refGene.tss.uniq.srt.bed | grep -v "\-AS" > hg19.refGene.tss.rand273.gene.bed
bedtools closest -D b -io -a hg19.refGene.tss.uniq.srt.bed -b hg19.refGene.tss.rand273.gene.bed | awk '$7!="." && (($13 > -200000 && $13 < -1000) || ($13 > 1000 && $13 < 200000))' > rand273.200kb.neighbors.txt
perl ~/myScripts/add_any_2files_together.pl <(cut -f 4 hg19.refGene.tss.rand273.gene.bed) rand273.200kb.neighbors.txt 0 3 | awk '$14!="/"' | cut -f 1-13 | uniq > rand273.200kb.pseudoncogene.neighbors.txt
##### remove miR
perl ~/myScripts/add_any_2files_together.pl <(cut -f 4 hg19.refGene.tss.rand273.gene.bed) rand273.200kb.neighbors.txt 0 3 | awk '$14=="/"' | cut -f 1-13 | uniq | grep -v MIR > rand273.200kb.nonpseudoncogene.neighbors.txt
###### merge alternative promoter to genes ######
paste <(cut -f 4 rand273.200kb.nonpseudoncogene.neighbors.txt | cut -f 1 -d "|") <(cut -f 10 rand273.200kb.nonpseudoncogene.neighbors.txt | cut -f 1 -d "|") | sort | uniq | cut -f 2 | sort | uniq -c | awk '{print $2"\t"$1}' | sort -k2,2nr > rand273.200kb.nonpseudoncogene.neighbors.stats
perl ~/myScripts/add_any_2files_together.pl 13cancers.pMUT rand273.200kb.nonpseudoncogene.neighbors.txt 4 3 > rand273.200kb.nonpseudoncogene.neighbors.pMUT
awk -vOFS="\t" '$17!="/"{print $14,$15,$16,$17,$4,$10}' rand273.200kb.nonpseudoncogene.neighbors.pMUT | sort | uniq > rand273.200kb.nonpseudoncogene.neighbors.pMUT.sim
cut -f 1-4,6 rand273.200kb.nonpseudoncogene.neighbors.pMUT.sim | sed 's/|/\t/g' | cut -f 1-4,6 | sort | uniq | cut -f 4-5 | sort | uniq -c | awk '{print $2"\t"$3"\t"$1}' | sort -k3,3nr > rand273.200kb.nonpseudoncogene.neighbors.pMUT.stats

####################################################################################
############################## get the motif information ##############################
# for pre in `cut -d " " -f2 mut2promoter/selected.cancers/all.cancer.wgs1M.cancertypes`
for pre in MELA-AU
do
	mut2tss=mut2promoter/selected.cancers/$pre.dis2tss.tsv
	# for motif in rand
	for motif in ESR1 FOXA1 GATA3 MA0139.1
	do
		cut -f 2-7 $mut2tss | bedtools closest -D b -a - -b hg19.fimo.$motif.mid.srt.bed > $pre.$motif.dis.txt &
	done
done
for pre in ESR1 FOXA1 GATA3 MA0139.1 rand
do
	cat *$pre.dis.txt | awk -F"\t" -vOFS="\t" '$7!="." && $13>-20 && $13<20{print $1,$2,$3,$4,$6}' | uniq > $pre.motifmut2tss.dis &
done
   # 515674 ESR1.motifmut2tss.dis
   # 788049 FOXA1.motifmut2tss.dis
   # 566188 GATA3.motifmut2tss.dis
   # 624196 MA0139.1.motifmut2tss.dis
   # 637906 rand.motifmut2tss.dis


for mutations in simple_somatic_mutation.open.*.tsv.gz.WGS.srt.bed
do
	outpre=`echo $mutations | sed 's/simple_somatic_mutation.open.//;s/.tsv.gz.WGS.srt.bed//;'`
	bedtools closest -D b -a <(uniq $mutations) -b $tssBED | awk -F "\t" -v var=$outpre -v OFS="\t" '$(NF-6)!="." {print var,$1,$2,$3,$4,$10,$13}' > ${outpre}.dis2tss.tsv &
done


# cat *.tsv > all.cancer.wgs.dis2tss.tsv
wc -l *.dis2tss.tsv | awk '$1 > 5000 {print $2}' | grep -v total | xargs -n 1 -I old cat old >> all.cancer.wgs.dis2tss.tsv
wc -l *.dis2tss.tsv | awk '$1 > 100000 {print $2}' | grep -v total | grep -v "all.cancer.wgs" | xargs -n 1 -I old cat old >> all.cancer.wgs100k.dis2tss.tsv
wc -l *.dis2tss.tsv | awk '$1 > 1000000 {print $2}' | grep -v total | grep -v "all.cancer.wgs" | xargs -n 1 -I old cat old >> all.cancer.wgs1M.dis2tss.tsv

cut -f 1 all.cancer.wgs1M.dis2tss.tsv | uniq -c > all.cancer.wgs1M.cancertypes

# 72 cancer type
# 63 with WGS
# 13 with WGS & mutation > 1M

#          5 simple_somatic_mutation.open.ALL-US.tsv.gz.WGS.srt.bed
#       6445 simple_somatic_mutation.open.BLCA-US.tsv.gz.WGS.srt.bed
#      36670 simple_somatic_mutation.open.BOCA-FR.tsv.gz.WGS.srt.bed
#     210646 simple_somatic_mutation.open.BOCA-UK.tsv.gz.WGS.srt.bed
#    3973947 simple_somatic_mutation.open.BRCA-EU.tsv.gz.WGS.srt.bed
#     622660 simple_somatic_mutation.open.BRCA-FR.tsv.gz.WGS.srt.bed
#     378575 simple_somatic_mutation.open.BRCA-UK.tsv.gz.WGS.srt.bed
#       9651 simple_somatic_mutation.open.BRCA-US.tsv.gz.WGS.srt.bed
#     134159 simple_somatic_mutation.open.BTCA-SG.tsv.gz.WGS.srt.bed
#       1588 simple_somatic_mutation.open.CESC-US.tsv.gz.WGS.srt.bed
#     415646 simple_somatic_mutation.open.CLLE-ES.tsv.gz.WGS.srt.bed
#      51144 simple_somatic_mutation.open.CMDI-UK.tsv.gz.WGS.srt.bed
#      56122 simple_somatic_mutation.open.COAD-US.tsv.gz.WGS.srt.bed
#     521271 simple_somatic_mutation.open.COCA-CN.tsv.gz.WGS.srt.bed
#       1427 simple_somatic_mutation.open.DLBC-US.tsv.gz.WGS.srt.bed
#     499691 simple_somatic_mutation.open.EOPC-DE.tsv.gz.WGS.srt.bed
#    8104957 simple_somatic_mutation.open.ESAD-UK.tsv.gz.WGS.srt.bed
#       1089 simple_somatic_mutation.open.ESCA-CN.tsv.gz.WGS.srt.bed
#     564986 simple_somatic_mutation.open.GACA-CN.tsv.gz.WGS.srt.bed
#       6684 simple_somatic_mutation.open.GBM-US.tsv.gz.WGS.srt.bed
#       8207 simple_somatic_mutation.open.HNSC-US.tsv.gz.WGS.srt.bed
#       1993 simple_somatic_mutation.open.KICH-US.tsv.gz.WGS.srt.bed
#       2622 simple_somatic_mutation.open.KIRC-US.tsv.gz.WGS.srt.bed
#       2704 simple_somatic_mutation.open.KIRP-US.tsv.gz.WGS.srt.bed
#      13479 simple_somatic_mutation.open.LAML-KR.tsv.gz.WGS.srt.bed
#        568 simple_somatic_mutation.open.LGG-US.tsv.gz.WGS.srt.bed
#      41689 simple_somatic_mutation.open.LIAD-FR.tsv.gz.WGS.srt.bed
#     422195 simple_somatic_mutation.open.LICA-CN.tsv.gz.WGS.srt.bed
#    1246616 simple_somatic_mutation.open.LICA-FR.tsv.gz.WGS.srt.bed
#       5978 simple_somatic_mutation.open.LIHC-US.tsv.gz.WGS.srt.bed
#     443958 simple_somatic_mutation.open.LINC-JP.tsv.gz.WGS.srt.bed
#    3883294 simple_somatic_mutation.open.LIRI-JP.tsv.gz.WGS.srt.bed
#    3668778 simple_somatic_mutation.open.LMS-FR.tsv.gz.WGS.srt.bed
#      12157 simple_somatic_mutation.open.LUAD-US.tsv.gz.WGS.srt.bed
#        423 simple_somatic_mutation.open.LUSC-CN.tsv.gz.WGS.srt.bed
#    1717259 simple_somatic_mutation.open.LUSC-KR.tsv.gz.WGS.srt.bed
#      17614 simple_somatic_mutation.open.LUSC-US.tsv.gz.WGS.srt.bed
#    3274936 simple_somatic_mutation.open.MALY-DE.tsv.gz.WGS.srt.bed
#   23536812 simple_somatic_mutation.open.MELA-AU.tsv.gz.WGS.srt.bed
#        145 simple_somatic_mutation.open.NBL-US.tsv.gz.WGS.srt.bed
#         23 simple_somatic_mutation.open.NKTL-SG.tsv.gz.WGS.srt.bed
#     246891 simple_somatic_mutation.open.ORCA-IN.tsv.gz.WGS.srt.bed
#    1070420 simple_somatic_mutation.open.OV-AU.tsv.gz.WGS.srt.bed
#       3854 simple_somatic_mutation.open.OV-US.tsv.gz.WGS.srt.bed
#    1366211 simple_somatic_mutation.open.PACA-AU.tsv.gz.WGS.srt.bed
#    1969896 simple_somatic_mutation.open.PACA-CA.tsv.gz.WGS.srt.bed
#     172034 simple_somatic_mutation.open.PAEN-AU.tsv.gz.WGS.srt.bed
#     154291 simple_somatic_mutation.open.PAEN-IT.tsv.gz.WGS.srt.bed
#    1143929 simple_somatic_mutation.open.PBCA-DE.tsv.gz.WGS.srt.bed
#     570176 simple_somatic_mutation.open.PRAD-CA.tsv.gz.WGS.srt.bed
#     124179 simple_somatic_mutation.open.PRAD-FR.tsv.gz.WGS.srt.bed
#     862099 simple_somatic_mutation.open.PRAD-UK.tsv.gz.WGS.srt.bed
#        451 simple_somatic_mutation.open.PRAD-US.tsv.gz.WGS.srt.bed
#      24365 simple_somatic_mutation.open.READ-US.tsv.gz.WGS.srt.bed
#     761813 simple_somatic_mutation.open.RECA-EU.tsv.gz.WGS.srt.bed
#       2057 simple_somatic_mutation.open.SARC-US.tsv.gz.WGS.srt.bed
#    3445643 simple_somatic_mutation.open.SKCA-BR.tsv.gz.WGS.srt.bed
#      36384 simple_somatic_mutation.open.SKCM-US.tsv.gz.WGS.srt.bed
#      20272 simple_somatic_mutation.open.STAD-US.tsv.gz.WGS.srt.bed
#      22445 simple_somatic_mutation.open.THCA-SA.tsv.gz.WGS.srt.bed
#        874 simple_somatic_mutation.open.THCA-US.tsv.gz.WGS.srt.bed
#      36136 simple_somatic_mutation.open.UCEC-US.tsv.gz.WGS.srt.bed
#      64252 simple_somatic_mutation.open.UTCA-FR.tsv.gz.WGS.srt.bed

# ALL-US 2
# BLCA-US 23
# BOCA-FR 98
# BOCA-UK 64
# BRCA-EU 569
# BRCA-FR 72
# BRCA-UK 45
# BRCA-US 91
# BTCA-SG 12
# CESC-US 20
# CLLE-ES 151
# CMDI-UK 30
# COAD-US 44
# COCA-CN 30
# DLBC-US 7
# EOPC-DE 202
# ESAD-UK 253
# ESCA-CN 17
# GACA-CN 37
# GBM-US 41
# HNSC-US 44
# KICH-US 45
# KIRC-US 37
# KIRP-US 33
# LAML-KR 8
# LGG-US 18
# LIAD-FR 5
# LICA-CN 25
# LICA-FR 39
# LIHC-US 54
# LINC-JP 31
# LIRI-JP 258
# LMS-FR 67
# LUAD-US 38
# LUSC-CN 10
# LUSC-KR 30
# LUSC-US 48
# MALY-DE 241
# MELA-AU 183
# NBL-US 42
# NKTL-SG 23
# ORCA-IN 25
# OV-AU 93
# OV-US 42
# PACA-AU 252
# PACA-CA 209
# PAEN-AU 50
# PAEN-IT 37
# PBCA-DE 457
# PRAD-CA 124
# PRAD-FR 25
# PRAD-UK 215
# PRAD-US 19
# READ-US 16
# RECA-EU 95
# SARC-US 34
# SKCA-BR 80
# SKCM-US 37
# STAD-US 38
# THCA-SA 129
# THCA-US 48
# UCEC-US 51
# UTCA-FR 7

####################################################################################
#################### get the motif for reccurent promoter mutation ####################
cat *.pMUT | cut -f 1-3,5-6 | sort | uniq -c | awk '$1>1' > all.cancer.pMUT.counts
sort -k1,1nr all.cancer.pMUT.counts  | grep -v MIR | more

awk '$1 >= 5 && $6 < 100 && $6 > -200' all.cancer.pMUT.counts | grep -v "MIR" | wc -l
awk -v OFS="\t" '$1 >= 5 && $6 < 100 && $6 > -200{print $2,$3,$4,$5,$6,$1}' all.cancer.pMUT.counts | grep -v "MIR" | bedtools sort -i > all.cancer.recurrentpmut.bed

# for pre in `s ../*-*.pMUT | cut -f 2 -d "/" | sed 's/.pMUT//'`
for pre in ALL-US BLCA-US BOCA-FR BOCA-UK BRCA-EU BRCA-FR BRCA-UK BRCA-US BTCA-SG CESC-US CLLE-ES CMDI-UK COAD-US COCA-CN DLBC-US EOPC-DE ESAD-UK ESCA-CN GACA-CN GBM-US HNSC-US KICH-US KIRC-US KIRP-US LAML-KR LGG-US LIAD-FR LICA-CN LICA-FR LIHC-US LINC-JP LIRI-JP LMS-FR LUAD-US LUSC-CN LUSC-KR LUSC-US MALY-DE MELA-AU NBL-US NKTL-SG ORCA-IN OV-AU OV-US PACA-AU PACA-CA PAEN-AU PAEN-IT PBCA-DE PRAD-CA PRAD-FR PRAD-UK PRAD-US READ-US RECA-EU SARC-US SKCA-BR SKCM-US STAD-US THCA-SA THCA-US UCEC-US UTCA-FR
do
	rawmutation=../../../simple_somatic_mutation.open.$pre.tsv.gz
	gunzip -c $rawmutation | awk -F"\t" -vOFS="\t" '$34=="WGS"' | cut -f 2-3,9-11,14,16-17 | uniq > $pre.WGS.rawmut.txt &
done

for pre in `ls *.WGS.rawmut.txt | sed 's/.WGS.rawmut.txt//'`
do
	perl ~/myScripts/add_any_2files_together.pl <(awk -F"\t" -vOFS="\t" '{print $3":"$5,$0}' $pre.WGS.rawmut.txt) <(awk -F"\t" -vOFS="\t" '{print $1":"$3,$0}' all.cancer.recurrentpmut.bed) 0 0 | awk '$8!="/"' > $pre.recurrentpmut.txt &
done
# more SKCA-BR.recurrentpmut.txt | grep -w "39636928"
# cat *.recurrentpmut.txt | cut -f 2-7,9-10,14- > all.cancer.recurrentpmut.info
# cut -f 1-5,10-11 all.cancer.recurrentpmut.info | sort | uniq > all.cancer.recurrentpmut.uniq.info
cat *.recurrentpmut.txt | awk -F"\t" -vOFS="\t" '{print $11,$12-1,$13,$5,$6,$7,$9,$10,$14,$15,$16}' | sort | uniq  > all.cancer.recurrentpmut.uniq.info

echo "Hugo_Symbol Chromosome Start_Position End_Position Variant_Classification Variant_Type Tumor_Sample_Barcode Gene_Type Reference_Allele Tumor_Seq_Allele1 Tumor_Seq_Allele2 Cancer_Type" | tr " " "\t" > all.cancer.recurrentpmut.u200d100.pmut.maf
paste <(cut -f 4 all.cancer.recurrentpmut.uniq.info | cut -f1 -d "|") <(awk -F"\t" -vOFS="\t" '{print $1,$2+1,$3,$9,$9,$7,"CorePromoter",$10,$10,$11,$8}' all.cancer.recurrentpmut.uniq.info) >> all.cancer.recurrentpmut.u200d100.pmut.maf


##############################
cut -f 1-4,10-11 all.cancer.recurrentpmut.u200d100.pmut.maf | uniq > all.cancer.recurrentpmut.u200d100.pmut.unique.info

pre=all.cancer.recurrentpmut.u200d100.pmut.unique.info
refseq=/home1/04935/shaojf/scratch/bwa-index/hg19.fa
JASPAR=/home1/04935/shaojf/stampede2/myTools/MEME/motif_databases/JASPAR/JASPAR_CORE_2016_vertebrates.meme
HOCOMOCO=/home1/04935/shaojf/stampede2/myTools/MEME/motif_databases/HUMAN/HOCOMOCOv10_HUMAN_mono_meme_format.meme

# bedtools getfasta -fi $refseq -bed <(echo "10 111797697 111797697" | awk -vOFS="\t" '{print "chr"$1,$2-1,$3}') 
# bedtools getfasta -fi $refseq -bed <(echo "10 111797697 111797697" | awk -vOFS="\t" '{print "chr"$1,$2-21,$3+20}') 
# >chr10:111797696-111797697
# c
# >chr10:111797676-111797717
# tagctgggtgtgtggtggcacatgcctgtaatcccagctac
# >chr10:111797676-111797696
# tagctgggtgtgtggtggca
# >chr10:111797697-111797717
#                      atgcctgtaatcccagctac
# flanking 20 bp
awk -vOFS="\t" '{print "chr"$2,$3-21,$4+20}' $pre | tail -n +2 | bedtools getfasta -fi $refseq -bed - | sed -n '2~2p' | paste <(awk -F"\t" '{print $1"|"$2":"$3"-"$4"|"$5"/"$6}' $pre | tail -n +2) - | awk '{print ">"$1"\n"$2}' > $pre.ref.fa

awk -vOFS="\t" '{print "chr"$2,$3-21,$3-1}' $pre | tail -n +2 | bedtools getfasta -fi $refseq -bed - > $pre.left.fa
awk -vOFS="\t" '{if($5=="-"){print "chr"$2,$4-1,$4+20}else{print "chr"$2,$4,$4+20}}' $pre | tail -n +2 | bedtools getfasta -fi $refseq -bed - > $pre.right.fa
awk -F"\t" '{print ">"$1"|"$2":"$3"-"$4"\n"$6}' $pre | tail -n +3 | paste $pre.left.fa - $pre.right.fa | sed 's/\t//g' | sed '2~2s/-//g' > $pre.mut.fa
sed -n '2~2p' $pre.mut.fa > $pre.seq
awk -F"\t" '{print ">"$1"|"$2":"$3"-"$4"|"$5"/"$6}' $pre | tail -n +2 | sed "R $pre.seq" > $pre.a
mv $pre.a $pre.mut.fa
rm $pre.seq $pre.left.fa $pre.right.fa
# ctcttttttttttttttttttAAACACATTTTTTTCCTGGC
# ctctttttttttttttttttTTTTTtAAACACATTTTTTTCCTGGC

# TTTATTTAttgtttttttttttgtttgtttttgtttttgag
# TTTATTTAttgtttttttttGTTttgtttgtttttgtttttgag

# fimo
~/.local/bin/fimo -oc JASPAR.$pre.ref $JASPAR $pre.ref.fa
~/.local/bin/fimo -oc JASPAR.$pre.mut $JASPAR $pre.mut.fa
~/.local/bin/fimo -oc HOCOMOCO.$pre.ref $HOCOMOCO $pre.ref.fa
~/.local/bin/fimo -oc HOCOMOCO.$pre.mut $HOCOMOCO $pre.mut.fa

sed -i 's/\t/|/' JASPAR.$pre.ref/fimo.txt
sed -i 's/\t/|/' JASPAR.$pre.mut/fimo.txt
sed -i 's/\t//' HOCOMOCO.$pre.ref/fimo.txt
sed -i 's/\t//' HOCOMOCO.$pre.mut/fimo.txt
cat JASPAR.$pre.ref/fimo.txt HOCOMOCO.$pre.ref/fimo.txt > $pre.ref.fimo.txt
cat JASPAR.$pre.mut/fimo.txt HOCOMOCO.$pre.mut/fimo.txt > $pre.mut.fimo.txt
## may need big memory
# cut -f 1-2 $pre.ref.fimo.txt | grep -vf - $pre.mut.fimo.txt > $pre.motif.gain
# cut -f 1-2 $pre.mut.fimo.txt | grep -vf - $pre.ref.fimo.txt > $pre.motif.lose

# rm $pre.mut.fa $pre.ref.fa 
# rm $pre.left.fa $pre.right.fa
# rm -rf JASPAR.$pre.ref/ JASPAR.$pre.mut/ HOCOMOCO.$pre.ref/ HOCOMOCO.$pre.mut/

mycodes=/home1/04935/shaojf/myTools/UTHealthMut2Loop/GTEx_codes
perl $mycodes/motif_compare.pl $pre

ln -s JASPAR.$pre.ref/fimo.txt all.cancer.recurrentpmut.u200d100.jaspar.ref.fimo.txt
ln -s JASPAR.$pre.mut/fimo.txt all.cancer.recurrentpmut.u200d100.jaspar.mut.fimo.txt
perl $mycodes/motif_compare.pl all.cancer.recurrentpmut.u200d100.jaspar
######## format a matrix ########
perl GTEx.snp2motif.makemat.pl all.cancer.recurrentpmut.u200d100.jaspar.motif.cmp > all.cancer.recurrentpmut.u200d100.jaspar.motif.mat
######## format a matrix ########
# library(pheatmap)
# x <- read.table("all.cancer.recurrentpmut.u200d100.jaspar.motif.mat", header = T, row.names = 1)
# zerocount <- apply(x, 2, function(x){length(x[x==0])})
# # zerocount[zerocount < summary(zerocount)[2]]
# a <- 1 - zerocount/nrow(x)
# data.x <- x[, a > 0.02]
# grep("CTCF", colnames(data.x))
# zerocount <- apply(data.x, 1, function(x){length(x[x==0])})
# data.y <- data.x[zerocount < ncol(data.x), ]

# colors <- colorRampPalette(c("blue", "white", "red"))(3)
# png(filename = paste("all.cancer.recurrentpmut.u200d100.jaspar", "pheatmap.cluster.png", sep = "."), width = 1000, height = 1800)
# myplot <- pheatmap(data.y, scale = "none", 
# 	show_rownames = F, show_colnames = T, color = colors,
# 	cluster_cols = T, cluster_rows = T)
# dev.off()

# pdf(file = paste("all.cancer.recurrentpmut.u200d100.jaspar", "pheatmap.cluster.pdf", sep = "."), width = 10, height = 18)
# myplot <- pheatmap(data.y, scale = "none", fontsize_row = 1,
# 	show_rownames = T, show_colnames = T, color = colors,
# 	cluster_cols = T, cluster_rows = T)
# dev.off()
#######################
findMotifsGenome.pl <(awk -vOFS="\t" '{print "chr"$2,$3-21,$4+20,$1"|"$2":"$3"-"$4"|"$5"/"$6,"1000","+"}' all.cancer.recurrentpmut.u200d100.pmut.unique.info | tail -n +2) hg19 all.cancer.recurrentpmut.u200d100.pmut.unique.homer.motifs -size given -mknown /home1/04935/shaojf/myTools/HOMER/data/knownTFs/vertebrates/known.motifs > all.cancer.recurrentpmut.u200d100.pmut.unique.findMotifsGenome.txt &
annotatePeaks.pl <(awk -vOFS="\t" '{print "chr"$2,$3-21,$4+20,$1"|"$2":"$3"-"$4"|"$5"/"$6,"1000","+"}' all.cancer.recurrentpmut.u200d100.pmut.unique.info | tail -n +2) hg19 -m /home1/04935/shaojf/myTools/HOMER/data/knownTFs/vertebrates/known.motifs -size given 1> all.cancer.recurrentpmut.u200d100.pmut.unique.motifs.txt 2> all.cancer.recurrentpmut.u200d100.pmut.unique.motifs.log &
findMotifsGenome.pl <(awk -vOFS="\t" '{print "chr"$2,$3-21,$3+20,$1"|"$2":"$3"-"$4,"1000","+"}' all.cancer.recurrentpmut.u200d100.pmut.unique.info | uniq | tail -n +2) hg19 all.cancer.recurrentpmut.u200d100.pmut.unique.homer.motifs.u -size given -mknown /home1/04935/shaojf/myTools/HOMER/data/knownTFs/vertebrates/known.motifs > all.cancer.recurrentpmut.u200d100.pmut.unique.u.findMotifsGenome.txt &
####################################################################################


