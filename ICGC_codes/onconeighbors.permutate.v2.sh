#!/bin/bash
# /home1/04935/shaojf/stampede2/mutations/wenbo.paper.check/ICGC.WGS.mut2oncogene
# only gene name
oncogenelist=Cosmic.CancerGeneCensus.oncogene.gene
reftss=hg19.refGene.tss.uniq.srt.bed
myperl=/home1/04935/shaojf/myTools/BioinformaticsDaily/textProcess/add_any_2files_together.pl

############ get the neighboring promoter mutations for nononcogene #############
# once done, no need to run any more
# grep -wf $oncogenelist $reftss | grep -v "\-AS" | \
# 	grep -v "\-IT1" | grep -v "\-OT" | \
# 	grep -v "HOXA10-HOXA9" > hg19.refGene.tss.Cosmic.CancerGeneCensus.oncogene.gene.bed
# # cut -f 4 hg19.refGene.tss.Cosmic.CancerGeneCensus.oncogene.gene.bed | sort | uniq | cut -f 1 -d"|" | sort | uniq > hg19.refGene.tss.Cosmic.CancerGeneCensus.oncogene.gene.withTSS
# # grep -vwf hg19.refGene.tss.Cosmic.CancerGeneCensus.oncogene.gene.withTSS <(awk '{print $1}' Cosmic.CancerGeneCensus.oncogene.gene)
# # MLLT4
# # WHSC1
# # WHSC1L1
# # cut -f 4 hg19.refGene.tss.Cosmic.CancerGeneCensus.oncogene.gene.bed | tr ";" "\n" | cut -f 1 -d "|" | sort | uniq | wc -l
# cut -f 4 hg19.refGene.tss.Cosmic.CancerGeneCensus.oncogene.gene.bed | \
# 	grep -wf /dev/stdin 13cancers.pMUT > 13cancers.oncogene.pMUT
# sed 's/|/\t/' 13cancers.oncogene.pMUT | cut -f 4,6 | \
# 	sort | uniq -c | awk '{print $2"\t"$3"\t"$1}' | sort -k3,3nr > 13cancers.oncogene.pMUT.stats

# ###### 200kb neighboring TSS (exclude other oncogenes and oncogene itself) to oncogene TSS
# bedtools closest -D b -io -a $reftss -b hg19.refGene.tss.Cosmic.CancerGeneCensus.oncogene.gene.bed | \
# 	awk '$7!="." && (($13 > -200000 && $13 < -1000) || ($13 > 1000 && $13 < 200000))' \
# 	> oncogene.200kb.neighbors.txt
# # grep -wf <(cut -f 4 hg19.refGene.tss.Cosmic.CancerGeneCensus.oncogene.gene.bed) <(cut -f 4 oncogene.200kb.neighbors.txt)
# perl $myperl <(cut -f 4 hg19.refGene.tss.Cosmic.CancerGeneCensus.oncogene.gene.bed) \
# 	oncogene.200kb.neighbors.txt 0 3 | awk '$14!="/"' | \
# 	cut -f 1-13 | uniq > oncogene.200kb.oncogene.neighbors.txt
# ##### remove miR
# perl $myperl <(cut -f 4 hg19.refGene.tss.Cosmic.CancerGeneCensus.oncogene.gene.bed) \
# 	oncogene.200kb.neighbors.txt 0 3 | awk '$14=="/"' | \
# 	cut -f 1-13 | uniq | grep -v MIR > oncogene.200kb.nononcogene.neighbors.txt
# # cut -f 10 oncogene.200kb.nononcogene.neighbors.txt | uniq -c | awk '{print $2"\t"$1}' | sort -k2,2nr > oncogene.200kb.nononcogene.neighbors.stats
# ###### merge alternative promoter to genes ######
# paste <(cut -f 4 oncogene.200kb.nononcogene.neighbors.txt | cut -f 1 -d "|") \
# 	<(cut -f 10 oncogene.200kb.nononcogene.neighbors.txt | cut -f 1 -d "|") | \
# 	sort | uniq | cut -f 2 | sort | uniq -c | awk '{print $2"\t"$1}' | \
# 	sort -k2,2nr > oncogene.200kb.nononcogene.neighbors.stats

# perl $myperl 13cancers.pMUT oncogene.200kb.nononcogene.neighbors.txt 4 3 \
# 	> oncogene.200kb.nononcogene.neighbors.pMUT
# # cut -f 4,10,14-17 
# awk -vOFS="\t" '$17!="/"{print $14,$15,$16,$17,$4,$10}' oncogene.200kb.nononcogene.neighbors.pMUT | \
# 	sort | uniq > oncogene.200kb.nononcogene.neighbors.pMUT.sim
# paste <(cut -f 1-4 oncogene.200kb.nononcogene.neighbors.pMUT.sim) \
# 	<(cut -f 5 oncogene.200kb.nononcogene.neighbors.pMUT.sim | cut -f1 -d"|") \
# 	<(cut -f 6 oncogene.200kb.nononcogene.neighbors.pMUT.sim | cut -f1 -d"|") | \
# 	sort | uniq | sed 's/|/\t/' | cut -f 4,6-7 | sort | uniq -c | \
# 	awk '{print $2"\t"$3"\t"$4"\t"$1}' | \
# 	sort -k4,4nr > oncogene.200kb.nononcogene.neighbors.pMUT.each.stats
# # sed 's/|/\t/' oncogene.200kb.nononcogene.neighbors.pMUT.sim | cut -f 4,7 | sort | uniq -c | awk '{print $1"\t"$2"\t"$3}' | sort -k1,1nr > oncogene.200kb.nononcogene.neighbors.pMUT.stats
# cut -f 1-4,6 oncogene.200kb.nononcogene.neighbors.pMUT.sim | \
# 	sed 's/|/\t/g' | cut -f 1-4,6 | sort | uniq | cut -f 4-5 | \
# 	sort | uniq -c | awk '{print $2"\t"$3"\t"$1}' | \
# 	sort -k3,3nr > oncogene.200kb.nononcogene.neighbors.pMUT.stats
# awk -vOFS="\t" '{print $1, "TOTAL", $2}' 13cancers.pMUT.stats >> oncogene.200kb.nononcogene.neighbors.pMUT.stats
# cut -f 4 hg19.refGene.tss.uniq.srt.bed | cut -f 1 -d "|" | \
# 	sort | uniq | wc -l | awk '{print "TOTAL\t"$1}' >> oncogene.200kb.nononcogene.neighbors.stats

# awk -vOFS="\t" '{print $1, "TOTAL", "NONE", $2/27502}' \
# 	13cancers.pMUT.stats >> oncogene.200kb.nononcogene.neighbors.pMUT.each.stats
# cut -f 2-3 oncogene.200kb.nononcogene.neighbors.pMUT.each.stats | \
# 	sort | uniq | sort -k2 | grep -v TOTAL > oncogene.200kb.nononcogene.neighbors.gene2gene
# perl $myperl oncogene.200kb.nononcogene.neighbors.stats \
# 	oncogene.200kb.nononcogene.neighbors.pMUT.stats 0 1 | \
# 	cut -f 1-3,5 | grep -v "TOTAL" \
# 	> oncogene.200kb.nononcogene.neighbors.pMUT.forPermutation.txt

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
echo "CancerType Oncogene NeighborspMUTCount NeighborspCount Permutaion" | tr " " "\t" > onconeighbors.forPermutation.txt
cat oncogene.200kb.nononcogene.neighbors.pMUT.forPermutation.txt | \
	awk '{print $0"\tOncogene"}' >> onconeighbors.forPermutation.txt
i=0
while [ $i -lt 1000 ]
do
	shuf $reftss | grep "NM" | \
		cut -f 4 | cut -f 1 -d "|" | \
		grep -wvf $oncogenelist | head -n 273 | sort | uniq > myrand273.gene
	grep -wf myrand273.gene $reftss | \
		grep -v "\-AS" > hg19.refGene.tss.rand273.gene.bed
	bedtools closest -D b -io -a $reftss -b hg19.refGene.tss.rand273.gene.bed | \
		awk '$7!="." && (($13 > -200000 && $13 < -1000) || ($13 > 1000 && $13 < 200000))' \
		> rand273.200kb.neighbors.txt
	perl $myperl <(cut -f 4 hg19.refGene.tss.rand273.gene.bed) rand273.200kb.neighbors.txt 0 3 | \
		awk '$14!="/"' | cut -f 1-13 | uniq > rand273.200kb.pseudoncogene.neighbors.txt
	##### remove miR
	perl $myperl <(cut -f 4 hg19.refGene.tss.rand273.gene.bed) rand273.200kb.neighbors.txt 0 3 | \
		awk '$14=="/"' | cut -f 1-13 | uniq | grep -v MIR > rand273.200kb.nonpseudoncogene.neighbors.txt
	###### merge alternative promoter to genes ######
	paste <(cut -f 4 rand273.200kb.nonpseudoncogene.neighbors.txt | cut -f 1 -d "|") \
		<(cut -f 10 rand273.200kb.nonpseudoncogene.neighbors.txt | cut -f 1 -d "|") | \
		sort | uniq | cut -f 2 | sort | uniq -c | awk '{print $2"\t"$1}' | \
		sort -k2,2nr > rand273.200kb.nonpseudoncogene.neighbors.stats
	perl $myperl 13cancers.pMUT rand273.200kb.nonpseudoncogene.neighbors.txt \
		4 3 > rand273.200kb.nonpseudoncogene.neighbors.pMUT
	awk -vOFS="\t" '$17!="/"{print $14,$15,$16,$17,$4,$10}' \
		rand273.200kb.nonpseudoncogene.neighbors.pMUT | sort | \
		uniq > rand273.200kb.nonpseudoncogene.neighbors.pMUT.sim
	cut -f 1-4,6 rand273.200kb.nonpseudoncogene.neighbors.pMUT.sim | \
		sed 's/|/\t/g' | cut -f 1-4,6 | sort | uniq | cut -f 4-5 | \
		sort | uniq -c | awk '{print $2"\t"$3"\t"$1}' | \
		sort -k3,3nr > rand273.200kb.nonpseudoncogene.neighbors.pMUT.stats
	perl $myperl rand273.200kb.nonpseudoncogene.neighbors.stats \
		rand273.200kb.nonpseudoncogene.neighbors.pMUT.stats 0 1 | \
		cut -f 1-3,5 | awk -v var=$i '{print $0"\tPermutation."var}' >> onconeighbors.forPermutation.txt
	i=$[$i+1]
done

