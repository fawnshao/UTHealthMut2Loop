#!/bin/sh
# input coordinates, get the flanking sequence
# scan motif with the sequences
# find if a mutation will gain to lose motifs
outpre=$1
motifDIR=/home1/04935/shaojf/stampede2/mutations/mut_exp_lm/motifdir
promoterMOTIF=/home1/04935/shaojf/stampede2/mutations/ICGC/p_motif/hg19.tss1k.motif.simcount

#####
## find if the mutation is associated with some motifs
echo +++++++++ checking if the mutation is in any motif  ++++++++
# echo "mutation-motif" > ${outpre}.LoopBroken.motif
for m in `ls $motifDIR/`
do
	bedtools intersect -wao -a ${outpre}.LoopBroken.bed -b $motifDIR/$m | \
	awk '$NF > 0' >> ${outpre}.LoopBroken.motif
done

echo "mutated.sites	gene	promoter.motif.count	promoter.motifs	mutated.motifs" > $outpre.motif.tsv
while read line
do
	mutsite=`echo $line | awk '{print $1}' | sed 's/"//g'`
	gene=`echo $line | awk '{print $2}' | sed 's/"//g'`
	motifcount=`echo "" | awk -v a=$gene '{print "\t"a"|\n;"a"|"}' | \
	grep -f - $promoterMOTIF | awk -F"\t" '{print $8}' |  tr ', ' '\n' | \
	grep -v "^$" | sort | uniq | wc -l`
	motifs=`echo "" | awk -v a=$gene '{print "\t"a"|\n;"a"|"}' | \
	grep -f - $promoterMOTIF | awk -F"\t" '{print $8}' |  tr ', ' '\n' | \
	grep -v "^$" | sort | uniq | tr '\n' ',' | sed 's/,$//'`
	mutmotif=`awk -v a=$mutsite '$4==a{print $8}' ${outpre}.LoopBroken.motif | \
	sort | uniq | tr '\n' ',' | sed 's/,$//'`
	echo $mutsite"	"$gene"	"$motifcount"	"$motifs"	"$mutmotif >> $outpre.motif.tsv
done < $outpre.exp2mut.lm.tsv