#!/bin/sh
# cd /home1/04935/shaojf/stampede2/mutations/mut_exp_lm/
module load Rstats
module load bedtools
np=7
pre=$1
promoterLEN=1000
winSIZE=500000
matrixdir=/home1/04935/shaojf/stampede2/mutations/ICGC/lm/mut_exp_matrix
bindir=/home1/04935/shaojf/myTools/UTHealthMut2Loop
cordir=/home1/04935/shaojf/stampede2/mutations/ICGC/correlation
tssBED=/home1/04935/shaojf/stampede2/refs/hg19.refGene.tss.uniq.srt.bed
motifDIR=/home1/04935/shaojf/stampede2/mutations/mut_exp_lm/motifdir
promoterMOTIF=/home1/04935/shaojf/stampede2/mutations/ICGC/p_motif/hg19.tss1k.motif.simcount

date
echo $pre
awk -F "\t" -v var=$winSIZE -v OFS="\t" \
'{a=$2-var;b=$3+var;if(a<0){a=0}print $1,a,b,$4,$5,$6}' $tssBED > $tssBED.win.4.$pre

cut -f 1 $matrixdir/$pre.bothWGS.mut.tsv | \
sed 's/:/\t/g' | sed -n '2,$p' | bedtools intersect -wao -a - -b $tssBED.win.4.$pre | \
awk '$NF > 0 {print $1":"$2":"$3"\t"$7}' | \
perl $bindir/relatedScripts/make_mut_exp_pairs.pl /dev/stdin > $pre.mw.mut.tsv

head -1 $matrixdir/$pre.bothWGS.mut.tsv > $pre.p.mut.tsv
cut -f 1 $matrixdir/$pre.bothWGS.mut.tsv | \
sed 's/:/\t/g' | sed -n '2,$p' | bedtools sort -i - | \
bedtools closest -D b -a - -b $tssBED | \
awk -F "\t" -v var=$promoterLEN -v OFS="\t" \
'$(NF-6)!="." && $NF > -1 * var && $NF < var {print $1":"$2":"$3,$7"~"$10}' | \
perl $bindir/relatedScripts/add_any_2files_together.pl $pre.mw.mut.tsv /dev/stdin  0 0 | \
awk -vOFS="\t" 'NF==4{print $1,$2,$4}' | \
perl $bindir/relatedScripts/add_any_2files_together.pl $matrixdir/$pre.bothWGS.mut.tsv /dev/stdin  0 0 | \
cut -f 1-3,5- | sed 's/\t/%/' | sed 's/\t/%/' >> $pre.p.mut.tsv

# Rscript $bindir/relatedScripts/zscore_lm_mut_exp_from_matrix.R \
# $matrixdir/$pre.bothWGS.exp.tsv $pre.p.mut.tsv $cordir/correlated.$pre.tsv $pre

count=`wc -l $pre.p.mut.tsv | awk '{print $1}'`
each=`echo $count/$np | bc`
sed -n '2,$p' $pre.p.mut.tsv | split -l $each /dev/stdin $pre.tmp.
for f in $pre.tmp.*
do
	head -1 $pre.p.mut.tsv | cat - $f > $f.title
	Rscript $bindir/relatedScripts/zscore_lm_mut_exp_from_matrix.R $matrixdir/$pre.bothWGS.exp.tsv $f.title $cordir/correlated.$pre.tsv $f &
done
wait
cat $pre.tmp.*.exp2mut.lm.tsv | sort -r | uniq > $pre.exp2mut.lm.tsv
rm $pre.tmp.*

# grep -v "mutation" loop-shift/res/$pre.exp2mut.lm.tsv | cut -f 1 | sed 's/"//g' | sort | uniq | awk -F":" -vOFS="\t" '{print $1,$2,$3,$0}' > $type.LoopBroken.bed
grep -v "mutation" $pre.exp2mut.lm.tsv | cut -f 1 | sed 's/"//g' | sort | uniq | awk -F":" -vOFS="\t" '{print $1,$2,$3,$0}' > $type.LoopBroken.bed
$bindir/ICGC_fimo_denovo_motif_in_LoopShiftMut.2.sh $type

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
	motifgained=`awk -v a=$mutsite '$2==a{print $1}' ${outpre}.LoopBroken.motif.gain | \
	sort | uniq | tr '\n' ',' | sed 's/,$//'`
	motiflost=`awk -v a=$mutsite '$2==a{print $1}' ${outpre}.LoopBroken.motif.lose | \
	sort | uniq | tr '\n' ',' | sed 's/,$//'`
	echo $mutsite"	"$gene"	"$motifcount"	"$motifs"	"$mutmotif"	"$motifgained"	"$motiflost >> $outpre.motif.tsv
done < $outpre.exp2mut.lm.tsv

echo ++++++++finished+++++++
date

