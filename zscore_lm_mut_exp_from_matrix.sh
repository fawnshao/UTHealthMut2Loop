#!/bin/sh
# cd /home1/04935/shaojf/stampede2/mutations/mut_exp_lm/
module load Rstats
module load bedtools
np=33
pre=$1
promoterLEN=1000
winSIZE=500000
matrixdir=/home1/04935/shaojf/stampede2/mutations/ICGC/lm/mut_exp_matrix
bindir=/home1/04935/shaojf/myTools/UTHealthMut2Loop
tssBED=/home1/04935/shaojf/stampede2/refs/hg19.refGene.tss.uniq.srt.bed

date
echo $pre
awk -F "\t" -v var=$winSIZE -v OFS="\t" \
'{a=$2-var;b=$3+var;if(a<0){a=0}print $1,a,b,$4,$5,$6}' $tssBED > $tssBED.win.4.$pre

cut -f 1 $matrixdir/$pre.bothWGS.mut.tsv | \
sed 's/:/\t/g' | sed -n '2,$p' | bedtools intersect -wao -a - -b $tssBED.win.4.$pre | \
awk '$NF > 0 {print $1":"$2":"$3"\t"$7}' | \
perl $bindir/relatedScripts/make_mut_exp_pairs.pl /dev/stdin > $pre.win.mut.tsv

head -1 $matrixdir/$pre.bothWGS.mut.tsv > $pre.p.mut.tsv
cut -f 1 $matrixdir/$pre.bothWGS.mut.tsv | \
sed 's/:/\t/g' | sed -n '2,$p' | bedtools sort -i - | \
bedtools closest -D b -a - -b $tssBED | \
awk -F "\t" -v var=$promoterLEN -v OFS="\t" \
'$(NF-6)!="." && $NF > -1 * var && $NF < var {print $1":"$2":"$3,$7"~"$10}' | \
perl $bindir/relatedScripts/add_any_2files_together.pl $matrixdir/$pre.bothWGS.mut.tsv /dev/stdin  0 0 | \
cut -f 1-2,4- | sed 's/\t/:/' >> $pre.p.mut.tsv









count=`wc -l $matrixdir/$pre.bothWGS.mut.tsv | awk '{print $1}'`
each=`echo $count/$np | bc`
sed -n '2,$p' $matrixdir/$pre.bothWGS.mut.tsv | split -l $each /dev/stdin $pre.tmp.
for f in $pre.tmp.*
do
	head -1 $matrixdir/$pre.bothWGS.mut.tsv | cat - $f > $f.title
	Rscript $bindir/relatedScripts/lm_mut_exp_from_matrix.R $matrixdir/$pre.bothWGS.exp.tsv $f.title $f &
done
wait
cat $pre.tmp.*.exp2mut.lm.tsv > $pre.exp2mut.lm.tsv
rm $pre.tmp.*
echo ++++++++finished+++++++
date

