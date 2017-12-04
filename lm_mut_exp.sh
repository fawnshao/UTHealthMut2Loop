#!/bin/sh
bindir=/home1/04935/shaojf/myTools/UTHealthMut2Loop
module load Rstats
# TSS should be 6-column sorted bed files
# the TAD file is 4-column bed files
# default parameters
#input
# tssBED=hg19.refGene.tss.uniq.srt.bed
# tadBED=hg19.GSE63525_GM12878.srt.bed
# promoterLEN=1000
pre=$1
matrixdir=/home1/04935/shaojf/stampede2/mutations/ICGC/lm/mut_exp_matrix
np=69

# head -1 $matrixdir/$pre.bothWGS.mut.tsv > $pre.p.mut.tsv
# cut -f 1 $matrixdir/$pre.bothWGS.mut.tsv | \
# sed 's/:/\t/g' | sed -n '2,$p' | bedtools sort -i - | \
# bedtools closest -D b -a - -b $tssBED | \
# awk -F "\t" -v var=$promoterLEN -v OFS="\t" \
# '$(NF-6)!="." && $NF > -1 * var && $NF < var {print $1":"$2":"$3,$7"~"$10}' | \
# perl $bindir/relatedScripts/add_any_2files_together.pl $matrixdir/$pre.bothWGS.mut.tsv /dev/stdin 0 0 | cut -f 2,4- >> $pre.p.mut.tsv

count=`wc -l ${pre}.IamGroot.Rinput | awk '{print $1}'`
echo +++++++++ Running Rscript to output loop translocate candidates  ++++++++
# use Z score to find the expression alteration direction in the TAD
if [ $count -lt 100 ]
	then
	echo Running in one piece
	Rscript $bindir/ICGC_mut2exp_lm.R $matrixdir/$pre.bothWGS.exp.tsv $matrixdir/$pre.bothWGS.mut.tsv $pre.IamGroot.Rinput $pre
else
	echo Running in many pieces
	date
	each=`echo $count/$np | bc`
	sed -n '2,$p' ${pre}.IamGroot.Rinput | split -l $each /dev/stdin ${pre}.IamGroot.Rinput.
	for f in ${pre}.IamGroot.Rinput.*
	do
		echo "TAD sample mutsample genes mutgene" | cat - $f > $f.title
		rm $f
		Rscript $bindir/ICGC_mut2exp_lm.R $matrixdir/$pre.bothWGS.exp.tsv $matrixdir/$pre.bothWGS.mut.tsv $f.title $f &
	done
	wait
	echo Finishing all pieces
	date
	cat ${pre}.IamGroot.Rinput.*.TAD.labeled.tsv > ${pre}.TAD.labeled.tsv
	for f in ${pre}.IamGroot.Rinput.*.TAD_*.tsv
	do
		new=`echo $f | awk -F"." -vOFS="." '{print $1,$5,$6,$7}'`
		mv $f $new
	done
	rm ${pre}.IamGroot.Rinput.*
fi