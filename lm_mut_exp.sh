#!/bin/sh
bindir=/home1/04935/shaojf/myTools/UTHealthMut2Loop
module load Rstats
# TSS should be 6-column sorted bed files
# the TAD file is 4-column bed files
# default parameters
#input
tssBED=hg19.refGene.tss.uniq.srt.bed
tadBED=hg19.GSE63525_GM12878.srt.bed
promoterLEN=1000
pre=COAD-US
matrixdir=/home1/04935/shaojf/stampede2/mutations/ICGC/lm/mut_exp_matrix

head -1 $matrixdir/$pre.bothWGS.mut.tsv > $pre.p.mut.tsv
cut -f 1 $matrixdir/$pre.bothWGS.mut.tsv | \
sed 's/:/\t/g' | sed -n '2,$p' | bedtools sort -i - | \
bedtools closest -D b -a - -b $tssBED | \
awk -F "\t" -v var=$promoterLEN -v OFS="\t" \
'$(NF-6)!="." && $NF > -1 * var && $NF < var {print $1":"$2":"$3,$7"~"$10}' | \
perl $bindir/relatedScripts/add_any_2files_together.pl $matrixdir/$pre.bothWGS.mut.tsv /dev/stdin 0 0 | cut -f 2,4- >> $pre.p.mut.tsv

Rscript $bindir/ICGC_mut2exp_lm.R $expMAT.WGS.sim $f.title $f &