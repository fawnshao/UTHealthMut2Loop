#!/bin/sh
# filter mutation TAD list by expression correlation
tsvpre=$1
CorrelationFile=$2
# f=COAD-US.TAD_16-29607520-30607521.DO8360.tsv
for f in $tsvpre.TAD_*.tsv
do
	echo "pearson	p.value	spearman	p.value" > $f.tmp
	mutgene=`awk -F"\t" '$11~/MutatedPromoter/{print $1}' $f | sed 's/"//g' | sed 's/^/-e /' | tr '\n' ' '`
	cat $f | sed -n '2,$p' | while read line
	do
		if [[ `echo $line | grep -v MutatedPromoter` ]]; then
			target=`echo $line | awk '{print $1}' | sed 's/"//g'`
			res=`grep -w $mutgene $CorrelationFile | grep $target | tr '\n' ';'`
			echo $res >> $f.tmp
		else
			echo "" >> $f.tmp
		fi
	done
	paste $f $f.tmp >> $tsvpre.allwithcor.tsv
	rm $f.tmp
done
