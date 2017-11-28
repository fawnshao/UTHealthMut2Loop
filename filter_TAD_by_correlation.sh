#!/bin/sh
# filter mutation TAD list by expression correlation
tsvpre=$1
CorrelationFile=$2
# f=COAD-US.TAD_16-29607520-30607521.DO8360.tsv
for f in $tsvpre.TAD_*.tsv
do
	echo "pearson	p.value	spearman	p.value" > $f.tmp
	disease=`echo $f | awk -F"." '{print $1}'`
	tad=`echo $f | awk -F"." '{print $2}'`
	patient=`echo $f | awk -F"." '{print $3}'`
	echo "Disease	TAD	Patient" > $f.tmp1
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
		echo $disease"	"$tad"	"$patient"	" >> $f.tmp1
	done
	paste $f.tmp1 $f $f.tmp >> $tsvpre.allwithcor.tsv
	rm $f.tmp $f.tmp1
done
