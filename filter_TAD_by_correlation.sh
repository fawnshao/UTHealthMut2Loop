#!/bin/sh
# filter mutation TAD list by expression correlation
tsvpre=$1
CorrelationFile=$2
# f=COAD-US.TAD_16-29607520-30607521.DO8360.tsv
echo -n "disease	TAD	patient	Gene	" > $tsvpre.allwithcor.tsv
echo -n "mut.mean	ctr.mean	ctr.sd	fc	outlier.flag	zscore.1	zscore.2	" >> $tsvpre.allwithcor.tsv
echo -n "mut.exp	ctr.exp	mut.flag	alt.flag	" >> $tsvpre.allwithcor.tsv
echo "pearson	p.value	spearman	p.value" >> $tsvpre.allwithcor.tsv
for f in $tsvpre.TAD_*.tsv
do
	disease=`echo $f | awk -F"." '{print $1}'`
	tad=`echo $f | awk -F"." '{print $2}'`
	patient=`echo $f | awk -F"." '{print $3}'`
	mutgene=`awk -F"\t" '$11~/MutatedPromoter/{print $1}' $f | sed 's/"//g' | sed 's/^/-e /' | tr '\n' ' '`
	cat $f | sed -n '2,$p' | while read line
	do
		echo -n $disease"	"$tad"	"$patient"	" >> $tsvpre.allwithcor.tsv
		echo -n $line | awk '{for(i=1;i<=NF;i++){printf "%s\t",$i}}' >> $tsvpre.allwithcor.tsv
		if [[ `echo $line | grep -v MutatedPromoter` ]]; then
			target=`echo $line | awk '{print $1}' | sed 's/"//g'`
			res=`grep -w $mutgene $CorrelationFile | grep $target | tr '\n' ';'`
			echo $res >> $tsvpre.allwithcor.tsv
		else
			echo "" >> $tsvpre.allwithcor.tsv
		fi
	done
done
