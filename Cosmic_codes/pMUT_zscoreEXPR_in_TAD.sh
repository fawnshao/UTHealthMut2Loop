#!/bin/sh
bindir=/home1/04935/shaojf/myTools/UTHealthMut2Loop
module load Rstats
# all input files should be 6-column sorted bed files
###
# the 4th column of mutation is named as: SampleName~PrimarySite.PrimaryHistology
# the TAD file is 3-column bed files
# default parameters
#input
mutationBED=CosmicNCV.TCGA.srt.bed
tssBED=clean_M2MT_hg38.refGene.tss.uniq.srt.bed
expMAT=CosmicCompleteGeneExpression.abnormal.tsv
tadBED=hg38.GSE63525_GM12878.srt.bed
outpre=TAD.exp
promoterLEN=1000

# usage function
function usage(){
   prog=`basename "$0"`
   cat << EOF

   Usage: $prog [-mut mutationBED] [-tss tssBED] [-exp expMAT] [-tad tadBED] [-o outpre] [-plen promoterLEN]

   optional arguments:
     -h            show this help message and exit
     -mut          mutation bed files
     -tss          tss bed files
     -exp          expression matrix files with gene as row and sample as column
     -tad          TAD bed files
     -o            output prefix
     -plen         up/down stream of TSS for promoter promoterSIZE
EOF
}

if [[ $# -eq 0 ]]
then
	usage
	exit
fi

while [[ $# -gt 0 ]]
do
	key="$1"
	case $key in
		-h)
			usage
			exit
			;;
		-mut)
			mutationBED="$2"
			shift
			;;
		-exp)
			expMAT="$2"
			shift
			;;
		-tad)
			tadBED="$2"
			shift
			;;
		-tss)
			tssBED="$2"
			shift
			;;
		-o)
			outpre="$2"
			shift
			;;
		-plen)
			promoterLEN="$2"
			shift
			;;
		# -fc)
		# 	foldchange="$2"
		# 	shift
		# 	;;
		*)
			usage
			break
			;;
	esac
	shift
done

echo "mutationBED="$mutationBED
echo "expMAT     ="$expMAT
echo "tadBED     ="$tadBED
echo "tssBED     ="$tssBED
echo "outpre     ="$outpre
echo "promoterLEN="$promoterLEN
# echo "foldchange ="$foldchange

###############
# mutation to expressed promoter
# find the closet TSS for each mutation.
bedtools closest -D b -a $mutationBED -b $tssBED | \
awk -F "\t" -v var=$promoterLEN -v OFS="\t" \
'$(NF-6)!="." && $NF > -1 * var && $NF < var {print $1,$2,$3,$4,$10,$13}' > ${outpre}.pMUT

# mutated promoter to TAD
# find the overlapped TAD for each mutated promoter.
bedtools intersect -wao -a ${outpre}.pMUT -b $tadBED | awk '$NF > 0' | cut -f 1-10 > ${outpre}.pMUT.TAD

# promoter to TAD
# find all the promoters in a TAD, and keeps those with more than 2 genes
bedtools intersect -wao -a $tssBED -b $tadBED | awk '$NF > 0' | cut -f 1-10 > ${outpre}.p.TAD
cut -f10 ${outpre}.p.TAD | sort | uniq -c | awk '$1>1{print $2}' > ${outpre}.multiple.p.TAD
grep -wf ${outpre}.multiple.p.TAD ${outpre}.pMUT.TAD > ${outpre}.multiple.pMUT.TAD

# TCGAsample="TCGA-A8-A09Z-01~breast.carcinoma.NS"
# TADid="TAD_517"
echo +++++++++ test each mutated promoters in a TAD ++++++++
echo "TAD TADgenes MUTsample MUTgene Flag" > ${outpre}.IamGroot
cut -f 10 ${outpre}.multiple.pMUT.TAD | sort | uniq | while read TADid
do
	echo "|--TAD: "$TADid

	genes=`awk -v tad=$TADid '$10==tad{print $4}' ${outpre}.p.TAD | \
	tr ';' '\n' | cut -d "|" -f 1 | sort | uniq | tr '\n' ',' | sed 's/,$//'`

	# if [ "$genes" != "" ]
	if [ `echo $genes | grep ","` ]
		then
		awk -v tad=$TADid '$10==tad{print $4}' ${outpre}.multiple.pMUT.TAD | \
		sort | uniq | while read TCGAsample
		do
			echo "|----sample: "$TCGAsample
			mutp=`awk -v tad=$TADid -v sam=$TCGAsample '$10==tad && $4==sam {print $5}' \
			${outpre}.multiple.pMUT.TAD | tr ';' '\n' | cut -d"|" -f1 | \
			sort | uniq | tr '\n' ',' | sed 's/,$//'`
			# mutsampleid=`echo $TCGAsample | cut -c1-15`
			# echo $mutp | cut -d"," | while read g
			# do
			# 	mutflag=`awk -v genename=$g -v samplename=$TCGAsample \
			# 	'$1==samplename && $2==genename {print $3}' $expMAT | tr '\n' ','`
			# done
			echo $TADid $genes $TCGAsample $mutp >> ${outpre}.IamGroot
		done
	fi

cut -d" " -f 3,4 ${outpre}.IamGroot | tr '~' ' ' | awk '{print $1"\t"$3}' | \
sort | uniq | grep -wf - $expMAT | cut -f1 | \
grep -wf - ${outpre}.IamGroot > ${outpre}.IamGroot.candidate

cat ${outpre}.IamGroot.candidate | tr '~' ' ' | cut -d" " -f3 | sort | uniq > ${outpre}.IamGroot.sample

cut -d" " -f2 ${outpre}.IamGroot.candidate | tr ',' '\n' | sort | \
uniq | grep -wf - $expMAT | grep -wf ${outpre}.IamGroot.sample > ${outpre}.IamGroot.expression


echo +++++++++ Running Rscript to output loop translocate candidates  ++++++++
# Rscript $bindir/read.in.expression.matrix.R $expMAT IamGroot.Rinput $foldchange
# use Z score to find the expression alteration direction in the TAD
Rscript $bindir/opposite_expr.R $expMAT ${outpre}.IamGroot.Rinput ${outpre}.mutated.sampleid ${outpre}






