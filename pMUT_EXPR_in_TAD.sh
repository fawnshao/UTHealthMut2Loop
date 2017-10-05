#!/bin/sh
bindir=/home1/04935/shaojf/myTools/UTHealthMut2Loop
module load Rstats
# all input files should be 6-column sorted bed files
###
# the 4th column of mutation is named as: SampleName~PrimarySite~PrimaryHistology
# the TAD file is 3-column bed files
# default parameters
#input
mutationBED=CosmicNCV.WGSv79.large_intestine.carcinoma.TCGA
tssBED=clean_M2MT_hg38.refGene.tss.uniq.srt.bed
expMAT=COAD.expr.txt
tadBED=hg38.GSE63525_GM12878.srt.bed
outpre=TAD.exp
promoterLEN=1000
foldchange=2
#output
expr_matrix=
promoter_mut_matrix=
gene_in_TAD_matrix=

# usage function
function usage(){
   prog=`basename "$0"`
   cat << EOF

   Usage: $prog [-mut mutationBED] [-tss tssBED] [-exp expMAT] [-tad tadBED] [-o outpre] [-plen promoterLEN] [-fc foldchange]

   optional arguments:
     -h            show this help message and exit
     -mut          mutation bed files
     -tss          tss bed files
     -exp          expression matrix files with gene as row and sample as column
     -tad          TAD bed files
     -o            output prefix
     -plen         up/down stream of TSS for promoter promoterSIZE
     -fc           fold change threhold for differentially expressed genes in the same TAD.
EOF
}

# get exact row and columns from a file
function get_row_column(){
	rowinput=$1
	colnum=$2
	matrixinput=$3
	matrixoutput=$4
	if [ -e "$matrixoutput" ]
	then
		rm $matrixoutput
	fi
	cat $rowinput | while read line
	do
		awk -v a=$line '$1==a' $matrixinput | cut -f $colnum >> $matrixoutput
	done
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
		-fc)
			foldchange="$2"
			shift
			;;
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
echo "foldchange ="$foldchange

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
# if [ -e "${outpre}.multiple.pMUT.TAD" ]
# then
# 	rm ${outpre}.multiple.pMUT.TAD
# fi
# cat ${outpre}.multiple.p.TAD | while read line
# do
# 	awk -v a=$line '$10==a' ${outpre}.pMUT.TAD >> ${outpre}.multiple.pMUT.TAD
# done

# find the mutated promoter and with neigbors in the same TAD, 
# and the corresponding sample should have expression
head -1 $expMAT | awk -F"\t" '{for(i=1;i<=NF;i++){print substr($i,1,15)}}' > ${expMAT}.id
cut -f 1 ${expMAT}.id | grep -wf - ${outpre}.multiple.pMUT.TAD > ${outpre}.multiple.pMUT.TAD.withexp

# TCGAsample="TCGA-AD-A5EJ-01"
# TADid="TAD_3"
echo +++++++++ test each mutated promoters in a TAD ++++++++
cut -f 10 ${outpre}.multiple.pMUT.TAD.withexp | sort | uniq | while read TADid
do
	echo "|--TAD: "$TADid
	awk -v tad=$TADid '$10==tad{print $4}' ${outpre}.p.TAD | \
	tr ';' '\n' | cut -d "|" -f 1 | sort | uniq > IamGroot.0
	
	awk -v tad=$TADid '$10==tad{print $4}' ${outpre}.multiple.pMUT.TAD.withexp | \
	cut -d"~" -f1 | sort | uniq | while read TCGAsample
	do
		echo "|----sample: "$TCGAsample
		col=`head -1 $expMAT | awk -F"\t" -v sample=$TCGAsample \
		'{for(i=1;i<=NF;i++){if(substr($i,1,15) == sample){print i}}}'`	
		get_row_column IamGroot.0 1 $expMAT ${outpre}.mutMAT
		# grep -wf IamGroot.0 $expMAT | cut -f 1 > ${outpre}.mutMAT
		cp ${outpre}.mutMAT ${outpre}.ctlMAT
		for i in $col
		do
			# grep -wf IamGroot.0 $expMAT | cut -f $i > IamGroot.$i
			get_row_column IamGroot.0 $i $expMAT IamGroot.$i
			paste ${outpre}.mutMAT IamGroot.$i > IamGroot.IamGroot
			mv IamGroot.IamGroot ${outpre}.mutMAT
			rm IamGroot.$i
		done

		awk -v tad=$TADid '$10==tad{print $4}' ${outpre}.pMUT.TAD | \
		cut -d"~" -f1 | sort | uniq > IamGroot.00
		head -1 $expMAT | awk -F"\t" '{{for(i=1;i<=NF;i++){o=substr($i,1,15);print i,o}}}' | \
		grep -vwf IamGroot.00 | awk '{print $1}' | sed -n '2,$p' > IamGroot.000
		# for i in `cat IamGroot.000`
		for i in `shuf IamGroot.000 | head -10`
		do
			# grep -wf IamGroot.0 $expMAT | cut -f $i > IamGroot.$i
			get_row_column IamGroot.0 $i $expMAT IamGroot.$i
			paste ${outpre}.ctlMAT IamGroot.$i > IamGroot.IamGroot
			mv IamGroot.IamGroot ${outpre}.ctlMAT
			rm IamGroot.$i
		done

		Rscript $bindir/Wilcox.test.R ${outpre}.ctlMAT ${outpre}.mutMAT $TADid $TCGAsample $foldchange
	done
done







