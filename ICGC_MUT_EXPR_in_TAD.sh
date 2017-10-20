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
mutationTSV=simple_somatic_mutation.open.COAD-US.tsv
expMAT=exp_seq.COAD-US.tsv
outpre=TAD.exp


# usage function
function usage(){
   prog=`basename "$0"`
   cat << EOF

   Usage: $prog [-mut mutationTSV] [-tss tssBED] [-exp expMAT] [-tad tadBED] [-o outpre] [-plen promoterLEN]

   optional arguments:
     -h            show this help message and exit
     -mut          mutation TSV files
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
			mutationTSV="$2"
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
		*)
			usage
			break
			;;
	esac
	shift
done

echo "mutationTSV="$mutationTSV
echo "expMAT     ="$expMAT
echo "tadBED     ="$tadBED
echo "tssBED     ="$tssBED
echo "outpre     ="$outpre
echo "promoterLEN="$promoterLEN

###############
# get the WGS mutations from simple_somatic_mutation.open.*.tsv
# gunzip -c $mutationTSV | awk -F"\t" -vOFS="\t" '$34=="WGS"{print $9,$10-1,$11,substr($7,1,15),".","+"}' | 
# uniq | bedtools sort -i - > $mutationTSV.WGS.srt.bed
gunzip -c $mutationTSV | awk -F"\t" -vOFS="\t" '$34=="WGS"{print $9,$10-1,$11,$2,".","+"}' | 
uniq | bedtools sort -i - > $mutationTSV.WGS.srt.bed

# mutation to expressed promoter
# find the closet TSS for each mutation.
# ???? seems not WGS....
bedtools closest -D b -a $mutationTSV.WGS.srt.bed -b $tssBED | \
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

# extend to TAD regions, and exclude all the sample with mutation in these extended regions
# the rest samples will be considered as a control
# and z scores are calulated as:
# Z-score (if comparing mt vs. wt) = [(value gene X in mt Y) - (mean gene X in wt)] / (standard deviation of gene X in wt)
awk -v OFS="\t" '{a = $2 - 10000; b = $3 + 10000;}{if(a < 0){a = 0;}if(b < 0){ b = 0;} \
{print $1, a, b, $4}}' $tadBED > ${outpre}.extended.TAD
bedtools intersect -wao -a $mutationTSV.WGS.srt.bed -b ${outpre}.extended.TAD | \
awk '$NF > 0' > ${outpre}.extended.TAD.mut

# find the TCGA id with WGS data and expression.
# extract the patient ID with WGS availble.
cut -f 4 $mutationTSV.WGS.srt.bed | sort | uniq > ${outpre}.WGS.sampleid
# gunzip -c $expMAT | awk -F"\t" -vOFS="\t" '{print substr($5,1,15), $8, $9}' | \
# grep -f ${outpre}.WGS.sampleid > $expMAT.WGS.sim
gunzip -c $expMAT | awk -F"\t" -vOFS="\t" '{print $1, $8, $9}' | \
grep -wf ${outpre}.WGS.sampleid > $expMAT.WGS.sim

# find the mutated promoter and with neigbors in the same TAD, 
# and the corresponding sample should have expression
# list the mutation with expression for the patient
cut -f1 $expMAT.WGS.sim | uniq | grep -wf - ${outpre}.multiple.pMUT.TAD > ${outpre}.multiple.pMUT.TAD.withexp

# TCGAsample="TCGA-AD-A5EJ-01"
# TADid="TAD_3"
echo +++++++++ test each mutated promoters in a TAD ++++++++
echo "TAD sample mutsample genes mutgene" > ${outpre}.IamGroot.Rinput
cut -f 10 ${outpre}.multiple.pMUT.TAD.withexp | sort | uniq | while read TADid
do
	echo "|--TAD: "$TADid

	genes=`awk -v tad=$TADid '$10==tad{print $4}' ${outpre}.p.TAD | \
	tr ';' '\n' | cut -d "|" -f 1 | sort | uniq | tr '\n' ',' | sed 's/,$//'`

	# if [ "$genes" != "" ]
	if [ `echo $genes | grep ","` ]
		then
		# find samples with any (non-coding) mutations in the extended regions
		# and then excluded them
		mutsample=`awk -v tad=$TADid '$10==tad{print $4}' ${outpre}.extended.TAD.mut | \
		cut -d"~" -f1 | sort | uniq | tr '\n' ',' | sed 's/,$//'`

		awk -v tad=$TADid '$10==tad{print $4}' ${outpre}.multiple.pMUT.TAD.withexp | \
		cut -d"~" -f1 | sort | uniq | while read TCGAsample
		do
			echo "|----sample: "$TCGAsample
			mutp=`awk -v tad=$TADid -v sam=$TCGAsample '$10==tad && $4==sam {print $5}' \
			${outpre}.multiple.pMUT.TAD.withexp | sort | uniq | tr '\n' ',' | sed 's/,$//'`
			echo $TADid $TCGAsample $mutsample $genes $mutp >> ${outpre}.IamGroot.Rinput
		done
	fi
done

echo +++++++++ Running Rscript to output loop translocate candidates  ++++++++
# use Z score to find the expression alteration direction in the TAD
Rscript $bindir/ICGC_expression.R $expMAT.WGS.sim ${outpre}.IamGroot.Rinput ${outpre}

