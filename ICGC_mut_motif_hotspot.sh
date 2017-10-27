#!/bin/sh
bindir=/home1/04935/shaojf/myTools/UTHealthMut2Loop
module load Rstats
# TSS should be 6-column sorted bed files
# the TAD file is 4-column bed files
# default parameters
#input
tssBED=hg19.refGene.tss.uniq.srt.bed
promoterLEN=1000
mutationTSV=simple_somatic_mutation.open.COAD-US.tsv
motifDIR=motifdir
outpre=TAD.exp


# usage function
function usage(){
   prog=`basename "$0"`
   cat << EOF

   Usage: $prog [-mut mutationTSV] [-tss tssBED] [-m motifBED] [-o outpre] [-plen promoterLEN]

   optional arguments:
     -h            show this help message and exit
     -mut          mutation TSV files
     -tss          tss bed files
     -m            motif dir for motif files
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
		-tss)
			tssBED="$2"
			shift
			;;
		-m)
			motifDIR="$2"
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
echo "tssBED     ="$tssBED
echo "motifDIR   ="$motifDIR
echo "outpre     ="$outpre
echo "promoterLEN="$promoterLEN

###############
# get the WGS mutations from simple_somatic_mutation.open.*.tsv
# gunzip -c $mutationTSV | awk -F"\t" -vOFS="\t" '$34=="WGS"{print $9,$10-1,$11,substr($7,1,15),".","+"}' | 
# uniq | bedtools sort -i - > $mutationTSV.WGS.srt.bed
######### skip for temporary
# echo +++++++++ get WGS mutations ++++++++
# gunzip -c $mutationTSV | awk -F"\t" -vOFS="\t" '$34=="WGS"{print $9,$10-1,$11,$2,".","+"}' | 
# uniq | bedtools sort -i - > $mutationTSV.WGS.srt.bed

# mutation to expressed promoter
# find the closet TSS for each mutation.
# ???? seems not WGS....
######### skip for temporary
# echo +++++++++ get promoter WGS mutations ++++++++
# bedtools closest -D b -a $mutationTSV.WGS.srt.bed -b $tssBED | \
# awk -F "\t" -v var=$promoterLEN -v OFS="\t" \
# '$(NF-6)!="." && $NF > -1 * var && $NF < var {print $1,$2,$3,$4,$10,$13}' > ${outpre}.pMUT

# mutation to motif
# echo +++++++++ find if the mutation is in any motifs ++++++++
# bedtools intersect -wao -a $mutationTSV.WGS.srt.bed -b $motifBED | \
# awk '$NF > 0' > ${outpre}.WGSmut2motif
echo -n "" > ${outpre}.WGS.motif
echo -n "" > ${outpre}.WGS.p.motif
for m in `ls $motifDIR/`
do
	bedtools intersect -wao -a $mutationTSV.WGS.srt.bed -b $motifDIR/$m | \
	awk -v lab=$m '$NF > 0 {print $0"\t"lab}' >> ${outpre}.WGS.motif
	bedtools intersect -wao -a ${outpre}.pMUT -b $motifDIR/$m | \
	awk -v lab=$m '$NF > 0 {print $0"\t"lab}' >> ${outpre}.WGS.p.motif
done


