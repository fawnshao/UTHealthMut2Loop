#!/bin/sh
# all input files should be 6-column sorted bed files
###
# the 4th column of mutation is named as: SampleName~PrimarySite~PrimaryHistology
# default parameters
mutationBED=CosmicNCV.WGSv79.srt.bed
motifBED=homo_sapiens.GRCh38.motiffeatures.20161111.srt.bed
tssBED=clean_M2MT_hg38.refGene.tss.uniq.srt.bed
outpre=ensembl.motif.CosmicNCV.WGSv79
###

# usage function
function usage()
{
   prog=`basename "$0"`
   cat << EOF

   Usage: $prog [-mut mutationBED] [-motif motifBED] [-tss tssBED] [-o outpre] 

   optional arguments:
     -h            show this help message and exit
     -mut          mutation bed files
     -motif        motif bed files
     -tss          tss bed files
     -o            output prefix
EOF
}

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
		-motif)
			motifBED="$2"
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
		*)
			usage
			break
			;;
	esac
	shift
done

echo mutationBED=$mutationBED
echo motifBED=$motifBED
echo tssBED=$tssBED
echo outpre=$outpre

echo "mutationBED="$mutationBED
echo "motifBED   ="$motifBED
echo "tssBED     ="$tssBED
echo "outpre     ="$outpre


# find the closet motif and TSS for each mutation.
bedtools closest -D b -a $mutationBED -b $motifBED | \
awk '$(NF-1)!="."' | \
bedtools closest -D b -a - -b $tssBED | \
awk '$(NF-1)!="."' >$outpre.closest

motifCOUNT=$motifBED".count"
if [ ! -e "$motifCOUNT" ]; then
	cut -f4 $motifBED | sort | uniq -c | awk -vOFS="\t" '{print $2,$1}' > $motifCOUNT
fi
mutmotifCOUNT=$mutationBED".motif.count"
if [ ! -e "mutmotifCOUNT" ]; then
	cut -f4 $motifBED | sort | uniq -c | awk -vOFS="\t" '{print $2,$1}' > $motifCOUNT
fi
