#!/bin/sh
# all input files should be 6-column sorted bed files
###
# the 4th column of mutation is named as: SampleName~PrimarySite~PrimaryHistology
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

echo "mutationBED="$mutationBED
echo "motifBED   ="$motifBED
echo "tssBED     ="$tssBED
echo "outpre     ="$outpre