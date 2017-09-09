#!/bin/sh
bindir=/home1/04935/shaojf/myTools/UTHealthMut2Loop
# all input files should be 6-column sorted bed files
###
# the 4th column of mutation is named as: SampleName~PrimarySite~PrimaryHistology
# default parameters
mutationBED=CosmicNCV.WGSv79.srt.bed
motifBED=homo_sapiens.GRCh38.motiffeatures.20161111.srt.bed
tssBED=clean_M2MT_hg38.refGene.tss.uniq.srt.bed
outpre=ensembl.motif.CosmicNCV.WGSv79
promoterLEN=1000
motifFLANK=10
###

# usage function
function usage()
{
   prog=`basename "$0"`
   cat << EOF

   Usage: $prog [-mut mutationBED] [-motif motifBED] [-tss tssBED] [-o outpre] [-plen promoterLEN] [-mlen motifFLANK]

   optional arguments:
     -h            show this help message and exit
     -mut          mutation bed files
     -motif        motif bed files
     -tss          tss bed files
     -o            output prefix
     -plen 		   up/down stream of TSS for promoter promoterSIZE
     -mlen		   flanking size for motif mutation.
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
		-plen)
			promoterLEN="$2"
			shift
			;;
		-mlen)
			motifFLANK="$2"
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
echo "promoterLEN="$promoterLEN
echo "motifFLANK ="$motifFLANK

# find the closet motif and TSS for each mutation.
bedtools closest -D b -a $mutationBED -b $motifBED | \
awk '$(NF-1)!="."' | \
bedtools closest -D b -a - -b $tssBED | \
awk '$(NF-1)!="."' >${outpre}.closest

# find the closet TSS for each motif.
prom_motif=$motifBED".tss"
if [ ! -e "$prom_motif" ]; then
	bedtools closest -D b -a $motifBED -b $tssBED | \
	awk '$(NF-1)!="."' > $prom_motif
fi

motifCOUNT=$motifBED".count"
if [ ! -e "$motifCOUNT" ]; then
	cut -f4 $motifBED | sort | uniq -c | \
	awk -vOFS="\t" '{print $2,$1}' > $motifCOUNT
fi

prom_motifCOUNT=$prom_motif".count"
if [ ! -e "$prom_motifCOUNT" ]; then
	awk -v var=promoterLEN '$13 > -1 * var && $13 < var {print $4}' $prom_motif | \
	sort | uniq -c | awk -vOFS="\t" '{print $2,$1}' > $prom_motifCOUNT
fi

mut_motifCOUNT=$mutationBED.$motifBED".motif.count"
if [ ! -e "$mut_motifCOUNT" ]; then
	awk -v var=motifFLANK '$13 > -1 * var && $13 < var {print $10}' ${outpre}.closest | \
	sort | uniq -c | awk -vOFS="\t" '{print $2,$1}' > $mut_motifCOUNT
fi

prom_mut_motifCOUNT=$mutationBED.$motifBED".motif.prom.count"
if [ ! -e "$prom_mut_motifCOUNT" ]; then
	awk -v var1=promoterLEN -v var2=motifFLANK \
	'$13 > -1 * var && $13 < var && $20 > -1 * var && $20 < var {print $10}' ${outpre}.closest | \
	sort | uniq -c | awk -vOFS="\t" '{print $2,$1}' > $prom_mut_motifCOUNT
fi

Rscript $bindir/HyperTest4MotifEnrichment.R \
	$mut_motifCOUNT $motifCOUNT ${outpre}.mut.enriched.motif
Rscript $bindir/HyperTest4MotifEnrichment.R \
	$prom_mut_motifCOUNT $prom_motifCOUNT ${outpre}.mut.enriched.motif.onlyprom

