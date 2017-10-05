#!/bin/sh
bindir=/home1/04935/shaojf/myTools/UTHealthMut2Loop
module load Rstats
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
     -plen         up/down stream of TSS for promoter promoterSIZE
     -mlen         flanking size for motif mutation.
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
awk -F "\t" '$(NF-6)!="."' | \
bedtools closest -D b -a - -b $tssBED | \
awk -F "\t" '$(NF-6)!="."' >${outpre}.closest

#########
# motif enrichment
#########
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
	awk -v var=$promoterLEN '$13 > -1 * var && $13 < var {print $4}' $prom_motif | \
	sort | uniq -c | awk -vOFS="\t" '{print $2,$1}' > $prom_motifCOUNT
fi

mut_motifCOUNT=$mutationBED.$motifBED".motif.count"
if [ ! -e "$mut_motifCOUNT" ]; then
	awk -v var=$motifFLANK -v OFS="\t" '$13 > -1 * var && $13 < var {print $7,$8,$9,$10}' ${outpre}.closest | \
	sort | uniq | cut -f 4 | sort | uniq -c | awk -vOFS="\t" '{print $2,$1}' > $mut_motifCOUNT
fi

prom_mut_motifCOUNT=$mutationBED.$motifBED".motif.prom.count"
if [ ! -e "$prom_mut_motifCOUNT" ]; then
	awk -v var1=$promoterLEN -v var2=$motifFLANK -v OFS="\t" \
	'$13 > -1 * var2 && $13 < var2 && $20 > -1 * var1 && $20 < var1 {print $7,$8,$9,$10}' ${outpre}.closest | \
	sort | uniq | cut -f 4 | sort | uniq -c | awk -vOFS="\t" '{print $2,$1}' > $prom_mut_motifCOUNT
fi

# hypergeometric test for motif enrichment
Rscript $bindir/HyperTest4MotifEnrichment.R \
	$mut_motifCOUNT $motifCOUNT ${outpre}.mut.enriched.motif
Rscript $bindir/HyperTest4MotifEnrichment.R \
	$prom_mut_motifCOUNT $prom_motifCOUNT ${outpre}.mut.enriched.motif.onlyprom

#########
# by disease
#########
mut2tss_by_disease=${mutationBED}.${tssBED}.distance.by.disease
if [ ! -e "$mut2tss_by_disease" ]; then
	cut -f 1-4,20 ${outpre}.closest | uniq | \
	cut -f 4-5 | sed 's/~/\t/' | cut -f 2-3 > $mut2tss_by_disease
fi

mut2motif_by_disease=${mutationBED}.${motifBED}.distance.by.disease
if [ ! -e "$mut2motif_by_disease" ]; then
	cut -f 1-4,13 ${outpre}.closest | uniq | \
	cut -f 4-5 | sed 's/~/\t/' | cut -f 2-3 > $mut2motif_by_disease
fi

motif2tss_by_motif=${motifBED}.${tssBED}.distance.by.motif
if [ ! -e "$motif2tss_by_motif" ]; then
	cut -f 1-4,13 $prom_motif | uniq | cut -f 4-5 > $motif2tss_by_motif
fi

motifmut2tss_by_disease=${mutationBED}.${motifBED}.${tssBED}.distance.by.disease
if [ ! -e "$motifmut2tss_by_disease" ]; then
	awk -v var=$motifFLANK -v OFS="\t" '$13 > -1 * var && $13 < var {print $1,$2,$3,$4,$20}' \
	${outpre}.closest | uniq | cut -f 4-5 | sed 's/~/\t/' | cut -f 2-3 > $motifmut2tss_by_disease
fi

Rscript $bindir/PlotDistance.R $mut2tss_by_disease $mut2motif_by_disease \
$motif2tss_by_motif $motifmut2tss_by_disease $outpre ${bindir}/ggplot2_multiple_plot.R




