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
motifDIR=motifdir
outpre=TAD.exp
id2name=
promoterMOTIF=/home1/04935/shaojf/stampede2/mutations/ICGC/p_motif/hg19.tss1k.motif.simcount

# usage function
function usage(){
   prog=`basename "$0"`
   cat << EOF

   Usage: $prog [-mut mutationTSV] [-tss tssBED] [-exp expMAT] [-tad tadBED] [-m motifBED] [-o outpre] [-plen promoterLEN] [-idconvert id2name] [-pmotif promoterMOTIF]

   optional arguments:
     -h            show this help message and exit
     -mut          mutation TSV files
     -tss          tss bed files
     -exp          expression matrix files with gene as row and sample as column
     -tad          TAD bed files
     -m            motif dir for motif files
     -o            output prefix
     -plen         up/down stream of TSS for promoter promoterSIZE
     -idconvert    ensembl ID to gene symbol file. do not set it if not needed
     -pmotif       promoter motif with counts
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
		-idconvert)
			id2name="$2"
			shift
			;;
		-pmotif)
			promoterMOTIF="$2"
			shift
			;;
		*)
			usage
			break
			;;
	esac
	shift
done

echo "mutationTSV  ="$mutationTSV
echo "expMAT       ="$expMAT
echo "tadBED       ="$tadBED
echo "tssBED       ="$tssBED
echo "motifDIR     ="$motifDIR
echo "outpre       ="$outpre
echo "promoterLEN  ="$promoterLEN
echo "id2name      ="$id2name
echo "promoterMOTIF="$promoterMOTIF

###############
# get the WGS mutations from simple_somatic_mutation.open.*.tsv
# gunzip -c $mutationTSV | awk -F"\t" -vOFS="\t" '$34=="WGS"{print $9,$10-1,$11,substr($7,1,15),".","+"}' | 
# uniq | bedtools sort -i - > $mutationTSV.WGS.srt.bed
echo +++++++++ get WGS mutations ++++++++
gunzip -c $mutationTSV | awk -F"\t" -vOFS="\t" '$34=="WGS"{print $9,$10-1,$11,$2,".","+"}' | 
uniq | bedtools sort -i - > $mutationTSV.WGS.srt.bed

# mutation to expressed promoter
# find the closet TSS for each mutation.
# ???? seems not WGS....
echo +++++++++ get promoter WGS mutations ++++++++
bedtools closest -D b -a $mutationTSV.WGS.srt.bed -b $tssBED | \
awk -F "\t" -v var=$promoterLEN -v OFS="\t" \
'$(NF-6)!="." && $NF > -1 * var && $NF < var {print $1,$2,$3,$4,$10,$13}' > ${outpre}.pMUT

# mutated promoter to TAD
# find the overlapped TAD for each mutated promoter.
echo +++++++++ get TAD promoter WGS mutations ++++++++
bedtools intersect -wao -a ${outpre}.pMUT -b $tadBED | awk '$NF > 0' | cut -f 1-10 > ${outpre}.pMUT.TAD

# promoter to TAD
# find all the promoters in a TAD, and keeps those with more than 2 genes
echo +++++++++ get TAd lists with many genes ++++++++
bedtools intersect -wao -a $tssBED -b $tadBED | awk '$NF > 0' | cut -f 1-10 > ${outpre}.p.TAD
cut -f10 ${outpre}.p.TAD | sort | uniq -c | awk '$1>1{print $2}' > ${outpre}.multiple.p.TAD
grep -wf ${outpre}.multiple.p.TAD ${outpre}.pMUT.TAD > ${outpre}.multiple.pMUT.TAD

# extend to TAD regions, and exclude all the sample with mutation in these extended regions
# the rest samples will be considered as a control
# and z scores are calulated as:
# Z-score (if comparing mt vs. wt) = [(value gene X in mt Y) - (mean gene X in wt)] / (standard deviation of gene X in wt)
echo +++++++++ get TAD mutations ++++++++
awk -v OFS="\t" '{a = $2 - 10000; b = $3 + 10000;}{if(a < 0){a = 0;}if(b < 0){ b = 0;} \
{print $1, a, b, $4}}' $tadBED > ${outpre}.extended.TAD
bedtools intersect -wao -a $mutationTSV.WGS.srt.bed -b ${outpre}.extended.TAD | \
awk '$NF > 0' > ${outpre}.extended.TAD.mut

# mutation to motif
# echo +++++++++ find if the mutation is in any motifs ++++++++
# bedtools intersect -wao -a $mutationTSV.WGS.srt.bed -b $motifBED | \
# awk '$NF > 0' > ${outpre}.WGSmut2motif

# find the TCGA id with WGS data and expression.
# extract the patient ID with WGS availble.
echo +++++++++ get WGS expressions ++++++++
cut -f 4 $mutationTSV.WGS.srt.bed | sort | uniq > ${outpre}.WGS.sampleid
# gunzip -c $expMAT | awk -F"\t" -vOFS="\t" '{print substr($5,1,15), $8, $9}' | \
# grep -f ${outpre}.WGS.sampleid > $expMAT.WGS.sim
gunzip -c $expMAT | awk -F"\t" -vOFS="\t" '{print $1, $8, $9}' | \
grep -wf ${outpre}.WGS.sampleid > $expMAT.WGS.sim

if [ ! -z "$id2name" ]
	then
	perl $bindir/relatedScripts/replace_ensembl_with_genesymbol.pl $id2name $expMAT.WGS.sim > $expMAT.WGS.sim.tmp
	mv $expMAT.WGS.sim.tmp $expMAT.WGS.sim
fi


# find the mutated promoter and with neigbors in the same TAD, 
# and the corresponding sample should have expression
# list the mutation with expression for the patient
cut -f1 $expMAT.WGS.sim | uniq | grep -wf - ${outpre}.multiple.pMUT.TAD > ${outpre}.multiple.pMUT.TAD.withexp

# TCGAsample="TCGA-AD-A5EJ-01"
# TADid="TAD_3"
echo +++++++++ test each mutated promoters in a TAD ++++++++
echo "TAD sample mutsample genes mutgene" > ${outpre}.IamGroot.Rinput
# count=0
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
			# count=$((count+1))
		done
	fi
done

count=`wc -l ${outpre}.IamGroot.Rinput | awk '{print $1}'`
echo +++++++++ Running Rscript to output loop translocate candidates  ++++++++
# use Z score to find the expression alteration direction in the TAD
if [ $count -lt 1000 ]
	then
	echo Running in one piece
	Rscript $bindir/ICGC_expression.R $expMAT.WGS.sim ${outpre}.IamGroot.Rinput ${outpre}
else
	echo Running in many pieces
	sed -n '2,$p' ${outpre}.IamGroot.Rinput | split -l 500 /dev/stdin ${outpre}.IamGroot.Rinput.
	# rm ${outpre}.IamGroot.Rinput.*.TAD.labeled.tsv
	np=`ls ${outpre}.IamGroot.Rinput.* | wc -l`
	for f in ${outpre}.IamGroot.Rinput.*
	do
		echo "TAD sample mutsample genes mutgene" | cat - $f > $f.title
		rm $f
		Rscript $bindir/ICGC_expression.R $expMAT.WGS.sim $f.title $f &
	done
	sleep 3m
	npfinished=`ls ${outpre}.IamGroot.Rinput.*.TAD.labeled.tsv | wc -l`
	while [ $npfinished -lt $np ]
	do
		sleep 30s
		npfinished=`ls ${outpre}.IamGroot.Rinput.*.TAD.labeled.tsv | wc -l`
	done
	echo Finishing all pieces
	cat ${outpre}.IamGroot.Rinput.*.TAD.labeled.tsv > ${outpre}.TAD.labeled.tsv
	for f in ${outpre}.IamGroot.Rinput.*.TAD_*.tsv
	do
		new=`echo $f | awk -F"." -vOFS="." '{print $1,$5,$6,$7}'`
		mv $f $new
	done
	rm ${outpre}.IamGroot.Rinput.*
fi

#####
## find if the mutation is associated with some motifs
echo +++++++++ checking if the mutation is in any motif  ++++++++
# echo "mutation-motif" > ${outpre}.LoopBroken.motif
echo -n "" > ${outpre}.LoopBroken.motif
echo -n "" > ${outpre}.LoopBroken.bed
grep LoopBroken ${outpre}.TAD.labeled.tsv | sed 's/"//g' | while read line
do
	tadid=`echo $line | awk '{print $2}'`
	sampleid=`echo $line | awk '{print $3}'`
	awk -v tad=$tadid -v sam=$sampleid -v OFS="\t" \
	'$10==tad && $4==sam {print $1,$2,$3,$10"~"$4"~"$5"~"$6}' \
	${outpre}.multiple.pMUT.TAD.withexp >> ${outpre}.LoopBroken.bed
done
for m in `ls $motifDIR/`
do
	bedtools intersect -wao -a ${outpre}.LoopBroken.bed -b $motifDIR/$m | \
	awk '$NF > 0' >> ${outpre}.LoopBroken.motif
done

rm ${outpre}.pMUT ${outpre}.p.TAD ${outpre}.pMUT.TAD
rm ${outpre}.multiple.p.TAD ${outpre}.multiple.pMUT.TAD # ${outpre}.multiple.pMUT.TAD.withexp
rm ${outpre}.extended.TAD ${outpre}.extended.TAD.mut
# rm ${outpre}.IamGroot.Rinput ${outpre}.WGS.sampleid

sh $bindir/ICGC_fimo_denovo_motif_in_LoopShiftMut.sh ${outpre}

# put all the results together
# ${outpre}.LoopBroken.motif
echo -n "disease	TAD	patient	Gene	" > ${outpre}.combined.tsv
echo -n "mut.mean	ctr.mean	ctr.sd	fc	outlier.flag	zscore.1	zscore.2	" >> ${outpre}.combined.tsv
echo -n "mut.exp	ctr.exp	mut.flag	alt.flag	" >> ${outpre}.combined.tsv
echo -n "promoter.motif.count	promoter.motifs	mutated.sites	mutated.motifs" >> ${outpre}.combined.tsv
echo "mut.motif.gain	mut.motif.lost" >> ${outpre}.combined.tsv
for f in ${outpre}.TAD_*.tsv
do
	disease=`echo $f | awk -F"." '{print $1}'`
	tad=`echo $f | awk -F"." '{print $2}'`
	patient=`echo $f | awk -F"." '{print $3}'`
	sed -n '2,$p' $f | while read line
	do
		gene=`echo $line | awk '{print $1}' | sed 's/"//g'`
		motifcount=`echo "" | awk -v a=$gene '{print "\t"a"|\n;"a"|"}' | grep -f - $promoterMOTIF | awk -F"\t" '{print $8}' |  tr ', ' '\n' | grep -v "^$" | sort | uniq | wc -l`
		motifs=`echo "" | awk -v a=$gene '{print "\t"a"|\n;"a"|"}' | grep -f - $promoterMOTIF | awk -F"\t" '{print $8}' |  tr ', ' '\n' | grep -v "^$" | sort | uniq | tr '\n' ',' | sed 's/,$//'`
		#echo $disease"	"$tad"	"$patient"	"$line"	"$motifcount"	"$motifs
		echo -n $disease"	"$tad"	"$patient"	" >> ${outpre}.combined.tsv
		echo -n $line | awk '{for(i=1;i<=NF;i++){printf "%s\t",$i}}' >> ${outpre}.combined.tsv
		echo -n $motifcount"	"$motifs"	" >> ${outpre}.combined.tsv
		if [[ `echo $line | grep MutatedPromoter` ]]; then
			position=`echo "" | awk -v a=$gene '{print "~"a"|\n;"a"|"}' | grep -f - ${outpre}.LoopBroken.bed | awk '{print $4}' | grep -w $tad | grep -w $patient | sort | uniq | tr '\n' '#' | sed 's/#$//'`
			mutmotif=`echo "" | awk -v a=$gene '{print "~"a"|\n;"a"|"}' | grep -f - ${outpre}.LoopBroken.motif | grep -w $tad | grep -w $patient | awk '{print $8}' | sort | uniq | tr '\n' ',' | sed 's/,$//'`
			tss=`echo "" | awk -v a=$gene '{print "~"a"|\n;"a"|"}' | grep -f - ${outpre}.LoopBroken.bed | awk '{print $4}' | grep -w $tad | grep -w $patient | sort | uniq | cut -d"~" -f3-4`
			gains=`echo $tss | grep -f - ${outpre}.LoopBroken.motif.gain | awk '{print $1}' | sort | uniq | tr '\n' ';' | sed 's/;$//'`
			losts=`echo $tss | grep -f - ${outpre}.LoopBroken.motif.lose | awk '{print $1}' | sort | uniq | tr '\n' ';' | sed 's/;$//'`
			echo -n $position"	"$mutmotif"	"$gains"	"$losts >> ${outpre}.combined.tsv
		else
			echo "" >> ${outpre}.combined.tsv
		fi
	done
done
