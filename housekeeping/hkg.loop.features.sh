#!/bin/sh
# http://epigenomegateway.wustl.edu/browser/?genome=hg19&session=1zqrsNwUIo
# rclone sync -L /home1/04935/shaojf/scratch/HiChIP.test/hkg2loops/ mygoogle:hkg_tsg/HiChIP.test/hkg2loops/
# bedtools sort *.bed
# bgzip *.bed
# tabix -p bed *.bed.gz

# for f in *.filt.intra.loop_counts.bedpe.washU.txt
# do
# 	acc=`echo $f | sed 's/.filt.intra.loop_counts.bedpe.washU.txt//'`
# 	awk -vOFS="\t" '{print $1,$2","$3,NR"_1\t.\n"$2,$1","$3,NR"_2\t."}' $f | sed 's/:/\t/;s/-/\t/' | bedtools sort -i - > $acc.interaction.txt
# 	bgzip $acc.interaction.txt
# 	tabix -p bed $acc.interaction.txt.gz
# done

# args <- c("v2", "GTEx_sample.tissue.txt", "0.15", "1.5", "2", "0.25")
# outputpre <- args[1]
# tau.threshold <- as.numeric(args[3])
# sd.threshold <- as.numeric(args[4])
# median.threshold <- as.numeric(args[5])
# tau.threshold2 <- as.numeric(args[6])
myperl=/home1/04935/shaojf/myTools/BioinformaticsDaily/textProcess/add_any_2files_together.pl
mymatrix=/home1/04935/shaojf/myTools/BioinformaticsDaily/textProcess/make_matrixwith3col_from_single_file.pl
mymerge=/home1/04935/shaojf/myTools/BioinformaticsDaily/textProcess/merge.rows.by.specificID.for.single.file.pl
analyzenetwork=/home1/04935/shaojf/myTools/UTHealthMut2Loop/housekeeping/analyze.looping.network.pl

# only these two are CTCF HiChIP:
# SRR5831509	SAMN07357487	GSM2705061	OTHER	152	SRX3008815	Illumina HiSeq 4000	45,947	18,503	Homo sapiens	male	--	Sorted primary Naive T cells	Naïve T primary cell
# SRR5831508	SAMN07357442	GSM2705060	OTHER	152	SRX3008814	Illumina HiSeq 4000	52,887	21,407	Homo sapiens	male	--	Sorted primary Naive T cells	Naïve T primary cell

for input in SRR5831489:GM12878_cell_line SRR5831490:GM12878_cell_line SRR5831492:K562_cell_line SRR5831493:K562_cell_line SRR5831494:My-La_cell_line SRR5831495:My-La_cell_line SRR5831496:Naive_T_primary_cell SRR5831497:Naive_primary_cell SRR5831498:Naive_primary_cell SRR5831499:Naive_primary_cell SRR5831501:Th17_primary_cell SRR5831503:Th17_primary_cell SRR5831504:Th17_primary_cell SRR5831505:Treg_primary_cell SRR5831505:Treg_primary_cell SRR5831507:Treg_primary_cell SRR5831509:Naive_primary_cell SRR5831511:HCASMC_cell_line SRR5831512:HCASMC_cell_line
do
	echo $input
	acc=`echo $input | awk -F":" '{print $1}'`
	cell=`echo $input | awk -F":" '{print $2}'`
	f=$acc.filt.intra.loop_counts.bedpe
	awk '{print $1":"$2"-"$3"\t"$4":"$5"-"$6"\t"$8}' $acc.filt.intra.loop_counts.bedpe > $acc.filt.intra.loop_counts.bedpe.washU.txt
	# awk -vOFS="\t" '$8 >= 5{print $1,$2,$3"\n"$4,$5,$6}' $f | bedtools sort -i - | uniq | awk -vOFS="\t" '{print $1,$2,$3,$1":"$2"-"$3}' | bedtools intersect -wo -a - -b <(awk -vOFS="\t" '{a=$2-1000;b=$2+1000;if($6=="-"){a=$3-1000;b=$3+1000}if(a<0){a=0;}print $1,a,b,$4,"1000",$6}' gencode.v19.gene.bed) > $f.gencode.txt
	###### $8 >= 5
	awk -vOFS="\t" '$8 >= 10{print $1,$2,$3"\n"$4,$5,$6}' $f | bedtools sort -i - | uniq | awk -vOFS="\t" '{print $1,$2,$3,$1":"$2"-"$3}' | bedtools intersect -wo -a - -b <(awk -vOFS="\t" '{a=$2-5000;b=$2+5000;if($6=="-"){a=$3-5000;b=$3+5000}if(a<0){a=0;}print $1,a,b,$4,"5000",$6}' gencode.v19.gene.bed) > $f.gencode.txt
	###### $3 >= 5
	perl $myperl <(cut -f 4,8 $f.gencode.txt) <(awk '$3 >= 10' $f.washU.txt) 0 0 | cut -f 1-3,5 | perl $myperl <(cut -f 4,8 $f.gencode.txt) /dev/stdin 0 1 | cut -f 1-4,6 > $f.anntotation.txt
	perl $myperl v2.hkg.tsg.anntotation.txt $f.anntotation.txt 0 3 | cut -f 1-5,7-8 | perl $myperl v2.hkg.tsg.anntotation.txt /dev/stdin 0 4 | cut -f 1-7,9-10 > $f.anntotation.hkg.tsg
	awk -F"\t" -vOFS="\t" '$4!="/" && $5!="/" {print $4,$5,$6,$7,$8,$9}' $f.anntotation.hkg.tsg > $f.pp.txt
	awk -F"\t" -vOFS="\t" '{if($4!="/" && $5=="/") {print $4,$2,$6,$7,$8,$9}}' $f.anntotation.hkg.tsg > $f.pe.txt
	awk -F"\t" -vOFS="\t" '{if($4=="/" && $5!="/") {print $5,$1,$8,$9,$6,$7}}' $f.anntotation.hkg.tsg >> $f.pe.txt

	echo "Class Total promoter-promoter.count promoter-enhancer.count promoter-promoter% promoter-enhancer%"
	for pre in HKG.1 HKG.2 Brain Liver Spleen EBV Testis
	do
		pcount=`awk -vvar=$pre -F"\t" '{if($3~var){print $1}if($4~var){print $2}}' $f.pp.txt | sort | uniq | wc -l`
		ecount=`cut -f 1,3 $f.pe.txt | awk -F"\t" '$2!="/"' | sort | uniq | grep $pre | wc -l`
		all=`grep $pre v2.hkg.tsg.anntotation.txt | wc -l`
		pperc=`echo $pcount/$all | bc -l`
		eperc=`echo $ecount/$all | bc -l`
		echo $pre $all $pcount $ecount $pperc $eperc
	done
	pcount=`awk -F"\t" '{print $1"\n"$2}' $f.pp.txt | sort | uniq | wc -l`
	ecount=`cut -f 1 $f.pe.txt | sort | uniq | wc -l`
	all=`cat gencode.v19.gene.bed | wc -l`
	pperc=`echo $pcount/$all | bc -l`
	eperc=`echo $ecount/$all | bc -l`
	echo gencode $all $pcount $ecount $pperc $eperc

	pre=protein_coding
	pcount=`awk -vvar=$pre -F"\t" '{if($4==var){print $1}if($6==var){print $2}}' $f.pp.txt | sort | uniq | wc -l`
	ecount=`cut -f 1,4 $f.pe.txt | awk -F"\t" '$2!="/"' | sort | uniq | grep $pre | wc -l`
	all=`grep protein_coding gencode.v19.gene.anntotation.txt | wc -l`
	pperc=`echo $pcount/$all | bc -l`
	eperc=`echo $ecount/$all | bc -l`
	echo $pre $all $pcount $ecount $pperc $eperc

	pcount=`awk -F"\t" '{print $1"\n"$2}' $f.pp.txt | sort | uniq | perl $myperl <(awk '$22>2{print $1}' v2.log2tpm.median.tsv) /dev/stdin 0 0 | awk '$2!="/"' | wc -l`
	ecount=`cut -f 1 $f.pe.txt | sort | uniq | perl $myperl <(awk '$22>2{print $1}' v2.log2tpm.median.tsv) /dev/stdin 0 0 | awk '$2!="/"' | wc -l`
	all=`awk '$22>2{print $1}' v2.log2tpm.median.tsv | wc -l`
	pperc=`echo $pcount/$all | bc -l`
	eperc=`echo $ecount/$all | bc -l`
	echo EBV.expressed $all $pcount $ecount $pperc $eperc
done
# more 5k.stats.txt | awk '{print $5"\t"$6}' | sed -n '3~12p' | tr "\n" "\t"  | awk '{print $0}' > x
# more 5k.stats.txt | awk '{print $5"\t"$6}' | sed -n '4~12p' | tr "\n" "\t"  | awk '{print $0}' >> x
# more 5k.stats.txt | awk '{print $5"\t"$6}' | sed -n '5~12p' | tr "\n" "\t"  | awk '{print $0}' >> x
# more 5k.stats.txt | awk '{print $5"\t"$6}' | sed -n '6~12p' | tr "\n" "\t"  | awk '{print $0}' >> x
# more 5k.stats.txt | awk '{print $5"\t"$6}' | sed -n '7~12p' | tr "\n" "\t"  | awk '{print $0}' >> x
# more 5k.stats.txt | awk '{print $5"\t"$6}' | sed -n '8~12p' | tr "\n" "\t"  | awk '{print $0}' >> x
# more 5k.stats.txt | awk '{print $5"\t"$6}' | sed -n '9~12p' | tr "\n" "\t"  | awk '{print $0}' >> x
# more 5k.stats.txt | awk '{print $5"\t"$6}' | sed -n '10~12p' | tr "\n" "\t"  | awk '{print $0}' >> x
# more 5k.stats.txt | awk '{print $5"\t"$6}' | sed -n '11~12p' | tr "\n" "\t"  | awk '{print $0}' >> x
# more 5k.stats.txt | awk '{print $5"\t"$6}' | sed -n '12~12p' | tr "\n" "\t"  | awk '{print $0}' >> x

echo "gene partner PETcount experiment" | tr " " "\t" > gencode.looping.txt
for f in *filt.intra.loop_counts.bedpe.anntotation.txt
do
	acc=`echo $f | awk -F"." '{print $1}'`
	awk -v OFS="\t" -v var=$acc '{if($4=="/" && $5!="/"){print $5,$1,$3,var}if($4!="/" && $5=="/"){print $4,$2,$3,var}if($4!="/" && $5!="/"){print $4,$5,$3,var}}' $f >> gencode.looping.txt
done
cut -f2 gencode.looping.txt | grep chr | sort | uniq | awk '{print $0"\t"$0}' | sed 's/:/\t/;s/-/\t/' | bedtools sort -i > candidate.enhancer.bed
bedtools merge -d 1000 -i candidate.enhancer.bed | awk -vOFS="\t" '{print $0"\t"$1":"$2"-"$3}' > candidate.enhancer.merged.bed 
# bedtools intersect -wo -a candidate.enhancer.merged.bed -b <(awk -vOFS="\t" '{a=$2-5000;b=$2+5000;if($6=="-"){a=$3-5000;b=$3+5000}if(a<0){a=0;}print $1,a,b,$4,"5000",$6}' gencode.v19.gene.bed) > candidate.enhancer.merged.gencodeTSS
bedtools intersect -wo -a candidate.enhancer.merged.bed -b <(awk -vOFS="\t" '{a=$2-1000;b=$2+1000;if($6=="-"){a=$3-1000;b=$3+1000}if(a<0){a=0;}print $1,a,b,$4,"1000",$6}' gencode.v19.all.transcript.tss.bed) > candidate.enhancer.merged.gencodeTSS
bedtools intersect -wo -a candidate.enhancer.merged.bed -b <(awk -vOFS="\t" '{a=$2-1000;b=$2+1000;if(a<0){a=0;}print "chr"$1,a,b,$4,"1000",$6}' hg19.refGene.tss.uniq.srt.bed) > candidate.enhancer.merged.refseqTSS
perl $mymerge <(cat candidate.enhancer.merged.gencodeTSS candidate.enhancer.merged.refseqTSS | cut -f 4,8) > candidate.enhancer.merged.map2tss
# bedtools intersect -wa -v -a candidate.enhancer.merged.bed -b <(cat gencode.v19.gene.bed hg19.RefSeq.gene.bed | awk -vOFS="\t" '{a=$2-5000;b=$2+5000;if($6=="-"){a=$3-5000;b=$3+5000}if(a<0){a=0;}print $1,a,b,$4,"5000",$6}') > candidate.enhancer.merged.notTSS.bed
bedtools intersect -wao -a candidate.enhancer.bed -b candidate.enhancer.merged.bed | cut -f 4,8 > candidate.enhancer.mapping.txt.tmp
paste <(cut -f 1 candidate.enhancer.merged.map2tss) <(cut -f 2 candidate.enhancer.merged.map2tss | cut -d":" -f 2)| perl $myperl /dev/stdin candidate.enhancer.mapping.txt.tmp 0 1 > candidate.enhancer.mapping.txt.tmp.1
awk -vOFS="\t" '{if($4!="/"){print $1,$4}else{print $1,$2}}' candidate.enhancer.mapping.txt.tmp.1 > candidate.enhancer.mapping.txt
rm candidate.enhancer.mapping.txt.tmp candidate.enhancer.mapping.txt.tmp.1
# bedtools intersect -wo -a candidate.enhancer.bed -b candidate.enhancer.merged.notTSS.bed | cut -f 4,8 > candidate.enhancer.mapping.txt
# perl $mymerge <(cat candidate.enhancer.merged.gencodeTSS candidate.enhancer.merged.refseqTSS | cut -f 4,8) >> candidate.enhancer.mapping.txt

grep -v chr gencode.looping.txt > gencode.looping.merged.txt
perl $myperl candidate.enhancer.mapping.txt <(grep chr gencode.looping.txt) 0 1 | awk -vOFS="\t" '{print $1,$6,$3,$4}' >> gencode.looping.merged.txt

for exp in `cut -f 4 gencode.looping.merged.txt | tail -n +2 | uniq | sort | uniq`
do
    awk -v var=$exp '$3>=10 && $4==var && $1!=$2{print $1"\t"$2}' gencode.looping.merged.txt | sort | uniq > interaction.$exp
    perl $analyzenetwork interaction.$exp
done

perl $mymatrix <(awk -vOFS="\t" '{print $1"%"$2,$4,$3}' gencode.looping.merged.txt | tail -n +2) > gencode.looping.merged.mat
# grep chr gencode.looping.merged.mat | 
head -1 gencode.looping.merged.mat | cut -f 2- | awk '{print "ID1%ID2\t"$0"\tExperimentCount\tID1.anntotation\tID2.anntotation"}' > gencode.looping.merged.mat.anntotation
#################### needs to set the column number
# awk '{sum=0;for(i=2;i<=NF;i++){if($i>0){sum+=1}}print $0"\t"sum}' <(tail -n +2 gencode.looping.merged.mat) | sed 's/%/\t/' | perl $myperl v2.hkg.tsg.anntotation.txt /dev/stdin 0 0 | cut -f 1-12,14 | perl $myperl v2.hkg.tsg.anntotation.txt /dev/stdin 0 1 | cut -f 1-13,15 | sed 's/\t/%/' >> gencode.looping.merged.mat.anntotation
awk '{sum=0;for(i=2;i<=NF;i++){if($i>9){sum+=1}}print $0"\t"sum}' <(tail -n +2 gencode.looping.merged.mat) | sed 's/%/\t/' | perl $myperl v2.hkg.tsg.anntotation.txt /dev/stdin 0 0 | cut -f 1-21,23 | perl $myperl v2.hkg.tsg.anntotation.txt /dev/stdin 0 1 | cut -f 1-22,24 | sed 's/\t/%/' >> gencode.looping.merged.mat.anntotation
####################

# head -1 gencode.looping.merged.mat.anntotation > gencode.looping.merged.mat.hkg.tsg.anntotation
awk -F"\t" '$(NF-1)!="/" || $NF!="/"' gencode.looping.merged.mat.anntotation > gencode.looping.merged.mat.hkg.tsg.anntotation
awk -F"\t" '$(NF-2) >= 8' gencode.looping.merged.mat.hkg.tsg.anntotation > hkg.tsg.commonloops.txt
awk '$1~/chr/{print $1}' hkg.tsg.commonloops.txt | awk -F"%" '{print $2"\t"$2}' | sed 's/:/\t/;s/-/\t/' | bedtools sort -i - | uniq > hkg.tsg.commonenhancer.bed
bedtools intersect -wo -a hkg.tsg.commonenhancer.bed -b hg19.tissue.roadmap.activeenhancer.bed  > hkg.tsg.commonenhancer.chromHMM.txt
perl $mymerge <(cut -f 4,9 hkg.tsg.commonenhancer.chromHMM.txt) > hkg.tsg.commonenhancer.map2chromHMM.txt
# split -l 500000 hg19.tissue.roadmap.25_imputed12marks.bed tissue.roadmap.
# for f in tissue.roadmap.??
# do
# 	bedtools intersect -wo -a hkg.tsg.commonenhancer.bed -b $f >> hkg.tsg.commonenhancer.chromHMM.txt
# done

# cut -f 1,11- gencode.looping.merged.mat.hkg.tsg.anntotation | awk '$2 > 3 && $3~/HKG/ && $4~/HKG/{print $1}' | sort | uniq | wc -l
# cut -f 1,11- gencode.looping.merged.mat.hkg.tsg.anntotation | awk '$2 > 3 && $3~/HKG/ && $1~/chr/{print $1}' | sort | uniq | wc -l
#################### do not need to set the column number
awk -F"\t" '$(NF-2) >= 10 && $0~/HKG/ && $1!~/chr/' gencode.looping.merged.mat.hkg.tsg.anntotation > common.hkg2promoter.txt
awk -F"\t" '$(NF-2) >= 10 && $0~/HKG/ && $1~/chr/' gencode.looping.merged.mat.hkg.tsg.anntotation > common.hkg2enhancer.txt
awk -F"\t" '$(NF-2) >= 10 && (($NF!="/" && $NF!~/HKG/) || ($(NF-1)!="/" && $(NF-1)!~/HKG/)) && $1!~/chr/' gencode.looping.merged.mat.hkg.tsg.anntotation > common.tsg2promoter.txt
awk -F"\t" '$(NF-2) >= 10 && (($NF!="/" && $NF!~/HKG/) || ($(NF-1)!="/" && $(NF-1)!~/HKG/)) && $1~/chr/' gencode.looping.merged.mat.hkg.tsg.anntotation > common.tsg2enhancer.txt
####################
# ENSG00000011052.17|NME2 does not have expression in GTEx data
cut -f 1 common.hkg2promoter.txt | sed 's/%/\t/' | awk '$1~/ENSG/ && $2~/ENSG/ && $1!=$2 && $2!~/;/' > common.hkg2promoter.sim.pairs
cut -f 1 common.tsg2promoter.txt | sed 's/%/\t/' | awk '$1~/ENSG/ && $2~/ENSG/ && $1!=$2 && $2!~/;/' > common.tsg2promoter.sim.pairs
perl $myperl GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct <(cut -f 1 common.hkg2promoter.sim.pairs) 0 0 | cut -f 2- > tmp.1
perl $myperl GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct <(cut -f 2 common.hkg2promoter.sim.pairs) 0 0 | cut -f 2- > tmp.2
Rscript ~/myTools/UTHealthMut2Loop/housekeeping/calc_cor_gene_expr_for_GTEx.R tmp.1 tmp.2 common.hkg2promoter.sim.pairs.cor.tsv
perl $myperl GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct <(cut -f 1 common.tsg2promoter.sim.pairs) 0 0 | cut -f 2- > tmp.1
perl $myperl GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct <(cut -f 2 common.tsg2promoter.sim.pairs) 0 0 | cut -f 2- > tmp.2
Rscript ~/myTools/UTHealthMut2Loop/housekeeping/calc_cor_gene_expr_for_GTEx.R tmp.1 tmp.2 common.tsg2promoter.sim.pairs.cor.tsv
rm tmp.1 tmp.2

#################### annotate all the loopping partner of HKGs
awk -F"\t" '($NF!="/" && $NF!~/HKG/) || ($(NF-1)!="/" && $(NF-1)!~/HKG/)' gencode.looping.merged.mat.hkg.tsg.anntotation > all.tsg.partner.txt
head -1 all.tsg.partner.txt > all.hkg.partner.txt
grep HKG gencode.looping.merged.mat.hkg.tsg.anntotation >> all.hkg.partner.txt

awk -F"\t" '$(NF-2) >= 12{print $1"\t"$(NF-2)}' all.hkg.partner.txt | sed 's/%/\t/' > all.hkg.gt12.interaction.txt
awk -F"\t" '$(NF-2) >= 12{print $1"\t"$(NF-1)"\t"$NF}' all.hkg.partner.txt | sed 's/%/\t/' | awk -F"\t" -vOFS="\t" '{print $1,$3"\n"$2,$4}' | tail -n +2 | sort | uniq | awk 'BEGIN{print "ID\tannotation"}{print $0}'> all.hkg.gt12.node.txt

