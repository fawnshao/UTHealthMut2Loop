#!/bin/sh
# rclone sync -L /home1/04935/shaojf/scratch/HiChIP.test/hkg2loops/ mygoogle:hkg_tsg/HiChIP.test/hkg2loops/
# args <- c("v2", "GTEx_sample.tissue.txt", "0.15", "1.5", "2", "0.25")
# outputpre <- args[1]
# tau.threshold <- as.numeric(args[3])
# sd.threshold <- as.numeric(args[4])
# median.threshold <- as.numeric(args[5])
# tau.threshold2 <- as.numeric(args[6])
myperl=/home1/04935/shaojf/myTools/BioinformaticsDaily/textProcess/add_any_2files_together.pl

for input in SRR5831489:GM12878_cell_line SRR5831492:K562_cell_line SRR5831493:K562_cell_line SRR5831496:Naive_T_primary_cell SRR5831505:Treg_primary_cell SRR5831507:Treg_primary_cell SRR5831509:Naive_primary_cell SRR5831511:HCASMC_cell_line
do
	echo $input
	acc=`echo $input | awk -F":" '{print $1}'`
	cell=`echo $input | awk -F":" '{print $2}'`
	f=$acc.filt.intra.loop_counts.bedpe
	awk '{print $1":"$2"-"$3"\t"$4":"$5"-"$6"\t"$8}' $acc.filt.intra.loop_counts.bedpe > $acc.filt.intra.loop_counts.bedpe.washU.txt
	awk -vOFS="\t" '$8 >= 5{print $1,$2,$3"\n"$4,$5,$6}' $f | bedtools sort -i - | uniq | awk -vOFS="\t" '{print $1,$2,$3,$1":"$2"-"$3}' | bedtools intersect -wo -a - -b <(awk -vOFS="\t" '{a=$2-1000;b=$2+1000;if($6=="-"){a=$3-1000;b=$3+1000}if(a<0){a=0;}print $1,a,b,$4,"1000",$6}' gencode.v19.gene.bed) > $f.gencode.txt
	perl $myperl <(cut -f 4,8 $f.gencode.txt) <(awk '$3 >= 5' $f.washU.txt) 0 0 | cut -f 1-3,5 | perl $myperl <(cut -f 4,8 $f.gencode.txt) /dev/stdin 0 1 | cut -f 1-4,6 > $f.anntotation.txt
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
