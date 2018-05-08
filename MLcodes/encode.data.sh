#!/bin/sh
myperl=/home1/04935/shaojf/myScripts/add_any_2files_together.pl
# rclone sync -L /home1/04935/shaojf/scratch/ENCODE_res/formated.files mygoogle:hkg_tsg/ENCODE_res/formated.files
### paste tpm value into one matrix 
cut -f 1-2 A172__transcript.quantifications.ENCFF185SEV.tsv | awk '{print $2"|"$1}' > encode.transcript.tsv
ls *transcript*.tsv | grep -v "dendritic.cell" | while read files
do
	echo $files > a
	cut -f 6 $files | tail -n +2 >> a
	paste encode.transcript.tsv a > b
	mv b encode.transcript.tsv
	rm a
done

### prepare the tss for each transcript
# gunzip -c gencode.v19.annotation.gtf.gz | awk '$3=="transcript"' | awk -vOFS="\t" '{print $1,$4,$5,$10"|"$12,".",$7}' | sed 's/"//g;s/;//g' > gencode.v19.transcript.bed
# login2.stampede2(1052)$ gunzip -c gencode.v19.annotation.gtf.gz | awk '$3=="transcript"' | grep "appris_principal" | grep "level 1" | wc -l
# 464
# login2.stampede2(1053)$ gunzip -c gencode.v19.annotation.gtf.gz | awk '$3=="transcript"' | grep "appris_principal" | grep -e "level 1" -e "level 2"| wc -l
# 21269
# login2.stampede2(1021)$ gunzip -c gencode.v19.annotation.gtf.gz | awk '$3=="transcript"' | grep -w -e "appris_principal" -e "appris_candidate" | wc -l
# 28036
# gunzip -c gencode.v19.annotation.gtf.gz | awk '$3=="transcript"' | grep "appris_principal" | grep -e "level 1" -e "level 2" | awk -vOFS="\t" '{print $1,$4,$5,$10"|"$12,".",$7}' | sed 's/"//g;s/;//g' > gencode.v19.transcript.appris_principal.level12.bed
# ln -s gencode.v19.transcript.appris_principal.level12.bed gencode.v19.transcript.bed
gunzip -c gencode.v19.annotation.gtf.gz | awk '$3=="transcript"' | grep -w -e "appris_principal" -e "appris_candidate" | awk -vOFS="\t" '{print $1,$4,$5,$10"|"$12,".",$7}' | sed 's/"//g;s/;//g' > gencode.v19.transcript.appris.bed
ln -s gencode.v19.transcript.appris.bed gencode.v19.transcript.bed
perl $myperl <(sed 's/|/\t/' gencode.v19.transcript.bed) <(sed 's/|/\t/' hkg.tsg.srtbyPCA.class) 3 0 | awk -F"\t" -vOFS="\t" '{print $1"|"$2,$3,$4,$5,$6,$7"|"$8,$9,$10}' | tail -n +2 | grep "/" > hkg.tsg.srtbyPCA.no.transcript
perl $myperl <(sed 's/|/\t/' gencode.v19.transcript.bed) <(sed 's/|/\t/' hkg.tsg.srtbyPCA.class) 3 0 | awk -F"\t" -vOFS="\t" '{print $1"|"$2,$3,$4,$5,$6,$7"|"$8,$9,$10}' | tail -n +2 | grep -v "/" > hkg.tsg.srtbyPCA.transcript.bed
cut -f 1-2,6 hkg.tsg.srtbyPCA.transcript.bed > hkg.tsg.srtbyPCA.transcript.class

### find the overlapped features for transcript
bedtools intersect -wo -a hg19.cpgIslandExt.withCount.bed -b <(awk -vOFS="\t" '{a=$2-300;b=$2+100;if($6=="-"){a=$3-100;b=$3+300}if(a<0){a=0;}print $1,a,b,$4,$5,$6}' gencode.v19.transcript.bed) > gencode.cpgIslandExt.txt

### for ENCODE ChIP-seq
ls ../hg19.released.files/*peaks.* | while read files
do
	exp=`echo $files | awk -F"/" '{print $3}' | sed 's/.bed.gz//'`
	bedtools intersect -wo -a <(gunzip -c $files | awk -vOFS="\t" -v name=$exp '{print $1,$2,$3,name}') -b <(awk -vOFS="\t" '{a=$2-300;b=$2+100;if($6=="-"){a=$3-100;b=$3+300}if(a<0){a=0;}print $1,a,b,$4,$5,$6}' gencode.v19.transcript.bed) > gencode.peaks.$exp &
	sleep 0.1s
done

# cat gencode.peaks.* | cut -f 4,8 | sort | uniq -c | awk -vOFS="\t" '{print $2,$3,$1}' | perl /home1/04935/shaojf/stampede2/myTools/UTHealthMut2Loop/relatedScripts/make.occurence_count.matrix.from.GTRD.pl /dev/stdin > ../formated.files/encode.peaks.tsv
# perl ../make.peaks.matrix.pl <(cat gencode.peaks.*) > ../formated.files/encode.peaks.tsv 

for f in `ls gencode.peaks.*`
do
	# cut -f 4,8 $f | sort | uniq -c | awk -vOFS="\t" '{print $2,$3,$1}' >> encode.peaks.txt
	cut -f 4,8 $f | sort | uniq -c | awk -vOFS="\t" '{print $2,$3,$1}' > tmp.$f &
	sleep 10s
done
cat tmp.gencode.peaks.* | perl /home1/04935/shaojf/stampede2/myTools/UTHealthMut2Loop/relatedScripts/make.occurence_count.matrix.from.GTRD.pl /dev/stdin > ../formated.files/encode.peaks.tsv


### ChIA-PET
# awk -F"\t" -v OFS="\t" '$2!="bam" && $2!~"bigBed" && $2!="fastq" && $2!="bigWig" && $39=="hg19" && $42=="released"{print $38,$7"_"$13"_"$3}' metadata.tsv | awk -F"/" '{print $7}' | sed 's/ /./g' | awk '{print $0"."$1}' |  while read line
# do
# 	old=`echo $line | awk '{print $1}'`
# 	new=`echo $line | awk '{print $2}'`
# 	if [ -e $old ]
# 	then
# 		mv $old $new
# 	fi
# done
ls *.bed.gz | while read files
do
	exp=`echo $files | sed 's/.bed.gz//'`
	bedtools intersect -wo -a <(gunzip -c $files | awk -vOFS="\t" -v name=$exp '{print $1,$2,$3,name}') -b <(awk -vOFS="\t" '{a=$2-300;b=$2+100;if($6=="-"){a=$3-100;b=$3+300}if(a<0){a=0;}print $1,a,b,$4,$5,$6}' gencode.v19.transcript.bed) > gencode.loop.$exp &
	sleep 0.1s
done
perl ../make.peaks.matrix.pl <(cat gencode.loop.*) > ../formated.files/encode.loop.tsv 



##### finalize the matrix
perl $myperl encode.transcript.tsv hkg.tsg.srtbyPCA.transcript.class 0 2 | cut -f 4- > hkg.tsg.srtbyPCA.transcript.tpm
perl $myperl encode.peaks.tsv hkg.tsg.srtbyPCA.transcript.class 0 2 | cut -f 4- > hkg.tsg.srtbyPCA.transcript.peaks
perl $myperl encode.loop.tsv hkg.tsg.srtbyPCA.transcript.class 0 2 | cut -f 4- > hkg.tsg.srtbyPCA.transcript.loop

