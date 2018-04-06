#!/bin/sh
myperl=/home1/04935/shaojf/myScripts/add_any_2files_together.pl
list=hkg.tsg.srtbyPCA.class
perl $myperl hkg.tsg.srtbyPCA.class hkg.tsg.srtbyPCA.class.tss.neighbors.sim.annotation 0 3 | cut -f 1-8,10 > a
mv a hkg.tsg.srtbyPCA.class.tss.neighbors.sim.annotation

# ENSG00000257017.4|HP	Liver	+	ENSG00000261701.2|HPR	+	31	protein_coding	protein_coding	Liver
# more gencode.v19.gene.bed | grep -e "ENSG00000257017.4|HP" -e "ENSG00000261701.2|HPR"
### gene class pair; pc vs nc?
##feture.bin_vis.R

### gene expression pair; high vs low? (1k)
head -1 v1.4.log2tpm.median.tsv | cut -f 2- > targets.left.tsv
head -1 v1.4.log2tpm.median.tsv | cut -f 2- > targets.right.tsv
perl $myperl v1.4.log2tpm.median.tsv <(awk -F"\t" '$6 > -1000 && $6 < 1000{print $1}' hkg.tsg.srtbyPCA.class.tss.neighbors.sim.annotation) 0 0 | cut -f 3- >> targets.left.tsv
perl $myperl v1.4.log2tpm.median.tsv <(awk -F"\t" '$6 > -1000 && $6 < 1000{print $4}' hkg.tsg.srtbyPCA.class.tss.neighbors.sim.annotation) 0 0 | cut -f 3- >> targets.right.tsv
paste <(awk -F"\t" -vOFS="\t" 'BEGIN{print "Gene1","Gene2","Type1","Type2","Anno1","Anno2"}{if($6 > -1000 && $6 < 1000){print $1,$4,$2,$9,$7,$8}}' hkg.tsg.srtbyPCA.class.tss.neighbors.sim.annotation) targets.left.tsv targets.right.tsv > targets.tsv

head -1 v1.4.log2tpm.median.tsv | cut -f 2- > targets.div.left.tsv
head -1 v1.4.log2tpm.median.tsv | cut -f 2- > targets.div.right.tsv
perl $myperl v1.4.log2tpm.median.tsv <(awk -F"\t" '$6 > -1000 && $6 < 1000 && $3!=$5 {print $1}' hkg.tsg.srtbyPCA.class.tss.neighbors.sim.annotation) 0 0 | cut -f 3- >> targets.div.left.tsv
perl $myperl v1.4.log2tpm.median.tsv <(awk -F"\t" '$6 > -1000 && $6 < 1000 && $3!=$5 {print $4}' hkg.tsg.srtbyPCA.class.tss.neighbors.sim.annotation) 0 0 | cut -f 3- >> targets.div.right.tsv
paste <(awk -F"\t" -vOFS="\t" 'BEGIN{print "Gene1","Gene2","Type1","Type2","Anno1","Anno2"}{if($6 > -1000 && $6 < 1000 && $3!=$5){print $1,$4,$2,$9,$7,$8}}' hkg.tsg.srtbyPCA.class.tss.neighbors.sim.annotation) targets.div.left.tsv targets.div.right.tsv > targets.div.tsv


head -1 v1.4.log2tpm.median.tsv | cut -f 2- > close.targets.div.left.tsv
head -1 v1.4.log2tpm.median.tsv | cut -f 2- > close.targets.div.right.tsv
perl $myperl v1.4.log2tpm.median.tsv <(awk -F"\t" '$6 > -100 && $6 < 100 && $3!=$5 {print $1}' hkg.tsg.srtbyPCA.class.tss.neighbors.sim.annotation) 0 0 | cut -f 3- >> close.targets.div.left.tsv
perl $myperl v1.4.log2tpm.median.tsv <(awk -F"\t" '$6 > -100 && $6 < 100 && $3!=$5 {print $4}' hkg.tsg.srtbyPCA.class.tss.neighbors.sim.annotation) 0 0 | cut -f 3- >> close.targets.div.right.tsv
paste <(awk -F"\t" -vOFS="\t" 'BEGIN{print "Gene1","Gene2","Type1","Type2","Anno1","Anno2"}{if($6 > -100 && $6 < 100 && $3!=$5){print $1,$4,$2,$9,$7,$8}}' hkg.tsg.srtbyPCA.class.tss.neighbors.sim.annotation) close.targets.div.left.tsv close.targets.div.right.tsv > close.targets.div.tsv


### gene TSS counts; more vs less?
# nearest and furthest distance?
gunzip -c gencode.v19.annotation.gtf.gz | awk -F"\t" '$3=="transcript"{print $1,$4,$5,$7,$9}' | awk -vOFS="\t" '{print $1,$2,$3,$6"|"$14,".",$4,$16"|"$20}' | sed 's/"//g;s/;//g' > gencode.v19.all.transcript.txt

#////////////appris_principal/////////////#
gunzip -c gencode.v19.annotation.gtf.gz | grep -v "level 3" | grep "appris_principal" | awk -F"\t" '$3=="transcript"{print $1,$4,$5,$7,$9}' | awk -vOFS="\t" '{print $1,$2,$3,$6"|"$14,".",$4,$16"|"$20}' | sed 's/"//g;s/;//g' > gencode.v19.all.transcript.appris_principal.leve12.txt


awk -vOFS="\t" '{a=$2-1;b=$2;if($6=="-"){a=$3-1;b=$3}if(a<0){a=0;}print $1,a,b,$4,$5,$6}' gencode.v19.all.transcript.txt | bedtools sort -i - | uniq > gencode.v19.all.transcript.tss.bed
awk -vOFS="\t" '{a=$2-1;b=$2;if($6=="-"){a=$3-1;b=$3}if(a<0){a=0;}print $1,a,b,$4,$5,$6}' gencode.v19.all.transcript.appris_principal.leve12.txt | bedtools sort -i - | uniq > gencode.v19.all.transcript.appris_principal.leve12.tss.bed

cut -f 4 gencode.v19.all.transcript.txt | sort | uniq -c | awk '{print $2"\t"$1}' | sort -k2,2nr > gencode.v19.all.transcript.count
cut -f 4 gencode.v19.all.transcript.tss.bed | sort | uniq -c | awk '{print $2"\t"$1}' | sort -k2,2nr > gencode.v19.all.transcript.tss.count
cut -f 4 gencode.v19.all.transcript.appris_principal.leve12.txt | sort | uniq -c | awk '{print $2"\t"$1}' | sort -k2,2nr > gencode.v19.all.transcript.appris_principal.leve12.count
cut -f 4 gencode.v19.all.transcript.appris_principal.leve12.tss.bed | sort | uniq -c | awk '{print $2"\t"$1}' | sort -k2,2nr > gencode.v19.all.transcript.appris_principal.leve12.tss.count



##### with appris_principal.leve12.tss, to find neibors
while read line
do
	echo $line | tr ' ' "\t" > a
	gene=`cut -f 4 a`
	grep -C 20 -wf a gencode.v19.all.transcript.appris_principal.leve12.tss.bed | grep -v $gene > b
	bedtools closest -D a -a a -b b >> gencode.appris_principal.leve12.tss.neighbors.txt
done < gencode.v19.all.transcript.appris_principal.leve12.tss.bed

perl ~/stampede2/myTools/UTHealthMut2Loop/relatedScripts/get.nearest.pairs.pl <(awk -F"\t" '$3!="/"' hkg.tsg.srtbyPCA.class.appris.tss.neighbors.sim.annotation) > uniq.tss.pairs

perl $myperl <(cut -f 4,6,10,12,13 gencode.appris_principal.leve12.tss.neighbors.txt) hkg.tsg.srtbyPCA.class 0 0 | cut -f 1-2,4- | perl $myperl genecode.gene.all.annotation.txt /dev/stdin 0 0 | cut -f 1-6,8 | perl $myperl genecode.gene.all.annotation.txt /dev/stdin 0 3 | cut -f 1-7,9 | perl $myperl hkg.tsg.srtbyPCA.class /dev/stdin 0 0 | cut -f 1-8,10 > hkg.tsg.srtbyPCA.class.appris.tss.neighbors.sim.annotation

head -1 v1.4.log2tpm.median.tsv | cut -f 2- > appris.close.targets.div.left.tsv
head -1 v1.4.log2tpm.median.tsv | cut -f 2- > appris.close.targets.div.right.tsv
perl $myperl v1.4.log2tpm.median.tsv <(awk -F"\t" '$6 > -100 && $6 < 100 && $3!=$5 {print $1}' hkg.tsg.srtbyPCA.class.tss.neighbors.sim.annotation) 0 0 | cut -f 3- >> appris.close.targets.div.left.tsv
perl $myperl v1.4.log2tpm.median.tsv <(awk -F"\t" '$6 > -100 && $6 < 100 && $3!=$5 {print $4}' hkg.tsg.srtbyPCA.class.tss.neighbors.sim.annotation) 0 0 | cut -f 3- >> appris.close.targets.div.right.tsv
paste <(awk -F"\t" -vOFS="\t" 'BEGIN{print "Gene1","Gene2","Type1","Type2","Anno1","Anno2"}{if($6 > -100 && $6 < 100 && $3!=$5){print $1,$4,$2,$9,$7,$8}}' hkg.tsg.srtbyPCA.class.tss.neighbors.sim.annotation) appris.close.targets.div.left.tsv appris.close.targets.div.right.tsv > appris.close.targets.div.tsv




### gene initiation pair(gro-cap or gro-seq on promoter region); high vs high?
