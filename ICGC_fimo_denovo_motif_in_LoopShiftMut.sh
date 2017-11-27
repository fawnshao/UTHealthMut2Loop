#!/bin/sh
# input coordinates, get the flanking sequence
# scan motif with the sequences
# find if a mutation will gain to lose motifs
pre=$1
refseq=/home1/04935/shaojf/scratch/bwa-index/hg19.fa
JASPAR=/home1/04935/shaojf/stampede2/myTools/MEME/motif_databases/JASPAR/JASPAR_CORE_2016_vertebrates.meme
HOCOMOCO=/home1/04935/shaojf/stampede2/myTools/MEME/motif_databases/HUMAN/HOCOMOCOv10_HUMAN_mono_meme_format.meme
coordinates=$pre.LoopBroken.bed
rawmutation=simple_somatic_mutation.open.$pre.tsv.gz
gunzip -c $rawmutation | cut -f 2,9-10,14,16-17 | uniq > $pre.rawmut.txt
sed 's/~/\t/g' $pre.LoopBroken.bed | sort | uniq | awk '{print $5"\t"$1"\t"$3}' | grep -wf - $pre.rawmut.txt > $pre.LoopBroken.mut.txt

# flanking 20 bp
awk -vOFS="\t" '{print "chr"$2,$3-21,$3-1}' $pre.LoopBroken.mut.txt | bedtools getfasta -fi $refseq -bed - > $pre.LoopBroken.left.fa
awk -vOFS="\t" '{print "chr"$2,$3,$3+20}' $pre.LoopBroken.mut.txt | bedtools getfasta -fi $refseq -bed - > $pre.LoopBroken.right.fa
awk -F"\t" '{print ">"$1"."$2"."$3"\n"$5}' $pre.LoopBroken.mut.txt | paste $pre.LoopBroken.left.fa - $pre.LoopBroken.right.fa | sed 's/\t//g' > $pre.ref.fa
awk -F"\t" '{print ">"$1"."$2"."$3"\n"$6}' $pre.LoopBroken.mut.txt | paste $pre.LoopBroken.left.fa - $pre.LoopBroken.right.fa | sed 's/\t//g' | sed '2~2s/-//g' > $pre.mut.fa

# put the gene name in the fasta head
grep ">" $pre.ref.fa  | cut -d">" -f 3 | awk -F"." '{print $2"\t"$3-1}' | grep -wf - $pre.LoopBroken.bed | sed 's/~/\t/g' | awk '{print $1"."$3"\t"$6"~"$7}' | uniq > $pre.a
echo -n "" > $pre.b
grep ">" $pre.ref.fa  | cut -d">" -f 3 | awk -F"." '{print $2"."$3}' | while read line
do
	grep -w $line $pre.a | awk '{print ">"$1"~"$2"\n"}' >> $pre.b
done
paste $pre.b $pre.ref.fa | sed 's/^\t//' | awk '{print $1}' > $pre.c
mv $pre.c $pre.ref.fa
paste $pre.b $pre.mut.fa | sed 's/^\t//' | awk '{print $1}' > $pre.c
mv $pre.c $pre.mut.fa
rm $pre.a $pre.b

# fimo
fimo -oc JASPAR.$pre.ref $JASPAR $pre.ref.fa
fimo -oc JASPAR.$pre.mut $JASPAR $pre.mut.fa
fimo -oc HOCOMOCO.$pre.ref $HOCOMOCO $pre.ref.fa
fimo -oc HOCOMOCO.$pre.mut $HOCOMOCO $pre.mut.fa

sed -i 's/\t/|/' JASPAR.$pre.ref/fimo.txt
sed -i 's/\t/|/' JASPAR.$pre.mut/fimo.txt
cat JASPAR.$pre.ref/fimo.txt HOCOMOCO.$pre.ref/fimo.txt > $pre.ref.fimo.txt
cat JASPAR.$pre.mut/fimo.txt HOCOMOCO.$pre.mut/fimo.txt > $pre.mut.fimo.txt
cut -f 1-2 $pre.ref.fimo.txt | grep -vf - $pre.mut.fimo.txt > $pre.LoopBroken.motif.gain
cut -f 1-2 $pre.mut.fimo.txt | grep -vf - $pre.ref.fimo.txt > $pre.LoopBroken.motif.lose

rm $pre.mut.fa $pre.ref.fa $pre.LoopBroken.left.fa $pre.LoopBroken.right.fa
rm $pre.rawmut.txt
rm -rf JASPAR.$pre.ref/ JASPAR.$pre.mut/ HOCOMOCO.$pre.ref/ HOCOMOCO.$pre.mut/
