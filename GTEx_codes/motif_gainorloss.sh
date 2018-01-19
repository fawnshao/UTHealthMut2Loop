#!/bin/sh
# input coordinates, get the flanking sequence
# scan motif with the sequences
# find if a mutation will gain to lose motifs
pre=$1
refseq=/home1/04935/shaojf/scratch/bwa-index/hg19.fa
JASPAR=/home1/04935/shaojf/stampede2/myTools/MEME/motif_databases/JASPAR/JASPAR_CORE_2016_vertebrates.meme
HOCOMOCO=/home1/04935/shaojf/stampede2/myTools/MEME/motif_databases/HUMAN/HOCOMOCOv10_HUMAN_mono_meme_format.meme
awk '{print $4}' $pre.bed | awk -F"_" '{print $0"\t"$3"\t"$4}' > $pre.txt

# flanking 20 bp
awk -vOFS="\t" '{print "chr"$1,$2-21,$2}' $pre.bed | bedtools getfasta -fi $refseq -bed - > $pre.left.fa
awk -vOFS="\t" '{print "chr"$1,$3,$3+20}' $pre.bed | bedtools getfasta -fi $refseq -bed - > $pre.right.fa
awk -F"\t" '{print ">"$1"\n"$2}' $pre.txt | paste $pre.left.fa - $pre.right.fa | sed 's/\t//g' > $pre.ref.fa
awk -F"\t" '{print ">"$1"\n"$3}' $pre.txt | paste $pre.left.fa - $pre.right.fa | sed 's/\t//g' | sed '2~2s/-//g' > $pre.mut.fa
sed -n '2~2p' $pre.ref.fa > $pre.seq
awk -F"\t" '{print ">"$1}' $pre.txt | sed "R $pre.seq" > $pre.a
mv $pre.a $pre.ref.fa
sed -n '2~2p' $pre.mut.fa > $pre.seq
awk -F"\t" '{print ">"$1}' $pre.txt | sed "R $pre.seq" > $pre.a
mv $pre.a $pre.mut.fa
rm $pre.seq

# fimo
fimo -oc JASPAR.$pre.ref $JASPAR $pre.ref.fa
fimo -oc JASPAR.$pre.mut $JASPAR $pre.mut.fa
fimo -oc HOCOMOCO.$pre.ref $HOCOMOCO $pre.ref.fa
fimo -oc HOCOMOCO.$pre.mut $HOCOMOCO $pre.mut.fa

sed -i 's/\t/|/' JASPAR.$pre.ref/fimo.txt
sed -i 's/\t/|/' JASPAR.$pre.mut/fimo.txt
sed -i 's/\t//' HOCOMOCO.$pre.ref/fimo.txt
sed -i 's/\t//' HOCOMOCO.$pre.mut/fimo.txt
cat JASPAR.$pre.ref/fimo.txt HOCOMOCO.$pre.ref/fimo.txt > $pre.ref.fimo.txt
cat JASPAR.$pre.mut/fimo.txt HOCOMOCO.$pre.mut/fimo.txt > $pre.mut.fimo.txt
## may need big memory
# cut -f 1-2 $pre.ref.fimo.txt | grep -vf - $pre.mut.fimo.txt > $pre.motif.gain
# cut -f 1-2 $pre.mut.fimo.txt | grep -vf - $pre.ref.fimo.txt > $pre.motif.lose

# rm $pre.mut.fa $pre.ref.fa 
rm $pre.left.fa $pre.right.fa
rm -rf JASPAR.$pre.ref/ JASPAR.$pre.mut/ HOCOMOCO.$pre.ref/ HOCOMOCO.$pre.mut/
