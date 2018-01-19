#!/bin/sh
myperl=/home1/04935/shaojf/myScripts/add_any_2files_together.pl
mycodes=/home1/04935/shaojf/stampede2/loop_shifting_variations/codes
cancergene=/home1/04935/shaojf/stampede2/refs/Oncogene_TumorSuppressor/Cosmic.CancerGeneCensus.all.gene.anno
tissue=$1
# tissue=Cells_EBV-transformed_lymphocytes

gunzip -c $tissue.v7.signif_variant_gene_pairs.txt.gz > $tissue.v7.signif_variant_gene_pairs.txt
cut -f 1 $tissue.v7.signif_variant_gene_pairs.txt | sort | uniq -c | sort -k1,1nr > $tissue.v7.eVariants.count.txt

########## single targeting promoter eVars #########
awk '$1==1{print $2}' $tissue.v7.eVariants.count.txt | perl $myperl $tissue.v7.signif_variant_gene_pairs.txt /dev/stdin 0 0 | cut -f 2- > $tissue.single-egenes.tsv
awk '$3>-1000 && $3<1000 {print $1}' $tissue.single-egenes.tsv | perl $myperl $tissue.single-egenes.tsv /dev/stdin 0 0 | cut -f 2- > $tissue.single-egenes.inpromoters.tsv
### pay attention: how to filter associated SNPs by linkage disequilibrium? ###
awk '$8 < 0' $tissue.single-egenes.inpromoters.tsv > $tissue.single-egenes.inpromoters.down.tsv
awk '$8 > 0' $tissue.single-egenes.inpromoters.tsv > $tissue.single-egenes.inpromoters.up.tsv


########## multiple targeting promoter eVars #########
awk '$1>1{print $2}' $tissue.v7.eVariants.count.txt | perl $myperl $tissue.v7.signif_variant_gene_pairs.txt /dev/stdin 0 0 | cut -f 2- > $tissue.multi-egenes.tsv
awk '$3>-1000 && $3<1000 {print $1}' $tissue.multi-egenes.tsv | perl $myperl $tissue.multi-egenes.tsv /dev/stdin 0 0 | cut -f 2- > $tissue.multi-egenes.inpromoters.tsv


########## gene name and coordinates of eGene for multiple targeting promoter eVars #########
gunzip -c $tissue.v7.egenes.txt.gz | sed -n '2,$p' | cut -f 1-6 | sort | uniq > $tissue.GTEx.gene
perl $myperl $tissue.GTEx.gene $tissue.multi-egenes.inpromoters.tsv 0 1 > $tissue.multi-egenes.inpromoters.tsv.gene
awk -vOFS="\t" '{print $15,$16,$17,$1">"$2">"$14">"$3">"$8,".",$18}' $tissue.multi-egenes.inpromoters.tsv.gene > $tissue.multi-egenes.inpromoters.eGene.bed

#### find any pair with contradict expression effect, one is promoter variation, the other are opposite to promoter variation

########## loop-shifting eVars #########
sed 's/>/\t/g'  $tissue.multi-egenes.inpromoters.eGene.bed | perl $mycodes/opposite_eQTL.pl /dev/stdin > $tissue.multi-egenes.inpromoters.eGene.flags

awk '$NF=="opposite"{print $4}' $tissue.multi-egenes.inpromoters.eGene.flags | sort | uniq | perl $myperl $tissue.multi-egenes.inpromoters.eGene.flags /dev/stdin 3 0 | awk '$NF!="/" || $(NF-1)=="promoter"' | cut -f 2-  > $tissue.multi-egenes.inpromoters.eGene.flags.loop-shifting.txt
awk -vOFS="\t" '{print $1,$2,$3,$4">"$5">"$6">"$7">"$8,$9,$10}' $tissue.multi-egenes.inpromoters.eGene.flags.loop-shifting.txt > $tissue.multi-egenes.inpromoters.eGene.flags.loop-shifting.bed

awk '$NF=="/"{print $4}' $tissue.multi-egenes.inpromoters.eGene.flags | sort | uniq | perl $myperl $tissue.multi-egenes.inpromoters.eGene.flags /dev/stdin 3 0 | awk '$NF!="/" || $(NF-1)=="promoter"' | cut -f 2-  > $tissue.multi-egenes.inpromoters.eGene.flags.loop-shifting.txt
awk -vOFS="\t" '{print $1,$2,$3,$4">"$5">"$6">"$7">"$8,$9,$10}' $tissue.multi-egenes.inpromoters.eGene.flags.loop-shifting.txt > $tissue.multi-egenes.inpromoters.eGene.flags.loop-shifting.bed

awk '$NF=="contradict"{print $4}' $tissue.multi-egenes.inpromoters.eGene.flags | sort | uniq | perl $myperl $tissue.multi-egenes.inpromoters.eGene.flags /dev/stdin 3 0 | awk '$NF!="/" || $(NF-1)=="promoter"' | cut -f 2-  > $tissue.multi-egenes.inpromoters.eGene.flags.contradict.txt
awk -vOFS="\t" '{print $1,$2,$3,$4">"$5">"$6">"$7">"$8,$9,$10}' $tissue.multi-egenes.inpromoters.eGene.flags.contradict.txt > $tissue.multi-egenes.inpromoters.eGene.flags.contradict.bed

########## loop-shifting eVars in TAD #########
awk '{print "chr"$0}' $tissue.multi-egenes.inpromoters.eGene.flags.loop-shifting.bed | bedtools intersect -wao -a - -b GM12878.TAD.bed > $tissue.multi-egenes.inpromoters.eGene.flags.loop-shifting.GM12878.TAD

awk -vOFS="\t" '{a=$2-1;b=$2;if($10=="-"){a=$3-1;b=$3} print $1,a,b,$5,$9,$10}' $tissue.multi-egenes.inpromoters.eGene.flags.loop-shifting.txt | sort | uniq | bedtools sort -i - | bedtools closest -D a -a - -b GM12878.bait.srt.bed | awk '$(NF-1)!="." && $NF > -1000 && $NF < 1000 {print $4"\t"$(NF-1)}' | perl $myperl /dev/stdin $tissue.multi-egenes.inpromoters.eGene.flags.loop-shifting.txt 0 4 > $tissue.multi-egenes.inpromoters.eGene.flags.loop-shifting.GM12878.interaction

#### motif ####
awk '{print $4}' $tissue.multi-egenes.inpromoters.eGene.flags.loop-shifting.txt | sed 's/_/\t/g' | awk -vOFS="\t" '{print $1,$2-1,$2+length($3)-1,$1"_"$2"_"$3"_"$4"_"$5}' | sort | uniq > $tissue.LoopBroken.bed
sh $mycodes/motif_gainorloss.sh $tissue.LoopBroken
perl $mycodes/motif_compare.pl $tissue.LoopBroken
##### backgroud motif #####
awk '{print $2}' $tissue.v7.eVariants.count.txt | sed 's/_/\t/g' | awk -vOFS="\t" '{print $1,$2-1,$2+length($3)-1,$1"_"$2"_"$3"_"$4"_"$5}' | sort | uniq | grep -v variant_id > $tissue.eVar.bed
sh $mycodes/motif_gainorloss.sh $tissue.eVar
perl $mycodes/motif_compare.pl $tissue.eVar

########## final #########
perl $myperl $tissue.LoopBroken.motif.cmp $tissue.multi-egenes.inpromoters.eGene.flags.loop-shifting.GM12878.interaction 0 3 > $tissue.LoopBroken.GM12878.interaction.motif

####### simple results ########
perl $myperl $tissue.LoopBroken.motif.cmp $tissue.multi-egenes.inpromoters.eGene.flags.loop-shifting.txt 0 3 > $tissue.LoopBroken.motif
perl $myperl $cancergene $tissue.LoopBroken.motif 0 5 > $tissue.LoopBroken.motif.cancergene
awk '$8<0 && $11=="promoter" && $12=="/"' $tissue.LoopBroken.motif > $tissue.LoopBroken.promoter.down.motif
awk '$8>0 && $11=="promoter" && $12=="/"' $tissue.LoopBroken.motif > $tissue.LoopBroken.promoter.up.motif
perl $mycodes/calc_neighboring_TSS_dis.pl $tissue.LoopBroken.motif > $tissue.LoopBroken.motif.tssdis
perl $myperl $cancergene $tissue.LoopBroken.motif.tssdis 0 6 > $tissue.LoopBroken.motif.tssdis.cancergene
