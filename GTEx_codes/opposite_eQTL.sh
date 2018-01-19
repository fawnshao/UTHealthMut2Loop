#!/bin/sh
# rclone copy mygoogle:NCmutation/GTEx/opposite_eQTL.pl .
# rclone copy mygoogle:NCmutation/GTEx/motif_compare.pl .
# expdir=/home1/04935/shaojf/stampede2/loop_shifting_variations/GTEx/GTEx_Analysis_v7_eQTL_expression_matrices
myperl=/home1/04935/shaojf/myScripts/add_any_2files_together.pl
myoppo=/home1/04935/shaojf/stampede2/loop_shifting_variations/codes/opposite_eQTL.pl
tissue=Whole_Blood

# https://github.com/intermine/intermine/issues/1607
#      *   File extension: *.egenes.txt.gz
#   Column headers in the file:

#                  gene_id:  GENCODE/Ensembl gene ID
#                gene_name:  GENCODE gene name
#                 gene_chr:  chromosome (gene)
#               gene_start:  gene start position (in base pairs; 1-based coordinates)
#                 gene_end:  gene end position (in base pairs; 1-based coordinates)
#                   strand:  genomic strand
#                  num_var:  number of variants in cis-window
#              beta_shape1:  1st shape parameter of the fitted Beta distribution: B(shape1, shape2)
#              beta_shape2:  2nd shape parameter of the fitted Beta distribution: B(shape1, shape2)
#                  true_df:  Effective degrees of freedom the Beta distribution approximation
#               variant_id:  variant ID in the format {chr}_{pos_first_ref_base}_{ref_seq}_{alt_seq}_b37
#             tss_distance:  Distance between variant and transcription start site. Positive when variant is downstream of the TSS, negative otherwise
#                      chr:  chromosome (variant; same as gene_chr for cis-eQTLs)
#                  snp_pos:  position of the first reference base of the variant
#                      ref:  reference sequence of the variant
#                      alt:  alternate sequence of the variant
# rs_id_dbSNP142_GRCh37p13:  dbSNP142 rsID
#         num_alt_per_site:  number of alternative alleles observed at this site
#     minor_allele_samples:  number of samples carrying the minor allele
#       minor_allele_count:  total number of minor alleles across individuals
#                      maf:  minor allele frequency observed in the set of donors for a given tissue
#               ref_factor:  '1', when the minor allele is the alt base, '-1' when the minor allele is the reference base
#             pval_nominal:  nominal p-value associated with the most significant variant for this gene
#                    slope:  regression slope
#                 slope_se:  standard error of the regression slope
#                pval_perm:  permutation p-value
#                pval_beta:  beta-approximated permutation p-value
#                     qval:  Storey q-value derived from pval_beta
#   pval_nominal_threshold:  nominal p-value threshold for calling a variant-gene pair significant for the gene


# #-----------------------------------------------
# # Significant SNP-gene pairs
# #-----------------------------------------------

#   File extension: *.signif_snpgene_pairs.txt.gz
#   Column headers in the file:

#               variant_id:  variant ID in the format {chr}_{pos_first_ref_base}_{ref_seq}_{alt_seq}_b37
#                  gene_id:  GENCODE/Ensembl gene ID
#             tss_distance:  distance between variant and transcription start site. Positive when variant is downstream of the TSS, negative otherwise
#             pval_nominal:  nominal p-value
#                    slope:  regression slope
#                 slope_se:  standard error of the regression slope
#               slope_fpkm:  regression slope in FPKM units, calculated using quantile normalized expression tables
#            slope_fpkm_se:  standard error of slope_fpkm
#   pval_nominal_threshold:  nominal p-value threshold for calling a variant-gene pair significant for the gene
#         min_pval_nominal:  smallest nominal p-value for the gene
#                pval_beta:  beta-approximated permutation p-value for the gene


# eGene and significant variant-gene associations based on permutations. The archive contains a *.egenes.txt.gz and *.signif_variant_gene_pairs.txt.gz file for each tissue. Note that the *.egenes.txt.gz files contain data for all genes tested; to obtain the list of eGenes, select the rows with 'qval' ≤ 0.05.

################
# gunzip -c Whole_Blood.v7.egenes.txt.gz | awk '$29<=0.05{print $12}' | sort | uniq > egene.sig.evar.txt


# gzip $tissue.v7.signif_variant_gene_pairs.txt.gz
# awk '$3>-1000 && $3<1000' $tissue.v7.signif_variant_gene_pairs.txt > $tissue.v7.promoter.pairs.txt
# cut -f 1 $tissue.v7.promoter.pairs.txt | sort | uniq -c | sort -k1,1nr > $tissue.v7.promoter.eVariants.count.txt

cut -f 1 $tissue.v7.signif_variant_gene_pairs.txt | sort | uniq -c | sort -k1,1nr > $tissue.v7.eVariants.count.txt


########## multiple targeting promoter eVars #########
awk '$1>1{print $2}' $tissue.v7.eVariants.count.txt | perl $myperl $tissue.v7.signif_variant_gene_pairs.txt /dev/stdin 0 0 | cut -f 2- > $tissue.multi-egenes.tsv
awk '$3>-1000 && $3<1000 {print $1}' $tissue.multi-egenes.tsv | perl $myperl $tissue.multi-egenes.tsv /dev/stdin 0 0 | cut -f 2- > $tissue.multi-egenes.inpromoters.tsv


# gunzip -c $expdir/$tissue.v7.normalized_expression.bed.gz | cut -f 1-4 > $tissue.GTEx.gene.tss.bed
gunzip -c $tissue.v7.egenes.txt.gz | sed -n '2,$p' | cut -f 1-6 | sort | uniq > $tissue.GTEx.gene

# perl $myperl $tissue.GTEx.gene.tss.bed $tissue.multi-egenes.inpromoters.tsv 3 1 > $tissue.multi-egenes.inpromoters.tsv.tss
# awk -vOFS="\t" '{print $13,$14,$15,$1">"$2">"$3">"$8}' $tissue.multi-egenes.inpromoters.tsv.tss > $tissue.multi-egenes.inpromoters.eTSS.bed
perl $myperl $tissue.GTEx.gene $tissue.multi-egenes.inpromoters.tsv 0 1 > $tissue.multi-egenes.inpromoters.tsv.gene
awk -vOFS="\t" '{print $15,$16,$17,$1">"$2">"$14">"$3">"$8,".",$18}' $tissue.multi-egenes.inpromoters.tsv.gene > $tissue.multi-egenes.inpromoters.eGene.bed

#### without any restriction, find any pair with contradict expression effect
# awk '$8<0{print $1}' $tissue.multi-egenes.inpromoters.tsv | sort | uniq > $tissue.multi-egenes.inpromoters.neg
# awk '$8>0{print $1}' $tissue.multi-egenes.inpromoters.tsv | sort | uniq > $tissue.multi-egenes.inpromoters.pos

# # cat $tissue.multi-egenes.inpromoters.neg $tissue.multi-egenes.inpromoters.pos | sort | uniq -c | awk '$1>1{print $2}' | wc -l
# cat $tissue.multi-egenes.inpromoters.neg $tissue.multi-egenes.inpromoters.pos | sort | uniq -c | awk '$1>1{print $2}' | perl $myperl $tissue.multi-egenes.inpromoters.tsv.tss /dev/stdin 0 0 | awk -vOFS="\t" '{print $14,$15,$16,$2">"$3">"$4">"$9}' > $tissue.multi-egenes.inpromoters.loop-shifting.tss.bed

#### find any pair with contradict expression effect, one is promoter variation, the other are opposite to promoter variation

########## loop-shifting eVars #########
sed 's/>/\t/g'  $tissue.multi-egenes.inpromoters.eGene.bed | perl $myoppo /dev/stdin > $tissue.multi-egenes.inpromoters.eGene.flags
awk '$NF=="contradict" || $NF=="opposite"{print $4}' $tissue.multi-egenes.inpromoters.eGene.flags | sort | uniq | perl $myperl $tissue.multi-egenes.inpromoters.eGene.flags /dev/stdin 3 0 | awk '$NF!="/" || $(NF-1)=="promoter"' | cut -f 2-  > $tissue.multi-egenes.inpromoters.eGene.flags.loop-shifting.txt

 awk -vOFS="\t" '{print $1,$2,$3,$4">"$5">"$6">"$7">"$8,$9,$10}' $tissue.multi-egenes.inpromoters.eGene.flags.loop-shifting.txt > $tissue.multi-egenes.inpromoters.eGene.flags.loop-shifting.bed

########## loop-shifting eVars in TAD #########
awk '{print "chr"$0}' $tissue.multi-egenes.inpromoters.eGene.flags.loop-shifting.bed | bedtools intersect -wao -a - -b GM12878.TAD.bed > $tissue.multi-egenes.inpromoters.eGene.flags.loop-shifting.GM12878.TAD
awk '{print "chr"$0}' $tissue.multi-egenes.inpromoters.eGene.flags.loop-shifting.bed | bedtools intersect -wao -a - -b CD34.TAD.bed > $tissue.multi-egenes.inpromoters.eGene.flags.loop-shifting.CD34.TAD

# bedtools sort -i $tissue.multi-egenes.inpromoters.eGene.flags.loop-shifting.bed | bedtools closest -D a -a - -b GM12878.bait.srt.bed > $tissue.multi-egenes.inpromoters.eGene.flags.loop-shifting.GM12878.other
# bedtools sort -i $tissue.multi-egenes.inpromoters.eGene.flags.loop-shifting.bed | bedtools closest -D a -a - -b CD34.bait.srt.bed > $tissue.multi-egenes.inpromoters.eGene.flags.loop-shifting.CD34.other

awk -vOFS="\t" '{a=$2-1;b=$2;if($10=="-"){a=$3-1;b=$3} print $1,a,b,$5,$9,$10}' $tissue.multi-egenes.inpromoters.eGene.flags.loop-shifting.txt | sort | uniq | bedtools sort -i - | bedtools closest -D a -a - -b GM12878.bait.srt.bed | awk '$(NF-1)!="." && $NF > -1000 && $NF < 1000 {print $4"\t"$(NF-1)}' | perl $myperl /dev/stdin $tissue.multi-egenes.inpromoters.eGene.flags.loop-shifting.txt 0 4 > $tissue.multi-egenes.inpromoters.eGene.flags.loop-shifting.GM12878.interaction
awk -vOFS="\t" '{a=$2-1;b=$2;if($10=="-"){a=$3-1;b=$3} print $1,a,b,$5,$9,$10}' $tissue.multi-egenes.inpromoters.eGene.flags.loop-shifting.txt | sort | uniq | bedtools sort -i - | bedtools closest -D a -a - -b CD34.bait.srt.bed | awk '$(NF-1)!="." && $NF > -1000 && $NF < 1000 {print $4"\t"$(NF-1)}' | perl $myperl /dev/stdin $tissue.multi-egenes.inpromoters.eGene.flags.loop-shifting.txt 0 4 > $tissue.multi-egenes.inpromoters.eGene.flags.loop-shifting.CD34.interaction

#### motif ####
awk '{print $4}' $tissue.multi-egenes.inpromoters.eGene.flags.loop-shifting.txt | sed 's/_/\t/g' | awk -vOFS="\t" '{print $1,$2-1,$2+length($3)-1,$1"_"$2"_"$3"_"$4"_"$5}' | sort | uniq > $tissue.LoopBroken.bed
sh motif.gainorloss.sh $tissue
perl motif_compare.pl $tissue
#11	117070546	117070550	11_117070547_TGAG_T_b37

########## loop-shifting promoter interaction #########
# awk '$(NF-1)!="." && $NF > -1000 && $NF < 1000' $tissue.multi-egenes.inpromoters.eGene.flags.loop-shifting.GM12878.other | sed 's/>/\t/g' | awk -vFOS="\t" '{print $4,$14}' | more
# awk '$(NF-1)!="." && $NF > -1000 && $NF < 1000' $tissue.multi-egenes.inpromoters.eGene.flags.loop-shifting.GM12878.other | sed 's/>/\t/g' | awk -vFOS="\t" '{print $5,$14}' | perl $myperl /dev/stdin $tissue.multi-egenes.inpromoters.eGene.flags.loop-shifting.txt 0 4 > $tissue.multi-egenes.inpromoters.eGene.flags.loop-shifting.GM12878.interaction
# awk '$(NF-1)!="." && $NF > -1000 && $NF < 1000' $tissue.multi-egenes.inpromoters.eGene.flags.loop-shifting.CD34.other | sed 's/>/\t/g' | awk -vFOS="\t" '{print $5,$14}' | perl $myperl /dev/stdin $tissue.multi-egenes.inpromoters.eGene.flags.loop-shifting.txt 0 4 > $tissue.multi-egenes.inpromoters.eGene.flags.loop-shifting.CD34.interaction


########## final #########
perl $myperl $tissue.LoopBroken.motif.cmp $tissue.multi-egenes.inpromoters.eGene.flags.loop-shifting.GM12878.interaction 0 3 > $tissue.LoopBroken.GM12878.interaction.motif
perl $myperl $tissue.LoopBroken.motif.cmp $tissue.multi-egenes.inpromoters.eGene.flags.loop-shifting.CD34.interaction 0 3 > $tissue.LoopBroken.CD34.interaction.motif


########## stats #########
gunzip -c $tissue.v7.signif_variant_gene_pairs.txt.gz | cut -f 1 | sort | uniq | wc -l
gunzip -c $tissue.v7.signif_variant_gene_pairs.txt.gz | cut -f 2 | sort | uniq | wc -l
gunzip -c $tissue.v7.signif_variant_gene_pairs.txt.gz | wc -l

cut -f 1 $tissue.multi-egenes.tsv | sort | uniq | wc -l
cut -f 2 $tissue.multi-egenes.tsv | sort | uniq | wc -l
wc -l $tissue.multi-egenes.tsv

cut -f 1 $tissue.multi-egenes.inpromoters.tsv | sort | uniq | wc -l
cut -f 2 $tissue.multi-egenes.inpromoters.tsv | sort | uniq | wc -l
wc -l $tissue.multi-egenes.inpromoters.tsv

cut -f 4 $tissue.multi-egenes.inpromoters.eGene.flags.loop-shifting.txt | sort | uniq | wc -l
cut -f 5 $tissue.multi-egenes.inpromoters.eGene.flags.loop-shifting.txt | sort | uniq | wc -l
wc -l $tissue.multi-egenes.inpromoters.eGene.flags.loop-shifting.txt

sed 's/>/\t/g' $tissue.multi-egenes.inpromoters.eGene.flags.loop-shifting.GM12878.TAD | awk -v OFS="\t" '$NF>0{print $4,$5,$11":"$12":"$13}' | cut -f 1 | sort | uniq | wc -l
sed 's/>/\t/g' $tissue.multi-egenes.inpromoters.eGene.flags.loop-shifting.GM12878.TAD | awk -v OFS="\t" '$NF>0{print $4,$5,$11":"$12":"$13}' | cut -f 2 | sort | uniq | wc -l
sed 's/>/\t/g' $tissue.multi-egenes.inpromoters.eGene.flags.loop-shifting.GM12878.TAD | awk -v OFS="\t" '$NF>0{print $4,$5,$11":"$12":"$13}' | cut -f 1,2 | sort | uniq | wc -l
sed 's/>/\t/g' $tissue.multi-egenes.inpromoters.eGene.flags.loop-shifting.GM12878.TAD | awk -v OFS="\t" '$NF>0{print $4,$5,$11":"$12":"$13}' | cut -f 3 | sort | uniq | wc -l

sed 's/>/\t/g' $tissue.multi-egenes.inpromoters.eGene.flags.loop-shifting.CD34.TAD | awk -v OFS="\t" '$NF>0{print $4,$5,$11":"$12":"$13}' | cut -f 1 | sort | uniq | wc -l
sed 's/>/\t/g' $tissue.multi-egenes.inpromoters.eGene.flags.loop-shifting.CD34.TAD | awk -v OFS="\t" '$NF>0{print $4,$5,$11":"$12":"$13}' | cut -f 2 | sort | uniq | wc -l
sed 's/>/\t/g' $tissue.multi-egenes.inpromoters.eGene.flags.loop-shifting.CD34.TAD | awk -v OFS="\t" '$NF>0{print $4,$5,$11":"$12":"$13}' | cut -f 1,2 | sort | uniq | wc -l
sed 's/>/\t/g' $tissue.multi-egenes.inpromoters.eGene.flags.loop-shifting.CD34.TAD | awk -v OFS="\t" '$NF>0{print $4,$5,$11":"$12":"$13}' | cut -f 3 | sort | uniq | wc -l

sed 's/>/\t/g' $tissue.multi-egenes.inpromoters.eGene.flags.loop-shifting.GM12878.TAD | awk -v OFS="\t" '$NF>0{print $4,$5,$11":"$12":"$13}' | cut -f 1,3 | sort | uniq | wc -l
sed 's/>/\t/g' $tissue.multi-egenes.inpromoters.eGene.flags.loop-shifting.CD34.TAD | awk -v OFS="\t" '$NF>0{print $4,$5,$11":"$12":"$13}' | cut -f 1,3 | sort | uniq | wc -l

awk '$NF!="/"{print $4}' $tissue.multi-egenes.inpromoters.eGene.flags.loop-shifting.GM12878.interaction | sort | uniq  |wc -l
awk '$NF!="/"{print $5}' $tissue.multi-egenes.inpromoters.eGene.flags.loop-shifting.GM12878.interaction | sort | uniq  |wc -l
awk '$NF!="/"{print $4,$5}' $tissue.multi-egenes.inpromoters.eGene.flags.loop-shifting.GM12878.interaction | sort | uniq  |wc -l
awk '$NF!="/"{print $NF}' $tissue.multi-egenes.inpromoters.eGene.flags.loop-shifting.GM12878.interaction | sort | uniq  |wc -l

awk '$NF!="/"{print $4}' $tissue.multi-egenes.inpromoters.eGene.flags.loop-shifting.CD34.interaction | sort | uniq  |wc -l
awk '$NF!="/"{print $5}' $tissue.multi-egenes.inpromoters.eGene.flags.loop-shifting.CD34.interaction | sort | uniq  |wc -l
awk '$NF!="/"{print $4,$5}' $tissue.multi-egenes.inpromoters.eGene.flags.loop-shifting.CD34.interaction | sort | uniq  |wc -l
awk '$NF!="/"{print $NF}' $tissue.multi-egenes.inpromoters.eGene.flags.loop-shifting.CD34.interaction | sort | uniq  |wc -l

