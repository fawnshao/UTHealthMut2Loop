#!/bin/sh
# download data from https://xenabrowser.net/datapages/
# https://xenabrowser.net/datapages/?cohort=TCGA%20PanCanAtlas
# dataset: gene expression RNAseq - Batch effects normalized mRNA datacohort
# TCGA PanCanAtlas

# dataset IDEB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena
# downloadhttps://pancanatlas.xenahubs.net/download/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz; Full metadata
# samples11060
# version2016-12-29
# hubhttps://pancanatlas.xenahubs.net
# type of datagene expression RNAseq
# unitlog2(norm_value+1)
# raw datahttps://www.synapse.org/#!Synapse:syn4976369.3
# input data formatROWs (identifiers) x COLUMNs (samples) (i.e. genomicMatrix)
# rclone sync /home1/04935/shaojf/stampede2/housekeeping_genes/both.pc.and.nc.genes/TCGA.expr/ mygoogle:hkg_tsg/both.pc.and.nc.genes/TCGA.expr/
myperl=/home1/04935/shaojf/myTools/BioinformaticsDaily/textProcess/add_any_2files_together.pl
myheatmap=/home1/04935/shaojf/myTools/BioinformaticsDaily/RVisualization/pheatmap.with.two.annotation.R
perl $myperl Survival_SupplementalTable_S1_20171025_xena_sp <(head_line EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena | awk '{print $2}') 0 0 | cut -f 1,3- > EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.info
cut -f 1-3 EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.info > EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.info.sim

paste <(head -1 hkg.tsg.srtbyPCA.class) <(head -1 EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena) > hkg.tsg.srtbyPCA.AdjustPANCAN_IlluminaHiSeq_RNASeqV2.mat
perl $myperl EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena <(sed 's/|/\t/' hkg.tsg.srtbyPCA.class | tail -n +2) 0 1 | sed 's/\t/|/' >> hkg.tsg.srtbyPCA.AdjustPANCAN_IlluminaHiSeq_RNASeqV2.mat

grep "/" hkg.tsg.srtbyPCA.AdjustPANCAN_IlluminaHiSeq_RNASeqV2.mat | cut -f 1-2 > hkg.tsg.srtbyPCA.noTCGA
grep -v "/" hkg.tsg.srtbyPCA.AdjustPANCAN_IlluminaHiSeq_RNASeqV2.mat > hkg.tsg.srtbyPCA.TCGA.mat
# ENSG00000147381.7|MAGEA4	Testis	MAGEA4	8.02	5.24	0.00	7.94	13.22	0.00	0.00	0.00	0.00	0.00	0.50	2.58	0.00	0.00	0.00	0.62	0.00	0.00	3.26	0.00	3.01	2.46	0.00	1.92	7.02	0.76	1.06	0.00	1.37	0.00	0.42	0.00	11.03	8.05	2.16	0.00	0.00	0.00	0.66	0.00	0.00	0.63	4.98	7.21	0.00	0.00	0.00	0.00	0.00	0.00	10.48	0.57	2.50	0.47	2.89	0.77	0.00	0.00	0.00	0.00	0.00	0.00	0.60	0.00	0.00	0.00	0.00	0.00	6.57	0.00	12.69	0.00	0.00	1.18	10.74	0.00	0.00	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	0.30	2.57	0.61	8.47	0.30	0.30	0.30	0.60	0.30	0.30	0.30	3.95
#############################################################################
# library(data.table)
# library(pheatmap)
# args <- c("hkg.tsg.srtbyPCA.TCGA.mat", "EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.info.sim")
# input <- fread(args[1], sep = "\t", header = T, na.strings = "/")
# class <- fread(args[2], sep = "\t", header = T, na.strings = "/")
# scores <- data.matrix(input[,4:11072])
# rownames(scores) <- as.matrix(input[,1])
# colnames(scores) <- 1:11069
# annosC <- class[,3]
# annosR <- input[,2]
# rownames(annosC) <- colnames(scores)
# rownames(annosR) <- rownames(scores)
# annosR[Type == "hkg1" | Type == "hkg2" | Type == "hkg3" | Type == "hkg4"] <- "HKG"
# annosR[Type != "HKG" & Type != "mixTSG"] <- "singleTSG"
# colors <- colorRampPalette(c("blue", "white", "red"))(100)
# scores[!is.na(scores) & scores > 15] <- 15
# png(filename = paste(args[1], "pheatmap.png", sep = "."), width = 1500, height = 1200)
# myplot <- pheatmap(scores, scale = "none", annotation_row = annosR, annotation_col = annosC,
# 	show_rownames = F, show_colnames = F, color = colors, 
# 	cluster_cols = F, cluster_rows = F)
# dev.off()
#############################################################################
# cut -f 1,228,230 hkg.pMUT.gt5.vert.motifs.txt | grep "(" > ../TCGA.expr/NRF.mut.hkg
# cut -f 1,103 hkg.pMUT.gt5.vert.motifs.txt | grep "(" > ../TCGA.expr/ETS.mut.hkg

for gene in `cut -f1 NRF.mut.hkg | cut -d"@" -f 2 | tr ";" "\n" | cut -f1 -d"|" | tail -n +2 | sort | uniq | grep -v "\-"`
do
	grep -w -e NRF1 -e $gene EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena > NRF1.$gene.tsv
	Rscript /home1/04935/shaojf/myTools/BioinformaticsDaily/RVisualization/scatterplot.for2rows.R NRF1.$gene.tsv
done

##### by mutation
# target=NRF1
for target in POLR2A EP300 MYC CTCF RAD21 TP53 ETV1 ETS1 ETS2 KLF9 YY1 SP1 ELK1 ELK4 NRF1
do
	gunzip -c /home1/04935/shaojf/stampede2/refs/xenahubs/cohort.TCGA.PanCanAtlas/mc3.v0.2.8.PUBLIC.xena.gz | grep -w $target | grep "probably_damaging" > $target.mut.sample
	head_line hkg.tsg.srtbyPCA.TCGA.mat | grep -wf <(cut -f 1 $target.mut.sample) | awk '{printf "%s,", $1}' | sed 's/,$//' | awk '{print "1-3,"$0}' | xargs -n 1 -I mycol cut -f mycol hkg.tsg.srtbyPCA.TCGA.mat > hkg.tsg.srtbyPCA.TCGA.$target.mut.sample
	count=`more $target.mut.sample | wc -l`
	head_line hkg.tsg.srtbyPCA.TCGA.mat | tail -n +4 | grep -vwf <(cut -f 1 $target.mut.sample) | shuf | head -n $count | awk '{printf "%s,", $1}' | sed 's/,$//' | awk '{print "1-3,"$0}' | xargs -n 1 -I mycol cut -f mycol hkg.tsg.srtbyPCA.TCGA.mat > hkg.tsg.srtbyPCA.TCGA.$target.nonmut.sample
	paste hkg.tsg.srtbyPCA.TCGA.$target.mut.sample <(cut -f 4- hkg.tsg.srtbyPCA.TCGA.$target.nonmut.sample) > hkg.tsg.srtbyPCA.TCGA.$target.cat.sample
	echo "ID MutStatus" | tr " " "\t" > hkg.tsg.srtbyPCA.TCGA.$target.cat.sample.annotation
	head -1 hkg.tsg.srtbyPCA.TCGA.$target.mut.sample | cut -f 4- | tr "\n" "\t" | sed 's?\t?\tmut\n?g' >> hkg.tsg.srtbyPCA.TCGA.$target.cat.sample.annotation
	head -1 hkg.tsg.srtbyPCA.TCGA.$target.nonmut.sample | cut -f 4- | tr "\n" "\t" | sed 's?\t?\tnotmut\n?g' >> hkg.tsg.srtbyPCA.TCGA.$target.cat.sample.annotation
	Rscript $myheatmap hkg.tsg.srtbyPCA.TCGA.$target.cat.sample hkg.tsg.srtbyPCA.TCGA.$target.cat.sample.annotation
done
for target in POLR2A EP300 MYC CTCF RAD21 TP53 ETV1 ETS1 ETS2 KLF9 YY1 SP1 ELK1 ELK4 NRF1
do
	awk -F"\t" '$2=="Type" || $2~/hkg/' hkg.tsg.srtbyPCA.TCGA.$target.cat.sample > hkg.TCGA.$target.cat.sample
	Rscript $myheatmap hkg.TCGA.$target.cat.sample hkg.tsg.srtbyPCA.TCGA.$target.cat.sample.annotation
done

target=KAT5
fcR=/home1/04935/shaojf/myTools/BioinformaticsDaily/gene_expression.mutation_rate.R/fc_for_mutvsnonmut.R
gunzip -c /home1/04935/shaojf/stampede2/refs/xenahubs/cohort.TCGA.PanCanAtlas/mc3.v0.2.8.PUBLIC.xena.gz | grep -w $target > $target.mut.sample
head_line EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena | grep -wf <(cut -f 1 $target.mut.sample) | awk '{printf "%s,", $1}' | sed 's/,$//' | awk '{print "1,"$0}' | xargs -n 1 -I mycol cut -f mycol EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena > TCGA.$target.mut.sample
count=`more $target.mut.sample | wc -l`
head_line EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena | tail -n +2 | grep -vwf <(cut -f 1 $target.mut.sample) | shuf | head -n $count | awk '{printf "%s,", $1}' | sed 's/,$//' | awk '{print "1,"$0}' | xargs -n 1 -I mycol cut -f mycol EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena > TCGA.$target.nonmut.sample
paste <(cut -f 1 TCGA.$target.mut.sample | awk '{print $1"\tType\tType"}') <(cut -f 2- TCGA.$target.mut.sample) <(cut -f 2- TCGA.$target.nonmut.sample) > TCGA.$target.cat.sample
echo "ID MutStatus" | tr " " "\t" > TCGA.$target.cat.sample.annotation
head -1 TCGA.$target.mut.sample | cut -f 2- | tr "\n" "\t" | sed 's?\t?\tmut\n?g' >> TCGA.$target.cat.sample.annotation
head -1 TCGA.$target.nonmut.sample | cut -f 2- | tr "\n" "\t" | sed 's?\t?\tnotmut\n?g' >> hkg.tsg.srtbyPCA.TCGA.$target.cat.sample.annotation
Rscript $fcR TCGA.$target.cat.sample TCGA.$target.cat.sample.annotation &
# Rscript $myheatmap TCGA.$target.cat.sample TCGA.$target.cat.sample.annotation &

##### by expression
for target in POLR2A EP300 MYC CTCF RAD21 TP53 ETV1 ETS1 ETS2 KLF9 YY1 SP1 ELK1 ELK4 NRF1
do
	grep -w $target EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena | head_line | grep -vw NA | sort -k2,2n | grep -vw $target | head -100 > $target.low100.sample
	grep -w $target EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena | head_line | grep -vw NA | sort -k2,2nr | head -100 > $target.top100.sample
	awk '{printf "%s,", $1}' $target.low100.sample | sed 's/,$//' | awk '{print "1-3,"$0}' | xargs -n 1 -I mycol cut -f mycol hkg.tsg.srtbyPCA.TCGA.mat > hkg.tsg.srtbyPCA.TCGA.$target.low100.sample
	awk '{printf "%s,", $1}' $target.top100.sample | sed 's/,$//' | awk '{print "1-3,"$0}' | xargs -n 1 -I mycol cut -f mycol hkg.tsg.srtbyPCA.TCGA.mat > hkg.tsg.srtbyPCA.TCGA.$target.top100.sample
	paste hkg.tsg.srtbyPCA.TCGA.$target.low100.sample <(cut -f 4- hkg.tsg.srtbyPCA.TCGA.$target.top100.sample) > hkg.tsg.srtbyPCA.TCGA.$target.100cat.sample
	echo "ID ExprStatus" | tr " " "\t" > hkg.tsg.srtbyPCA.TCGA.$target.100cat.sample.annotation
	head -1 hkg.tsg.srtbyPCA.TCGA.$target.low100.sample | cut -f 4- | tr "\n" "\t" | sed 's?\t?\tlow\n?g' >> hkg.tsg.srtbyPCA.TCGA.$target.100cat.sample.annotation
	head -1 hkg.tsg.srtbyPCA.TCGA.$target.top100.sample | cut -f 4- | tr "\n" "\t" | sed 's?\t?\ttop\n?g' >> hkg.tsg.srtbyPCA.TCGA.$target.100cat.sample.annotation
	Rscript $myheatmap hkg.tsg.srtbyPCA.TCGA.$target.100cat.sample hkg.tsg.srtbyPCA.TCGA.$target.100cat.sample.annotation
done

for target in POLR2A EP300 MYC CTCF RAD21 TP53 ETV1 ETS1 ETS2 KLF9 YY1 SP1 ELK1 ELK4 NRF1
do
	grep -w $target EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena | head_line | grep -vw NA | sort -k2,2n | grep -vw $target | head -30 > $target.low30.sample
	grep -w $target EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena | head_line | grep -vw NA | sort -k2,2nr | head -30 > $target.top30.sample
	awk '{printf "%s,", $1}' $target.low30.sample | sed 's/,$//' | awk '{print "1-3,"$0}' | xargs -n 1 -I mycol cut -f mycol hkg.tsg.srtbyPCA.TCGA.mat > hkg.tsg.srtbyPCA.TCGA.$target.low30.sample
	awk '{printf "%s,", $1}' $target.top30.sample | sed 's/,$//' | awk '{print "1-3,"$0}' | xargs -n 1 -I mycol cut -f mycol hkg.tsg.srtbyPCA.TCGA.mat > hkg.tsg.srtbyPCA.TCGA.$target.top30.sample
	paste hkg.tsg.srtbyPCA.TCGA.$target.low30.sample <(cut -f 4- hkg.tsg.srtbyPCA.TCGA.$target.top30.sample) > hkg.tsg.srtbyPCA.TCGA.$target.30cat.sample
	echo "ID ExprStatus" | tr " " "\t" > hkg.tsg.srtbyPCA.TCGA.$target.30cat.sample.annotation
	head -1 hkg.tsg.srtbyPCA.TCGA.$target.low30.sample | cut -f 4- | tr "\n" "\t" | sed 's?\t?\tlow\n?g' >> hkg.tsg.srtbyPCA.TCGA.$target.30cat.sample.annotation
	head -1 hkg.tsg.srtbyPCA.TCGA.$target.top30.sample | cut -f 4- | tr "\n" "\t" | sed 's?\t?\ttop\n?g' >> hkg.tsg.srtbyPCA.TCGA.$target.30cat.sample.annotation
	Rscript $myheatmap hkg.tsg.srtbyPCA.TCGA.$target.30cat.sample hkg.tsg.srtbyPCA.TCGA.$target.30cat.sample.annotation
done
# target=ETV1
# awk '{printf "%s,", $1}' $target.top30.sample | sed 's/,$//'
# head -1 EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena | cut -f 2307,1916,7693,2272,3059,7544,1803,2085,8060,2139,7647,7770,1850,7864,1821,7934,7736,7611,7574,7968,8008,7883,7942,8049,7877,8018,7536,7661,7697,7959 | tr "\t" "\n" | grep -wf /dev/stdin EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.info.sim 

awk -F"\t" '$3=="BRCA"{print $1}' EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.info.sim > BRCA.id
head_line EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena | grep -wf BRCA.id 
cut -f 1,5413-6630 EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena > EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.BRCA
post=BRCA
for target in POLR2A EP300 MYC CTCF RAD21 TP53 ETV1 ETS1 ETS2 KLF9 YY1 SP1 ELK1 ELK4 NRF1
do
	grep -w $target EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.$post | head_line | grep -vw NA | sort -k2,2n | grep -vw $target | head -30 > $target.low30.sample.$post
	grep -w $target EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.$post | head_line | grep -vw NA | sort -k2,2nr | head -30 > $target.top30.sample.$post
	awk '{printf "%s,", $1}' $target.low30.sample.$post | sed 's/,$//' | awk '{print "1-3,"$0}' | xargs -n 1 -I mycol cut -f mycol hkg.tsg.srtbyPCA.TCGA.mat > hkg.tsg.srtbyPCA.TCGA.$target.low30.sample.$post
	awk '{printf "%s,", $1}' $target.top30.sample.$post | sed 's/,$//' | awk '{print "1-3,"$0}' | xargs -n 1 -I mycol cut -f mycol hkg.tsg.srtbyPCA.TCGA.mat > hkg.tsg.srtbyPCA.TCGA.$target.top30.sample.$post
	paste hkg.tsg.srtbyPCA.TCGA.$target.low30.sample.$post <(cut -f 4- hkg.tsg.srtbyPCA.TCGA.$target.top30.sample.$post) > hkg.tsg.srtbyPCA.TCGA.$target.30cat.sample.$post
	echo "ID ExprStatus" | tr " " "\t" > hkg.tsg.srtbyPCA.TCGA.$target.30cat.sample.annotation.$post
	head -1 hkg.tsg.srtbyPCA.TCGA.$target.low30.sample.$post | cut -f 4- | tr "\n" "\t" | sed 's?\t?\tlow\n?g' >> hkg.tsg.srtbyPCA.TCGA.$target.30cat.sample.annotation.$post
	head -1 hkg.tsg.srtbyPCA.TCGA.$target.top30.sample.$post | cut -f 4- | tr "\n" "\t" | sed 's?\t?\ttop\n?g' >> hkg.tsg.srtbyPCA.TCGA.$target.30cat.sample.annotation.$post
	Rscript $myheatmap hkg.tsg.srtbyPCA.TCGA.$target.30cat.sample.$post hkg.tsg.srtbyPCA.TCGA.$target.30cat.sample.annotation.$post
done
myheatmap2=/home1/04935/shaojf/myTools/BioinformaticsDaily/RVisualization/pheatmap.with.one.annotation.one.cluster.R
post=BRCA
for target in POLR2A EP300 MYC CTCF RAD21 TP53 ETV1 ETS1 ETS2 KLF9 YY1 SP1 ELK1 ELK4 NRF1
do
	awk -F"\t" '$2=="Type" || $2~/hkg/' hkg.tsg.srtbyPCA.TCGA.$target.30cat.sample.$post > hkg.TCGA.$target.30cat.sample.$post
	Rscript $myheatmap2 hkg.TCGA.$target.30cat.sample.$post hkg.tsg.srtbyPCA.TCGA.$target.30cat.sample.annotation.$post &
done
for target in POLR2A EP300 MYC CTCF RAD21 TP53 ETV1 ETS1 ETS2 KLF9 YY1 SP1 ELK1 ELK4 NRF1
do
	Rscript $myheatmap hkg.TCGA.$target.30cat.sample.$post hkg.tsg.srtbyPCA.TCGA.$target.30cat.sample.annotation.$post
done
for target in POLR2A EP300 MYC CTCF RAD21 TP53 ETV1 ETS1 ETS2 KLF9 YY1 SP1 ELK1 ELK4 NRF1
do
	awk '{printf "%s,", $1}' $target.low30.sample.$post | sed 's/,$//' | awk '{print "1,"$0}' | xargs -n 1 -I mycol cut -f mycol EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.$post > TCGA.$target.low30.sample.$post
	awk '{printf "%s,", $1}' $target.top30.sample.$post | sed 's/,$//' | awk '{print "1,"$0}' | xargs -n 1 -I mycol cut -f mycol EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.$post > TCGA.$target.top30.sample.$post
	paste <(cut -f 1 TCGA.$target.low30.sample.$post | awk '{print $1"\tType\tType"}') <(cut -f 2- TCGA.$target.low30.sample.$post) <(cut -f 2- TCGA.$target.top30.sample.$post) > TCGA.$target.30cat.sample.$post
	echo "ID ExprStatus" | tr " " "\t" > TCGA.$target.30cat.sample.annotation.$post
	head -1 TCGA.$target.low30.sample.$post | cut -f 2- | tr "\n" "\t" | sed 's?\t?\tlow\n?g' >> TCGA.$target.30cat.sample.annotation.$post
	head -1 TCGA.$target.top30.sample.$post | cut -f 2- | tr "\n" "\t" | sed 's?\t?\ttop\n?g' >> TCGA.$target.30cat.sample.annotation.$post
	Rscript $myheatmap2 TCGA.$target.30cat.sample.$post TCGA.$target.30cat.sample.annotation.$post &
done



########## top mutated TFs
echo "gene mutsamples" | tr " " "\t" > mutsamples.txt
gunzip -c mc3.v0.2.8.PUBLIC.nonsilentGene.xena.gz | tail -n +2 | awk '{sum=0;for(i=2;i<=NF;i++){sum+=$i}print $1"\t"sum}' >> mutsamples.txt
# sort -k2,2nr mutsamples.txt 
perl $myperl mutsamples.txt <(sed 's/|/\t/' hkg.tsg.srtbyPCA.class | tail -n +2) 0 1 | sed 's/\t/|/' > hkg.tsg.srtbyPCA.mutsamples.txt
for target in POLR2A EP300 MYC CTCF RAD21 TP53 ETV1 ETS1 ETS2 KLF9 YY1 SP1 ELK1 ELK4 NRF1
do
	grep -w $target mutsamples.txt
done

fcR=/home1/04935/shaojf/myTools/BioinformaticsDaily/gene_expression.mutation_rate.R/fc_for_mutvsnonmut.R
### one column is actually one gene, should be t()
for target in `awk -F"\t" '$4!="/" && $4 > 10{print $2}' dnabingding.from.go0003677.mut`
do
	gunzip -c /home1/04935/shaojf/stampede2/refs/xenahubs/cohort.TCGA.PanCanAtlas/mc3.v0.2.8.PUBLIC.xena.gz | grep -w $target > $target.mut.sample
	# gunzip -c /home1/04935/shaojf/stampede2/refs/xenahubs/cohort.TCGA.PanCanAtlas/mc3.v0.2.8.PUBLIC.xena.gz | grep -w $target | grep "probably_damaging" > $target.mut.sample
	head_line hkg.tsg.srtbyPCA.TCGA.mat | grep -wf <(cut -f 1 $target.mut.sample) | awk '{printf "%s,", $1}' | sed 's/,$//' | awk '{print "1-3,"$0}' | xargs -n 1 -I mycol cut -f mycol hkg.tsg.srtbyPCA.TCGA.mat > hkg.tsg.srtbyPCA.TCGA.$target.mut.sample
	count=`head -1 hkg.tsg.srtbyPCA.TCGA.$target.mut.sample | cut -f 4- | awk '{print NF}'`
	head_line hkg.tsg.srtbyPCA.TCGA.mat | tail -n +4 | grep -vwf <(cut -f 1 $target.mut.sample) | shuf | head -n $count | awk '{printf "%s,", $1}' | sed 's/,$//' | awk '{print "1-3,"$0}' | xargs -n 1 -I mycol cut -f mycol hkg.tsg.srtbyPCA.TCGA.mat > hkg.tsg.srtbyPCA.TCGA.$target.nonmut.sample
	paste hkg.tsg.srtbyPCA.TCGA.$target.mut.sample <(cut -f 4- hkg.tsg.srtbyPCA.TCGA.$target.nonmut.sample) > hkg.tsg.srtbyPCA.TCGA.$target.cat.sample
	# echo "ID MutStatus" | tr " " "\t" > hkg.tsg.srtbyPCA.TCGA.$target.cat.sample.annotation
	# head -1 hkg.tsg.srtbyPCA.TCGA.$target.mut.sample | cut -f 4- | tr "\n" "\t" | sed 's?\t?\tmut\n?g' >> hkg.tsg.srtbyPCA.TCGA.$target.cat.sample.annotation
	# head -1 hkg.tsg.srtbyPCA.TCGA.$target.nonmut.sample | cut -f 4- | tr "\n" "\t" | sed 's?\t?\tnotmut\n?g' >> hkg.tsg.srtbyPCA.TCGA.$target.cat.sample.annotation
	# Rscript $myheatmap hkg.tsg.srtbyPCA.TCGA.$target.cat.sample hkg.tsg.srtbyPCA.TCGA.$target.cat.sample.annotation
	awk -F"\t" '$2=="Type" || $2~/hkg/' hkg.tsg.srtbyPCA.TCGA.$target.cat.sample > hkg.TCGA.$target.cat.sample
	Rscript $fcR hkg.TCGA.$target.cat.sample
done
for f in *.fc1.tsv
do
	samples=`head -1 $f | awk '{print NF}'`
	genes=`tail -n +2 $f | wc -l`
	if [ $samples -gt 50 ] && [ $genes -gt 50 ]
	then
		echo $f $samples $genes | tr " " "\t"
	fi
done

for target in POLR2A EP300 MYC CTCF RAD21 TP53 ETV1 ETS1 ETS2 KLF9 YY1 SP1 ELK1 ELK4 NRF1
do
	awk -F"\t" '$2=="Type" || $2~/hkg/' hkg.tsg.srtbyPCA.TCGA.$target.30cat.sample > hkg.TCGA.$target.30cat.sample
	Rscript $fcR hkg.TCGA.$target.30cat.sample
done
