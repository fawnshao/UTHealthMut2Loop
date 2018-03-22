library(data.table)
args <- c("GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct", "GTEx_sample.tissue.txt", "GTEx.gene.class.byGTExorder")
# args <- c("GTEx.log2tpm.median.tsv", "GTEx_sample.tissue.txt", "GTEx.gene.class.byGTExorder")
# cut -f 1 v1.5.HKG.tsv | awk '{print $0"\thkg"}' > GTEx.gene.class
# cut -f 1,156 v1.5.TSG.tsv | tail -n +2 > tsg.tsv
# grep -v "," tsg.tsv | sort -k2 >> GTEx.gene.class
# grep "," tsg.tsv | awk '{print $1"\tmixTSG"}' >> GTEx.gene.class
# perl ~/myTools/UTHealthMut2Loop/relatedScripts/add_any_2files_together.pl <(cut -f 1 GTEx.gene.class) <(cut -f 1 GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct | tail -n +2) 0 0 > rest
# awk '$2=="/"{print $1"\tother"}' rest >> GTEx.gene.class 
# perl ~/myTools/UTHealthMut2Loop/relatedScripts/add_any_2files_together.pl GTEx.gene.class <(cut -f 1 GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct)  0 0 | cut -f 1,3 > GTEx.gene.class.byGTExorder

## read in data
# tpm <- as.matrix(read.table(args[1], sep = "\t", row.names = 1, header = T))
print("Reading Data")
raw.table <- fread(args[1], sep = "\t", header = T)
# Read 56202 rows and 11689 (of 11689) columns from 2.622 GB file in 00:06:16
tpm <- log2(data.matrix(raw.table[,-1]) + 1)
rownames(tpm) <- as.matrix(raw.table[,1])
tissues <- fread(args[2], sep = "\t", header = T)
tissues.lab <- tissues[,2]
# labels <- as.matrix(read.table(args[3], sep = "\t"))
class <- fread(args[3], sep = "\t", header = T)
class.lab <- class[,2]
# save.image("raw.data.RData")

## t-SNE
## conda.R
# .libPaths("/work/04935/shaojf/stampede2/myTools/anaconda2/lib/R/library")
library(Rtsne)
# load("raw.data.RData")
set.seed(123)
tsne <- Rtsne(tpm, check_duplicates = FALSE, dims = 2, perplexity = 30, verbose = TRUE, max_iter = 500)
tsne.scale <- Rtsne(tpm, check_duplicates = FALSE, dims = 2, perplexity = 30, verbose = TRUE, max_iter = 500, pca_center = T, pca_scale = T)
t.tsne <- Rtsne(t(tpm), check_duplicates = FALSE, dims = 2, perplexity = 30, verbose = TRUE, max_iter = 500)
t.tsne.scale <- Rtsne(t(tpm), check_duplicates = FALSE, dims = 2, perplexity = 30, verbose = TRUE, max_iter = 500, pca_center = T, pca_scale = T)
# save.image("median.tsne.RData")
# save.image("tsne.RData")
# load("tsne.RData")

## Plotting
library(ggplot2)
class.lab <- as.list(class.lab)
tissues.lab <- as.list(tissues.lab)
v.tsne <- as.data.frame(tsne$Y)
colnames(v.tsne) <- c("tSNE.1", "tSNE.2")
# png(filename = paste(args[1], "subsets.noscale.tsne.png", sep = "."), width = 1500, height = 1200)
pdf(file = paste(args[1], "tsne.pdf", sep = "."), width = 15, height = 12)
ggplot(data = v.tsne, aes(x = tSNE.1, y = tSNE.2)) + 
	geom_point(aes(colour = class.lab, alpha = 1/10)) +
	ggtitle(args[1])
dev.off()
v.tsne <- as.data.frame(tsne.scale$Y)
colnames(v.tsne) <- c("tSNE.1", "tSNE.2")
# png(filename = paste(args[1], "subsets.noscale.tsne.png", sep = "."), width = 1500, height = 1200)
pdf(file = paste(args[1], "tsne.scale.pdf", sep = "."), width = 15, height = 12)
ggplot(data = v.tsne, aes(x = tSNE.1, y = tSNE.2)) + 
	geom_point(aes(colour = class.lab, alpha = 1/10)) +
	ggtitle(args[1])
dev.off()

v.tsne <- as.data.frame(t.tsne$Y)
colnames(v.tsne) <- c("tSNE.1", "tSNE.2")
# png(filename = paste(args[1], "subsets.noscale.tsne.png", sep = "."), width = 1500, height = 1200)
pdf(file = paste(args[1], "t.tsne.pdf", sep = "."), width = 15, height = 12)
ggplot(data = v.tsne, aes(x = tSNE.1, y = tSNE.2)) + 
	geom_point(aes(colour = tissues.lab, alpha = 1/10)) +
	ggtitle(args[1])
dev.off()
v.tsne <- as.data.frame(t.tsne.scale$Y)
colnames(v.tsne) <- c("tSNE.1", "tSNE.2")
# png(filename = paste(args[1], "subsets.noscale.tsne.png", sep = "."), width = 1500, height = 1200)
pdf(file = paste(args[1], "t.tsne.scale.pdf", sep = "."), width = 15, height = 12)
ggplot(data = v.tsne, aes(x = tSNE.1, y = tSNE.2)) + 
	geom_point(aes(colour = tissues.lab, alpha = 1/10)) +
	ggtitle(args[1])
dev.off()
