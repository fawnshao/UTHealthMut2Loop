args <- commandArgs(TRUE)
expr <- args[1]
outpre <- args[2]

# can not as.matrix, otherwise the expression will be strings, not numeric
data <- read.table(expr, sep = "\t")
# sample.name <- as.vector(data[, 1])
gene.name <- as.vector(data[, 2])
gene.name.uniq <- unique(gene.name)

res <- data.frame()
for(i in 1:(length(gene.name.uniq) - 1)){
	for(j in (i + 1):length(gene.name.uniq)){
		gene1 <- gene.name.uniq[i]
		gene2 <- gene.name.uniq[j]
		all1 <- as.numeric(as.matrix(data[gene.name == gene1, 3]))
		all2 <- as.numeric(as.matrix(data[gene.name == gene2, 3]))
		cor.x <- cor.test(x[,1],x[,2],method = "pearson")
	}
}

colnames(res) <- c("gene1", "gene2", "pearson.cor", "pearson.pvalues", "spearman.cor", "spearman.pvalues")