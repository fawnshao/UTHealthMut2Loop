args <- commandArgs(TRUE)
expr <- args[1]
out <- args[2]

# mat_expr_$f.tsv
data <- as.matrix(read.table(expr, sep = "\t", row.names = 1))

res <- data.frame()
for(i in 1:(nrow(data) - 1)){
	for(j in (i + 1):nrow(data)){
		print(paste(i,j,sep=" vs "))
		gene1 <- rownames(data)[i]
		gene2 <- rownames(data)[j]
		cor.x <- cor.test(data[i,],data[j,],method = "pearson")
		cor.y <- cor.test(data[i,],data[j,],method = "spearman")
		res <- rbind.data.frame(res, c(gene1, gene2, 
			cor.x$estimate[[1]], cor.x$p.value, cor.y$estimate[[1]], cor.y$p.value))
	}
}
colnames(res) <- c("gene1", "gene2", 
	"pearson.cor", "pearson.pvalues", "spearman.cor", "spearman.pvalues")
write.table(res, file = out, sep = "\t")
