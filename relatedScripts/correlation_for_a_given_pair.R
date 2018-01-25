args <- commandArgs(TRUE)
pairs <- read.table(args[1], sep = "\t")
exprs <- as.matrix(read.table(args[2], sep = "\t", row.names = 1))
genes <- rownames(exprs)
res <- matrix(ncol = 8, nrow = nrow(pairs))
for(i in 1:nrow(pairs)){
	x1 <- exprs[genes==pairs[i,1],]
	x2 <- exprs[genes==pairs[i,2],]
	cor.x <- cor.test(x1, x2, method = "pearson")
	cor.y <- cor.test(x1, x2, method = "spearman")
	res[i,] <- c(pairs[i,], 
		cor.x$estimate[[1]], cor.x$p.value, 
		cor.y$estimate[[1]], cor.y$p.value,
		mean(x1), mean(x2))
}
colnames(res) <- c("gene.rev","gene.fwd",
	"pearson.cor", "pearson.pvalues", 
	"spearman.cor", "spearman.pvalues",
	"ave.gene.rev","ave.gene.fwd")
write.table(res, file = paste(args[1], "cor.tsv", sep = "."), sep = "\t", row.names = FALSE, quote = FALSE)