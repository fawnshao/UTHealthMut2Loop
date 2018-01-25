args <- commandArgs(TRUE)
pairs <- as.matrix(read.table(args[1], sep = "\t"))
exprs <- as.matrix(read.table(args[2], sep = "\t", row.names = 1))
genes <- rownames(exprs)
res <- matrix(ncol = 4, nrow = nrow(pairs))
ave <- matrix(ncol = 2, nrow = nrow(pairs))
for(i in 1:nrow(pairs)){
	x1 <- exprs[genes==pairs[i,1],]
	x2 <- exprs[genes==pairs[i,2],]
	if(length(x1) == length(x2) && length(x1) == ncol(exprs)){
		cor.x <- cor.test(x1, x2, method = "pearson")
		cor.y <- cor.test(x1, x2, method = "spearman")
		res[i,] <- c(cor.x$estimate[[1]], cor.x$p.value, 
			cor.y$estimate[[1]], cor.y$p.value)
	}
	ave[i,] <- c(mean(x1), mean(x2))
}
out <- data.frame(pairs, ave, res)
colnames(out) <- c("gene.rev","gene.fwd",
	"ave.gene.rev","ave.gene.fwd",
	"pearson.cor", "pearson.pvalues", 
	"spearman.cor", "spearman.pvalues")
write.table(out, file = paste(args[1], "cor.tsv", sep = "."), sep = "\t", row.names = FALSE, quote = FALSE)