args <- commandArgs(TRUE)
expr <- args[1]
out <- args[2]
# from <- args[3]
# to <- args[4]

library(parallel)
# Calculate the number of cores
no_cores <- detectCores() - 1
# Initiate cluster
cl <- makeCluster(no_cores)

# mat_expr_$f.tsv
data <- as.matrix(read.table(expr, sep = "\t", row.names = 1))
mycor <- function(i, m){
	print(i)
	gene1 <- rownames(m)[i]
	res <- matrix(ncol = 6, nrow = (nrow(m) - i))
	c <- 1
	for(j in (i + 1):nrow(m)){
		gene2 <- rownames(m)[j]
		cor.x <- cor.test(m[i,],m[j,], method = "pearson")
		cor.y <- cor.test(m[i,],m[j,], method = "spearman")
		if(abs(cor.x$estimate[[1]]) > 0.6 || abs(cor.y$estimate[[1]]) > 0.6){
			res[c,] <- c(gene1, gene2, cor.x$estimate[[1]], cor.x$p.value, cor.y$estimate[[1]], cor.y$p.value)
			c <- c + 1
		}
	}
	colnames(res) <- c("gene1", "gene2", "pearson.cor", "pearson.pvalues", "spearman.cor", "spearman.pvalues")
	return(res)
}

clusterExport(cl, "data")
clusterExport(cl, "mycor")
# a <- parLapply(cl, 1:nrow(data), function(x) {mycor(x,data)})
# out.data <- do.call("rbind", a)
# write.table(out.data, file = out, sep = "\t")
for(s in 1:(floor(nrow(data)/10) - 1)){
	start <- (s - 1) * 10 + 1
	end <- s * 10
	a <- parLapply(cl, start:end, function(x) {mycor(x,data)})
	out.data <- do.call("rbind", a)
	write.table(out.data, file = out, sep = "\t", append = TRUE, quote = FALSE)
}
a <- parLapply(cl, (end + 1):(nrow(data)-1), function(x) {mycor(x,data)})
out.data <- do.call("rbind", a)
write.table(out.data, file = out, sep = "\t", append = TRUE, quote = FALSE)
stopCluster(cl)
save.image(file = paste(out,"RData", sep = "."))
