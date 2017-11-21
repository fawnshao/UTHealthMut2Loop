args <- commandArgs(TRUE)
expr <- args[1]
out <- args[2]
from <- args[3]
# to <- args[4]

library(parallel)
# Calculate the number of cores
no_cores <- detectCores() - 1
# Initiate cluster
cl <- makeCluster(no_cores)

# mat_expr_$f.tsv
data <- as.matrix(read.table(expr, sep = "\t", row.names = 1))
data <- data[rowSums(data) > 0, ]
# if rowsums==0, spearman correlation will be stopped, and errors caused.
mycor <- function(i, m){
	print(i)
	gene1 <- rownames(m)[i]
	res <- matrix(ncol = 6, nrow = (nrow(m) - i))
	c <- 0
	for(j in (i + 1):nrow(m)){
		gene2 <- rownames(m)[j]
		cor.x <- cor.test(m[i,],m[j,], method = "pearson")
		cor.y <- cor.test(m[i,],m[j,], method = "spearman")
		if(abs(cor.x$estimate[[1]]) > 0.6 || abs(cor.y$estimate[[1]]) > 0.6){
			c <- c + 1
			res[c,] <- c(gene1, gene2, cor.x$estimate[[1]], cor.x$p.value, cor.y$estimate[[1]], cor.y$p.value)
		}
	}
	if(c > 0){
		colnames(res) <- c("gene1", "gene2", "pearson.cor", "pearson.pvalues", "spearman.cor", "spearman.pvalues")
		return(res[1:c,])
	}
	else{
		return()
	}
}

clusterExport(cl, "data")
clusterExport(cl, "mycor")
# a <- parLapply(cl, 1:nrow(data), function(x) {mycor(x,data)})
# out.data <- do.call("rbind", a)
# write.table(out.data, file = out, sep = "\t")
for(s in from:(floor(nrow(data)/10) - 1)){
	start <- (s - 1) * 10 + 1
	end <- s * 10
	print(paste(start, end, sep=" : "))
	a <- parLapply(cl, start:end, function(x) {mycor(x,data)})
	if(!is.null(a)){
		out.data <- do.call("rbind", a)
		write.table(out.data, file = out, sep = "\t", row.names = FALSE, append = TRUE, quote = FALSE)
	}
}
print(paste((end + 1), (nrow(data)-1), sep=" : "))
a <- parLapply(cl, (end + 1):(nrow(data)-1), function(x) {mycor(x,data)})
if(!is.null(a)){
	out.data <- do.call("rbind", a)
	write.table(out.data, file = out, sep = "\t", row.names = FALSE, append = TRUE, quote = FALSE)
}
stopCluster(cl)
save.image(file = paste(out,"RData", sep = "."))
