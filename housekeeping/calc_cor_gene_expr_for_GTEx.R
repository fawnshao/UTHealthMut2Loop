args <- commandArgs(TRUE)
# args <- c("tmp.1", "tmp.2","out.correlation.tsv")
datax <- read.table(args[1], header = F, na.strings = "/")
datay <- read.table(args[2], header = F, na.strings = "/")
valuex <- log2(as.matrix(datax[,-1] + 1))
valuey <- log2(as.matrix(datay[,-1] + 1))
res <- c()
for(i in 1:nrow(datax)){
	if(!is.na(valuex[i,]) && !is.na(valuey[i,])){
		res[i] <- cor(valuex[i,], valuey[i,])
	}
	else{
		res[i] <- NA
	}
}
out <- data.frame(datax[,1], datay[,1], res)
colnames(out) <- c("gene1", "gene2", "correlation")
write.table(out, file = args[3], sep = "\t", row.names = F, quote = F)
