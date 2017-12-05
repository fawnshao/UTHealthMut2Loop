args <- commandArgs(TRUE)
expr <- args[1]
mutmat <- args[2]
outpre <- args[3]

# some samples will have duplication, so not use header = T
data <- as.matrix(read.table(expr, sep = "\t", row.names = 1, skip = 1))
gene.name <- rownames(data)
sample.name <- as.matrix(read.table(expr, sep = "\t", row.names = 1, nrow = 1))
m.data <- as.matrix(read.table(mutmat, sep = "\t", row.names = 1, skip = 1))
m.gene.name <- rownames(m.data)
m.sample.name <- as.matrix(read.table(mutmat, sep = "\t", row.names = 1, nrow = 1))
# m.data[rowSums(m.data)>5,]
for(i in 1:nrow(m.data)){
	lm.res <- data.frame()
	for(j in 1:nrow(data)){
		lm.fit <- summary(lm(data[j,] ~ m.data[i,]))
		if(lm.fit$coefficients[2,4] < 0.1){
			lm.res <- rbind.data.frame(lm.res, c(m.gene.name[i],gene.name[j],lm.fit$coefficients[2,c(1,4)]))
		}
	}
	write.table(lm.res, file = paste(outpre, "exp2mut.lm", "tsv", sep = "."), sep= "\t", append = T)
}
