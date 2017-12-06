args <- commandArgs(TRUE)
expr <- args[1]
mutmat <- args[2]
outpre <- args[3]

# some samples will have duplication, so not use header = T
data <- as.matrix(read.table(expr, sep = "\t", row.names = 1, skip = 1))
gene.name <- rownames(data)
# sample.name <- as.matrix(read.table(expr, sep = "\t", row.names = 1, nrow = 1))
m.data <- as.matrix(read.table(mutmat, sep = "\t", row.names = 1, skip = 1))
m.gene.name <- rownames(m.data)
m.sample.name <- as.matrix(read.table(mutmat, sep = "\t", row.names = 1, nrow = 1))
# m.data[rowSums(m.data)>5,]
for(i in 1:nrow(m.data)){
	mutation <- c()
	gene <- c()
	slop <- c()
	pvalue <- c()
	sample <- c()
	mutsample <- paste(m.sample.name[m.data[i,]==1], collapse=",")
	for(j in 1:nrow(data)){
		lm.fit <- summary(lm(data[j,] ~ m.data[i,]))
		if(lm.fit$coefficients[2,4] < 0.1 && !is.na(lm.fit$coefficients[2,4])){
			mutation <- c(mutation, m.gene.name[i])
			gene <-c(gene, gene.name[j])
			slop <- lm.fit$coefficients[2,1]
			pvalue <- lm.fit$coefficients[2,4]
			sample <- c(sample, mutsample)
		}
	}
	lm.res <- data.frame(mutation, gene, slop, pvalue, sample)
	write.table(lm.res, file = paste(outpre, "exp2mut.lm", "tsv", sep = "."), sep= "\t", append = T, row.names = F)
}
