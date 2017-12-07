args <- commandArgs(TRUE)
expr <- args[1]
mutmat <- args[2]
cor <- args[3]
outpre <- args[4]

# expr <- "mut_exp_matrix/COAD-US.bothWGS.exp.tsv"
# mutmat <- "COAD-US.p.mut.tsv"
# cor <- "correlation/correlated.COAD-US.tsv"
# outpre <- "COAD-US"

# some samples will have duplication, so not use header = T
data <- as.matrix(read.table(expr, sep = "\t", row.names = 1, skip = 1))
gene.name <- rownames(data)
# sample.name <- as.matrix(read.table(expr, sep = "\t", row.names = 1, nrow = 1))
m.data <- as.matrix(read.table(mutmat, sep = "\t", row.names = 1, skip = 1))
m.info <- matrix(unlist(strsplit(rownames(m.data), split="%")), ncol = 3, byrow = T)
m.sample.name <- as.matrix(read.table(mutmat, sep = "\t", row.names = 1, nrow = 1))
c.data <- read.table(cor, sep = "\t")
# c.data[c.data[,1]=="APC" & c.data[,2]=="NIN",]
# m.data[rowSums(m.data)>5,]
for(i in 1:nrow(m.data)){
	print(i)
	mutation <- c()
	gene <- c()
	slop <- c()
	pvalue <- c()
	sample <- c()
	zscore <- c()
	mutexp <- c()
	ctrexp <- c()
	cor.mat <- c()
	mutsample <- paste(m.sample.name[m.data[i,]==1], collapse=",")
	mutpos <- paste(m.info[i,1:2], collapse=":")
	mutgene <- unlist(strsplit(m.info[i,2], split="\\|"))[1]
	tadgene <- intersect(unlist(strsplit(m.info[i,3], split=",")), gene.name)
	for(j in 1:length(tadgene)){
		lm.fit <- summary(lm(data[gene.name==tadgene[j],] ~ m.data[i,]))
		wt <- sd(data[gene.name==tadgene[j], m.sample.name!=mutsample])
		if(wt == 0){
			wt <- 1
		}
		z <- (mean(data[gene.name==tadgene[j],m.sample.name==mutsample]) - mean(data[gene.name==tadgene[j], m.sample.name!=mutsample])) / wt
		if((lm.fit$coefficients[2,4] < 0.1 && !is.na(lm.fit$coefficients[2,4])) || abs(z) > 1.5){
			mutation <- c(mutation, mutpos)
			gene <-c(gene, tadgene[j])
			slop <- c(slop, lm.fit$coefficients[2,1])
			pvalue <- c(pvalue,lm.fit$coefficients[2,4])
			sample <- c(sample, mutsample)
			zscore <- c(zscore, z)
			mutexp <- c(mutexp, paste(data[gene.name==tadgene[j], m.sample.name==mutsample], collapse = ","))
			ctrexp <- c(ctrexp, paste(data[gene.name==tadgene[j], m.sample.name!=mutsample], collapse = ","))
			corval <- c.data[(c.data[,1]==tadgene[j] & c.data[,2]==mutgene) | (c.data[,2]==tadgene[j] & c.data[,1]==mutgene), 3:6]
			if(nrow(corval) == 0){
				corval <- rep("/", 4)
			}
			# dimnames(cor.mat) <- c()
			cor.mat <- rbind(cor.mat, unlist(corval))
		}
	}
	if(length(gene) > 1 && is.element(mutgene, gene) && length(zscore > 0) > 0 && length(zscore < 0) > 0){
		mutflag <- rep("", length(gene))
		mutflag[gene==mutgene] <- "MutatedPromoter"
		expfc <- zscore / zscore[gene==mutgene]
		expflag <- rep("", length(gene))
		expflag[expfc < 0] <- "Opposite"
		expflag[expfc > 0] <- ""
		# dimnames(cor.mat) <- c()
		rownames(cor.mat) <- c()
		colnames(cor.mat) <- c("pearson", "pearson.pvalue", "spearman", "spearman.pvalue")
		lm.res <- data.frame(mutation, gene, slop, pvalue, sample, zscore, mutexp, ctrexp, cor.mat, mutflag, expflag)
		write.table(lm.res, file = paste(outpre, "exp2mut.lm", "tsv", sep = "."), sep= "\t", append = T, row.names = F)
	}
}
