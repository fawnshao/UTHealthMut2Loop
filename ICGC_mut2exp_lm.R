# expr <- "exp_seq.COAD-US.tsv.WGS.sim"
# list <- "TAD.exp.IamGroot.Rinput"
# outpre <- "TAD.exp.COAD"
args <- commandArgs(TRUE)
expr <- args[1]
# mutmat <- args[2]
list <- args[2]
outpre <- args[3]

# some samples will have duplication, so not use header = T
data <- as.matrix(read.table(expr, sep = "\t", row.names = 1, skip = 1))
gene.name <- rownames(data)
sample.name <- as.matrix(read.table(expr, sep = "\t", row.names = 1, nrow = 1))
# m.data <- as.matrix(read.table(mutmat, sep = "\t", row.names = 1, skip = 1))
# m.gene.name <- rownames(m.data)
# m.sample.name <- as.matrix(read.table(mutmat, sep = "\t", row.names = 1, nrow = 1))

toprocess <- read.table(list, header = T)
ifLoop <- rep(" ", nrow(toprocess))
# expr.cutoff <- summary(data[data>0])[2]
# print(paste("expression levels:", expr.cutoff, sep = "    "))
for(i in 1:nrow(toprocess)){
	tad.id <- unlist(strsplit(as.vector(toprocess[i,1]), ","))
	p.mut.s <- unlist(strsplit(as.vector(toprocess[i,2]), ","))
	tad.mut.s <- unlist(strsplit(as.vector(toprocess[i,3]), ","))
	tad.gene <- unlist(strsplit(as.vector(toprocess[i,4]), ","))
	p.mut <- as.vector(toprocess[i,5])
	tad.gene.withexpr <- intersect(gene.name, tad.gene)
	tad.gene.withexpr <- tad.gene.withexpr[!is.na(tad.gene.withexpr)]
	if(length(tad.gene.withexpr) > 1){
		x1 <- is.element(sample.name, p.mut.s)
		x2 <- is.element(sample.name, tad.mut.s)
		lm.res <- data.frame()
		tad.expr <- data.frame()
		for(j in 1:length(tad.gene.withexpr)){
			lm.fit <- summary(lm(data[gene.name == tad.gene.withexpr[j], x1==TRUE | x2==FALSE] ~ x1[x1==TRUE | x2==FALSE]))
			# summary(lm(data[gene.name == tad.gene.withexpr[j], ] ~ x1))
			lm.res <- rbind.data.frame(lm.res, lm.fit$coefficients[2,c(1,4)])
			mut <- paste(format(data[gene.name == tad.gene.withexpr[j], x1==TRUE], format = "e", digits = 2), collapse=";")
			ctr <- paste(format(data[gene.name == tad.gene.withexpr[j], x1==FALSE & x2==FALSE], format = "e", digits = 2), collapse=";")
			tad.expr <- rbind.data.frame(tad.expr, cbind(mut, ctr))
		}
		rownames(tad.expr) <- tad.gene.withexpr
		colnames(tad.expr) <- c("mut","non-mut")
		rownames(lm.res) <- tad.gene.withexpr
		colnames(lm.res) <- colnames(lm.fit$coefficients)[c(1,4)]
		lm.res <- data.frame(lm.res, tad.expr)
		lm.res <- na.omit(lm.res[lm.res[,2] < 0.1,])
		# some are FPKM, and some are not. so use median value as cut off
		if(sum(lm.res[,1]>0) > 0 && sum(lm.res[,1]<0) > 0 && is.element(p.mut, rownames(lm.res))){
			zs.mut <- lm.res[is.element(rownames(lm.res),p.mut), 1]
			mut.flag <- rep("", nrow(lm.res))
			mut.flag[is.element(rownames(lm.res),p.mut)] <- "MutatedPromoter"
			alt.flag <- (lm.res[,1] / zs.mut < 0)
			alt.flag[alt.flag == TRUE] <- "Opposite" 
			alt.flag[alt.flag == FALSE] <- "" 
			out.data <- data.frame(lm.res, mut.flag, alt.flag)
			write.table(out.data[out.data[,5]=="MutatedPromoter" | out.data[,6]=="Opposite", ], 
				file = paste(outpre, tad.id, p.mut.s, "tsv", sep = "."), 
				sep = "\t")
			print(paste("Shoot:", outpre, tad.id, p.mut.s, sep = "    "))
			# write.table(data.frame(out, mut.flag), 
			# 	file = paste(x, y, "csv", sep = "."), 
			# 	sep = "\t", quote = TRUE, col.names = TRUE, row.names = TRUE)
			ifLoop[i] <- "LoopBroken"
			# flag <- 0
			# zs.mut <- 1
			# mut.flag <- rep("", nrow(lm.res))
			# for(j in 1:nrow(lm.res)){
			# 	# if(length(grep(pattern = paste("^", rownames(lm.res)[j], "\\|", sep = ""), x = p.mut)) > 0 || 
			# 	# 	length(grep(pattern = paste(";", rownames(lm.res)[j], "\\|", sep = ""), x = p.mut)) > 0){
			# 		if(rownames(lm.res)[j] == p.mut){
			# 		flag <- 1
			# 		mut.flag[j] <- "MutatedPromoter"
			# 		zs.mut <- lm.res[j,1]
			# 	}
			# }
			# if(flag == 1){
			# 	alt.flag <- (lm.res[,1] / zs.mut < 0)
			# 	alt.flag[alt.flag == TRUE] <- "Opposite" 
			# 	alt.flag[alt.flag == FALSE] <- "" 
			# 	out.data <- data.frame(lm.res, mut.flag, alt.flag)
			# 	write.table(out.data[out.data[,5]=="MutatedPromoter" | out.data[,6]=="Opposite", ], 
			# 		file = paste(outpre, tad.id, p.mut.s, "tsv", sep = "."), 
			# 		sep = "\t")
			# 	print(paste("Shoot:", outpre, tad.id, p.mut.s, sep = "    "))
			# 	# write.table(data.frame(out, mut.flag), 
			# 	# 	file = paste(x, y, "csv", sep = "."), 
			# 	# 	sep = "\t", quote = TRUE, col.names = TRUE, row.names = TRUE)
			# 	ifLoop[i] <- "LoopBroken"
			# }
		}
	}
}
write.table(data.frame(toprocess, ifLoop), file = paste(outpre, "TAD.labeled", "tsv", sep = "."), sep= "\t")
