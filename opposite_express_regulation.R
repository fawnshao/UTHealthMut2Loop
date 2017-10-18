expr <- "COAD.expr.txt"
list <- "IamGroot.Rinput"
wgsID <- "TAD.exp.mutated.sampleid"
# fc <- 2
args <- commandArgs(TRUE)
expr <- args[1]
list <- args[2]
# fc <- as.numeric(args[3])
wgsID <- args[3]
outpre <- args[4]

data <- read.table(expr, sep = "\t")
gene.name <- as.vector(data[-1, 1])
sample.name <- data[1, -1]
sample.name.sim <- apply(sample.name, 2, function(x){substr(x, start = 1, stop = 15)})
expr.value <- matrix(as.numeric(as.matrix(data[-1, -1])), byrow = F, nrow = length(gene.name))

toprocess <- read.table(list, header = T)
ifLoop <- rep(" ", nrow(toprocess))
id.with.WGS <- read.table(wgsID)
# TADid promoter_mutated_TCGAsample TCGAsample_with_NCV_in_extended_regions genes_in_the_TAD promoter_mutated_gene
# Z-score (if comparing mt vs. wt) = [(value gene X in mt Y) - (mean gene X in wt)] / (standard deviation of gene X in wt)
# https://www.easycalculation.com/statistics/p-value-for-z-score.php
# p(z=1.645)=0.05 
for(i in 1:nrow(toprocess)){
	x <- unlist(strsplit(as.vector(toprocess[i,1]), ","))
	y <- unlist(strsplit(as.vector(toprocess[i,2]), ","))
	z <- unlist(strsplit(as.vector(toprocess[i,3]), ","))
	u <- unlist(strsplit(as.vector(toprocess[i,4]), ","))
	p <- as.vector(toprocess[i,5])
	un <- length(intersect(gene.name, u))
	if(un > 1){
		ctr <- expr.value[is.element(gene.name, u), !is.element(sample.name.sim, z)  & is.element(sample.name.sim, as.character(id.with.WGS[,1]))]
		mut <- matrix(expr.value[is.element(gene.name, u), is.element(sample.name.sim, y)], 
			byrow = F, nrow = un)
		ctr <- log(ctr+1, 2)
		mut <- log(mut+1, 2)
		ctr.mean <- apply(ctr, 1, mean)
		ctr.sd <- apply(ctr, 1, sd)
		mut.mean <- apply(mut, 1, mean)
		# outlier is too harsh
		outlier.flag <- c()
		for(j in 1:nrow(ctr)){
			outlier.flag[j] <- paste(is.element(mut[j,], boxplot.stats(c(ctr[j,], mut[j,]))$out), collapse=",") 
		}
		# z score
		zscore <- c()
		# zscore.flag <- c()
		for(j in 1:un){
			zscore[j] <- (mut.mean[j] - ctr.mean[j]) / ctr.sd[j] # sd(ctr[j,])
			# zscore[j] <- paste(zs, collapse=",")
		}
		# v <- wilcox.test(ctr.mean, mut.mean)$p.value
		w <- mut.mean - ctr.mean
		# t <- data.frame(w, mut.mean, ctr.mean, outlier.flag, zscore)
		t <- data.frame(mut.mean, ctr.mean, ctr.sd, w, outlier.flag, zscore)
		rownames(t) <- gene.name[is.element(gene.name, u)]
		# colnames(t) <- c("Fold-Change", paste(y, v, sep = ": "), "Average-nonmut-Tumor-Sample", "is.outlier")
		# colnames(t) <- c("Fold-Change", paste(y, v, sep = ": "), "Average-nonmut-Tumor-Sample", "is.outlier", "Z score")
		colnames(t) <- c(y, "Average-nonmut-Tumor-Sample", "SD-nonmut-Tumor-Sample", "log2(Fold-Change)", "is.outlier", "Z score")
		# if(nrow(t[t[,1] > fc & t[,2] > 10, ]) > 0 && nrow(t[t[,1] < 1/fc & t[,3] > 10, ]) > 0){
		if(nrow(t[t[,6] > 1.645 & t[,1] > 3, ]) > 0 && nrow(t[t[,6] < -1.645 & t[,2] > 3, ]) > 0){
			# out <- rbind.data.frame(t[t[,1] > fc & t[,1] > 10, ], t[t[,1] < 1/fc & t[,2] > 10, ])
			out <- rbind.data.frame(t[t[,6] > 1.645 & t[,1] > 3, ], t[t[,6] < -1.645 & t[,2] > 3, ])
			flag <- 0
			zs.mut <- 1
			mut.flag <- rep("", nrow(out))
			for(j in 1:nrow(out)){
				if(length(grep(pattern = paste("^", rownames(out)[j], "\\|", sep = ""), x = p) > 0) || 
					length(grep(pattern = paste(";", rownames(out)[j], "\\|", sep = ""), x = p) > 0)){
					flag <- 1
					mut.flag[j] <- "MutatedPromoter"
					zs.mut <- out[j,6]
				}
			}
			alt.flag <- out[,6]/zs.mut < 0
			alt.flag[alt.flag == TRUE] <- "Opposite" 
			alt.flag[alt.flag == FALSE] <- "" 
			if(flag == 1){
				write.table(data.frame(out, mut.flag, alt.flag), 
					file = paste(outpre, x, y, "tsv", sep = "."), 
					sep = "\t")
				# write.table(data.frame(out, mut.flag), 
				# 	file = paste(x, y, "csv", sep = "."), 
				# 	sep = "\t", quote = TRUE, col.names = TRUE, row.names = TRUE)
				ifLoop[i] <- "LoopBroken"
			}
		}
	}
}
write.table(data.frame(toprocess, ifLoop), file = paste(outpre, "TAD.labeled", "tsv", sep = "."), sep= "\t")
