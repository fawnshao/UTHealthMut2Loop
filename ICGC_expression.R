# expr <- "exp_seq.COAD-US.tsv.WGS.sim"
# list <- "TAD.exp.IamGroot.Rinput"
# outpre <- "TAD.exp.COAD"
args <- commandArgs(TRUE)
expr <- args[1]
list <- args[2]
outpre <- args[3]

data <- read.table(expr, sep = "\t")
sample.name <- as.vector(data[, 1])
gene.name <- as.vector(data[, 2])
# expr.value <- data[,3]

toprocess <- read.table(list, header = T)
ifLoop <- rep(" ", nrow(toprocess))
# TADid promoter_mutated_TCGAsample TCGAsample_with_NCV_in_extended_regions genes_in_the_TAD promoter_mutated_gene
# Z-score (if comparing mt vs. wt) = [(value gene X in mt Y) - (mean gene X in wt)] / (standard deviation of gene X in wt)
# https://www.easycalculation.com/statistics/p-value-for-z-score.php
# p(z=1.645)=0.05 
expr.cutoff <- summary(data[data[,3] > 0,3])[3]
for(i in 1:nrow(toprocess)){
	# print(toprocess[i,1:2])
	tad.id <- unlist(strsplit(as.vector(toprocess[i,1]), ","))
	p.mut.s <- unlist(strsplit(as.vector(toprocess[i,2]), ","))
	tad.mut.s <- unlist(strsplit(as.vector(toprocess[i,3]), ","))
	tad.gene <- unlist(strsplit(as.vector(toprocess[i,4]), ","))
	p.mut <- as.vector(toprocess[i,5])
	tad.gene.withexpr <- intersect(gene.name, tad.gene)
	if(length(tad.gene.withexpr) > 1){
		outlier.flag <- c()
		ctr.mean <- c()
		mut.mean <- c()
		ctr.sd <- c()
		for(j in 1:length(tad.gene.withexpr)){
			ctr <- data[gene.name == tad.gene.withexpr[j] & !is.element(sample.name, tad.mut.s), 3]
			mut <- data[gene.name == tad.gene.withexpr[j] & is.element(sample.name, p.mut.s), 3]
			mut.mean[j] <- mean(mut)
			ctr.mean[j] <- mean(ctr)
			ctr.sd[j] <- sd(ctr)
			outlier.flag[j] <- paste(is.element(mut, boxplot.stats(c(ctr, mut))$out), collapse=",") 
		}
		zscore <- (mut.mean - ctr.mean) / ctr.sd
		fc <- mut.mean / ctr.mean
		t <- data.frame(mut.mean, ctr.mean, ctr.sd, fc, outlier.flag, zscore)
		rownames(t) <- tad.gene.withexpr
		# colnames(t) <- c(p.mut.s, "Average-nonmut-Tumor-Sample", "SD-nonmut-Tumor-Sample", "Fold-Change", "is.outlier", "Z score")
		# summary(data[data[,3] > 0,3])
		# some are FPKM, and some are not. so use median value as cut off
		if(nrow(t[t[,6] > 1.645 & t[,1] > expr.cutoff, ]) > 0 && 
			nrow(t[t[,6] < -1.645 & t[,2] > expr.cutoff, ]) > 0){
			out <- rbind.data.frame(t[t[,6] > 1.645 & t[,1] > expr.cutoff, ], 
				t[t[,6] < -1.645 & t[,2] > expr.cutoff, ])
			flag <- 0
			zs.mut <- 1
			mut.flag <- rep("", nrow(out))
			for(j in 1:nrow(out)){
				if(length(grep(pattern = paste("^", rownames(out)[j], "\\|", sep = ""), x = p.mut)) > 0 || 
					length(grep(pattern = paste(";", rownames(out)[j], "\\|", sep = ""), x = p.mut)) > 0){
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
					file = paste(outpre, tad.id, p.mut.s, "tsv", sep = "."), 
					sep = "\t")
				print(paste("Shoot:", outpre, tad.id, p.mut.s, sep = "    "))
				# write.table(data.frame(out, mut.flag), 
				# 	file = paste(x, y, "csv", sep = "."), 
				# 	sep = "\t", quote = TRUE, col.names = TRUE, row.names = TRUE)
				ifLoop[i] <- "LoopBroken"
			}
		}
	}
}
write.table(data.frame(toprocess, ifLoop), file = paste(outpre, "TAD.labeled", "tsv", sep = "."), sep= "\t")
