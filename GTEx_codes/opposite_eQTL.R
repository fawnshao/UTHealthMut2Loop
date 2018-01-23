args <- commandArgs(TRUE)
evarfile <- paste(args[1], "eVar.counts.tsv", sep = ".")
motiffile <- paste(args[1], "motif.cmp.txt", sep = ".")

# $tissue.motif.cmp.txt
# $tissue.eVar.counts.tsv
# do hypergeometric exact test for motif count
# Usage:

#      dhyper(x, m, n, k, log = FALSE)
#      phyper(q, m, n, k, lower.tail = TRUE, log.p = FALSE)
#      qhyper(p, m, n, k, lower.tail = TRUE, log.p = FALSE)
#      rhyper(nn, m, n, k)
     
# Arguments:

#     x, q: vector of quantiles representing the number of white balls
#           drawn without replacement from an urn which contains both
#           black and white balls.

#        m: the number of white balls in the urn.

#        n: the number of black balls in the urn.

#        k: the number of balls drawn from the urn.
evarcount <- read.table(evarfile, sep = "\t", row.names = 1)
usedevarcount <- evarcount[c(9,9,10,10,11,11,12,12,5,5,6,6),]
motifcount <- read.table(motiffile, sep = "\t", header = T, row.names = 1, na.strings = "/")
motifpvalue <- matrix(nrow = nrow(motifcount), ncol = 8)
motiffisher <- matrix(nrow = nrow(motifcount), ncol = 8)
for(i in 1:nrow(motifcount)){
	for(j in 1:4){
		motifpvalue[i, j] <- phyper(motifcount[i, j], motifcount[i, j+4], 
			usedevarcount[j+4] - motifcount[i, j+4], usedevarcount[j], 
			lower.tail = FALSE)
		motifpvalue[i, j+4] <- phyper(motifcount[i, j], motifcount[i, j+8], 
			usedevarcount[j+8] - motifcount[i, j+8], usedevarcount[j], 
			lower.tail = FALSE)
		cmat <- matrix(c(motifcount[i, j],usedevarcount[j]-motifcount[i, j],
			motifcount[i, j+4],usedevarcount[j+4] - motifcount[i, j+4]), nrow = 2, byrow = TRUE)
		if(sum(is.na(cmat)) == 0){
			motiffisher[i, j] <- fisher.test(cmat)$p.value
		}
		cmat <- matrix(c(motifcount[i, j],usedevarcount[j]-motifcount[i, j],
			motifcount[i, j+8],usedevarcount[j+8] - motifcount[i, j+8]), nrow = 2, byrow = TRUE)
		if(sum(is.na(cmat)) == 0){
			motiffisher[i, j+4] <- fisher.test(cmat)$p.value
		}
	}
}
colnames(motifpvalue) <- colnames(motifcount)[5:12]
colnames(motiffisher) <- colnames(motifcount)[5:12]
rownames(motifpvalue) <- rownames(motifcount)
rownames(motiffisher) <- rownames(motifcount)
write.table(rbind(usedevarcount,motifcount), file = paste("final.count", args[1], "tsv", sep = "."), quote = F, sep = "\t")
write.table(motifpvalue, file = paste("final.phyper", args[1], "tsv", sep = "."), quote = F, sep = "\t")
write.table(motiffisher, file = paste("final.fisher", args[1], "tsv", sep = "."), quote = F, sep = "\t")

