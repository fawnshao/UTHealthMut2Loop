library(data.table)
library(ggplot2)
library(pheatmap)
args <- c("housekeepinggene.GTRD.mat", 
           "tissuespecificgene.GTRD.mat",
           "GTEx.GTRD.mat","GTRD.maxcells.GTEx.hkg.tsg.txt")
# args <- c("housekeepinggene.GTRD.count.mat", 
#           "tissuespecificgene.GTRD.count.mat",
#           "GTEx.GTRD.count.mat","GTRD.count.GTEx.hkg.tsg.txt")
# hkg <- read.table(args[1], sep = "\t", header = T, row.names = 1)
hkg <- fread(args[1], sep = "\t", header = T)
tsg <- fread(args[2], sep = "\t", header = T)
all <- fread(args[3], sep = "\t", header = T)
anno <- as.matrix(read.table(args[4], sep = "\t", header = T, row.names = 1))

hkg.score <- hkg[,-1]
tsg.score <- tsg[,-1]
all.score <- all[,-1]
rownames(hkg.score) <- as.matrix(hkg[,1])
rownames(tsg.score) <- as.matrix(tsg[,1])
rownames(all.score) <- as.matrix(all[,1])

hkg.withTF <- apply(hkg.score, 2, function(x) {length(x[x>0])})
tsg.withTF <- apply(tsg.score, 2, function(x) {length(x[x>0])})
all.withTF <- apply(all.score, 2, function(x) {length(x[x>0])})

# pop size : 5260
# sample size : 131
# Number of items in the pop that are classified as successes : 1998
# Number of items in the sample that are classified as successes : 62
# phyper(62-1, 1998, 5260-1998, 131, lower.tail=FALSE)
# pvalues[1] <- phyper(count[1,2], count[1,1], count[2,1]-count[1,1], count[2,2], lower.tail = FALSE, log.p = FALSE)
pvalues.toall <- c()
pvalues.totsg <- c()
for(i in 1:nrow(hkg.withTF)){
	pvalues.toall[i] <- phyper(hkg.withTF[i], 
		lower.tail = FALSE, log.p = FALSE)
}
