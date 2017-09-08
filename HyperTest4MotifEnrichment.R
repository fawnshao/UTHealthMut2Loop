args <- commandArgs(TRUE)

x <- read.table(args[1], sep = "\t")
mut <- unique(x[,1:3])
promoter_mut <- unique(x[abs(x[,20]) < 1000,1:3])
promoter_motif_mut <- unique(x[abs(x[,20]) < 1000 & abs(x[,13]) < 5,1:3])


a <- read.table("promoter_motif_mutations.motif.count", sep = "\t", row.names = 1)
b <- read.table("homo_sapiens.GRCh38.motiffeatures.20161111.motif.count", sep = "\t", row.names = 1)

pv <- a
for(i in 1:nrow(a)){pv[i,1] <- phyper(a[i,1]-1, b[rownames(b) == rownames(a)[i],1], 562356-b[rownames(b) == rownames(a)[i],1], 955, lower.tail=FALSE)}
write.table(pv, file = "phyper.txt", quote = F)
