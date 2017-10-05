args <- commandArgs(TRUE)
ctr <- args[1]
mut <- args[2]
tad <- args[3]
sam <- args[4]
fc <- as.numeric(args[5])

a <- read.table(ctr, row.names=1)
b <- read.table(mut, row.names=1)
aa <- apply(a, 1, mean)
bb <- apply(b, 1, mean)
c <- wilcox.test(aa, bb)$p.value
d <- bb/aa
dd <- data.frame(d, bb, aa)
colnames(dd) <- c("Fold-Change", sam, "Average-Tumor-Sample")
if(nrow(dd[(dd[,1] > fc | dd[,1] < 1/fc) & (dd[,3] > 10 | dd[,3] > 10),]) > 0 && c < 0.05){
	write.csv(dd, file = paste(tad, sam, "csv", sep = "."))
}