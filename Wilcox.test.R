args <- commandArgs(TRUE)
ctr <- args[1]
mut <- args[2]
tad <- args[3]
sam <- args[4]
fc <- as.numeric(args[5])

a <- read.table(ctr, row.names=1)
b <- read.table(mut, row.names=1)
aa <- apply(a, 1, mean)
c <- wilcox.test(aa, b[,1])$p.value
d <- b[,1]/aa
dd <- data.frame(d, b[,1], aa)
colnames(dd) <- c("Fold-Change", mut, "Average-Tumor-Sample")
if(nrow(dd[(dd[,1] > fc | dd[,1] < 1/fc) & (dd[,3] > 10 | dd[,3] > 10),]) > 0 && c$p.value < 0.05){
	write.csv(dd, file = paste(tad, sam, "csv", sep = "."))
}