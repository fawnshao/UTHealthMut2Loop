expr <- read.table("GSM1536837_06_01_15_TCGA_24.tumor_Rsubread_TPM.txt", header = T, row.names = 1, sep = "\t")

a <- read.table("TAD.exp.ctlMAT", row.names=1)
b <- read.table("TAD.exp.mutMAT", row.names=1)
aa <- apply(a, 1, mean)
c <- wilcox.test(a[,i], b[,1])$p.value
d <- b[,1]/aa
dd <- data.frame(d, b[,1], aa)
colnames(dd) <- c("Fold-Change", "Mut-Sample", "Average-Tumor-Sample")
if(nrow(dd[(dd[,1] > 2 | dd[,1] < 0.5) & dd[,3] > 10,]) > 0 && c$p.value < 0.05){
	write.csv(dd, file = paste("TAD_3","TCGA-AD-A5EJ-01","csv", sep = "."))
}