args <- commandArgs(TRUE)

usrlist <- args[1]
backgroud <- args[2]
outprefix <- args[3]

a <- read.table(usrlist, sep = "\t", row.names = 1)
b <- read.table(backgroud, sep = "\t", row.names = 1)
pv <- a
pv <- data.frame(pv, a[,1], a[,1])
colnames(pv) <- c("hyper.pv", "count", "backgroud.count")
usrlistsum <- sum(a[,1])
backgroudsum <- sum(b[,1])

for(i in 1:nrow(a)){
	pv[i,1] <- phyper(a[i,1] - 1, b[rownames(b) == rownames(a)[i],1], 
		backgroudsum - b[rownames(b) == rownames(a)[i],1], usrlistsum, lower.tail = FALSE)
	pv[i, 2:3] <- c(a[i,1], b[rownames(b) == rownames(a)[i],1])
}
write.csv(pv, file = paste(outprefix, "phyper.csv", sep = "."), quote = F)