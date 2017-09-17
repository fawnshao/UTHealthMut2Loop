args <- commandArgs(TRUE)

usrlist <- args[1]
backgroud <- args[2]
outprefix <- args[3]

ax <- read.table(usrlist, sep = "\t")#, row.names = 1)
bx <- read.table(backgroud, sep = "\t")#, row.names = 1)
pv <- ax
a <- as.matrix(ax[,2])
b <- as.matrix(bx[,2])
pv <- data.frame(pv, a[,1], a[,1], a[,1], a[,1])
colnames(pv) <- c("hyper.pv", "count", "backgroud.count", "count.total", "backgroud.count.total")
usrlistsum <- sum(a[,1])
backgroudsum <- sum(b[,1])

for(i in 1:nrow(a)){
	pv[i,1] <- phyper(a[i,1] - 1, b[bx[,1] == ax[i,1],1], 
		backgroudsum - b[bx[,1] == ax[i,1],1], usrlistsum, lower.tail = FALSE)
	pv[i, 3:6] <- c(a[i,1], b[bx[,1] == ax[i,1],1], usrlistsum, backgroudsum)
}
# for(i in 1:nrow(a)){
# 	pv[i,1] <- phyper(a[i,1] - 1, b[rownames(b) == rownames(a)[i],1], 
# 		backgroudsum - b[rownames(b) == rownames(a)[i],1], usrlistsum, lower.tail = FALSE)
# 	pv[i, 2:5] <- c(a[i,1], b[rownames(b) == rownames(a)[i],1], usrlistsum, backgroudsum)
# }
# write.csv(pv, file = paste(outprefix, "phyper.csv", sep = "."), quote = F)
write.table(pv, file = paste(outprefix, "phyper.tsv", sep = "."), quote = TRUE, sep = "\t")
