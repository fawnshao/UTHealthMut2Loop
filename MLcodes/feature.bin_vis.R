########
### cage
library(ggplot2)
library(data.table)

args <- c("hkg.tsg.srtbyPCA.class.cage.sim")
input <- fread(args[1], sep = "\t", header = F, na.strings = "/")
types.sim <- as.matrix(input[,2])[,1]
types.sim[grep("hkg",types.sim)] <- "HKG"
types.sim[types.sim!="HKG" & types.sim!="mixTSG"] <- "singleTSG"
# types.sim[grep("hkg",types.sim)] <- "HKG: 2551"
# types.sim[types.sim!="HKG: 2551" & types.sim!="mixTSG"] <- "singleTSG: 2735"
# types.sim[types.sim=="mixTSG"] <- "mixTSG: 1574"
datax <- data.frame(types.sim, input[,3])
colnames(datax) <- c("class", "cage.dis2tss")
sums <- c(2551, 1574, 2735)
names(sums) <- c("HKG", "mixTSG", "singleTSG")
datax <- as.data.table(datax)
datax[class == 'HKG', y := sums[1]]
datax[class == 'mixTSG', y := sums[2]]
datax[class == 'singleTSG', y := sums[3]]

myplot <- ggplot(data = datax, aes(x = cage.dis2tss, colour = class)) + 
		geom_freqpoly(binwidth = 10, position = "identity", data = datax[class == 'HKG',], aes(y = ..count.. / 2551)) +
		geom_freqpoly(binwidth = 10, position = "identity", data = datax[class == 'mixTSG',], aes(y = ..count.. / 1574)) +
		geom_freqpoly(binwidth = 10, position = "identity", data = datax[class == 'singleTSG',], aes(y = ..count.. / 2735)) +
		scale_x_continuous(limits = c(-100,2000)) +
		labs(x = "relative postion of CAGE to TSS", y = "percentage") +
		ggtitle(args[1])
png(filename = paste(args[1], "cage.hist.png", sep = "."), width = 1000, height = 500)
print(myplot)
dev.off()

########
### CpG
library(ggplot2)
library(data.table)

args <- c("hkg.tsg.srtbyPCA.class.cpg.sim.1")
input <- fread(args[1], sep = "\t", header = F, na.strings = "/")
types.sim <- as.matrix(input[,2])[,1]
types.sim[grep("hkg",types.sim)] <- "HKG"
types.sim[types.sim!="HKG" & types.sim!="mixTSG"] <- "singleTSG"

datax <- data.frame(input[,3], types.sim)
colnames(datax) <- c("CpG.pos", "class")
# datax <- as.data.table(datax)

myplot <- ggplot(data = datax, aes(x = CpG.pos, colour = class)) + 
		geom_freqpoly(binwidth = 10, position = "identity", data = datax[types.sim == "HKG",], aes(y = ..count.. / 2551)) +
		geom_freqpoly(binwidth = 10, position = "identity", data = datax[types.sim == "mixTSG",], aes(y = ..count.. / 1574)) +
		geom_freqpoly(binwidth = 10, position = "identity", data = datax[types.sim == "singleTSG",], aes(y = ..count.. / 2735)) +
		scale_x_continuous(limits = c(-999,999)) +
		labs(x = "relative postion of CpG to TSS", y = "percentage") +
		ggtitle(args[1])
png(filename = paste(args[1], "CpG.hist.png", sep = "."), width = 1000, height = 500)
print(myplot)
dev.off()

# x <- as.matrix(input)
# y <- data.frame()
# for(i in 1:nrow(x)){
# 	print(i)
# 	base <- seq(x[i,3], x[i,4])
# 	labs <- rep(x[i,2], length(base))
# 	y <- rbind.data.frame(y, data.frame(base, labs))
# }

########
### TATA and INR
library(ggplot2)
library(data.table)

args <- c("hkg.promotermotif.dis.txt", "mixTSG.promotermotif.dis.txt", 
	"singleTSG.promotermotif.dis.txt")
input1 <- fread(args[1], sep = "\t", header = F, na.strings = "/")
input2 <- fread(args[2], sep = "\t", header = F, na.strings = "/")
input3 <- fread(args[3], sep = "\t", header = F, na.strings = "/")

types.sim <- c(rep("HKG",nrow(input1)), rep("mixTSG",nrow(input2)), rep("singleTSG",nrow(input3)))
datax <- data.frame(rbind.data.frame(input1, input2, input3), types.sim)
colnames(datax) <- c("motif", "dis2tss", "class")

myplot <- ggplot(data = datax, aes(x = dis2tss, colour = class)) + 
		geom_freqpoly(binwidth = 10, position = "identity", data = datax[types.sim == "HKG",], aes(y = ..count.. / 2551)) +
		geom_freqpoly(binwidth = 10, position = "identity", data = datax[types.sim == "mixTSG",], aes(y = ..count.. / 1574)) +
		geom_freqpoly(binwidth = 10, position = "identity", data = datax[types.sim == "singleTSG",], aes(y = ..count.. / 2735)) +
		facet_wrap( ~ motif, nrow = 2) + 
		scale_x_continuous(limits = c(-999,4999)) +
		labs(x = "relative postion of motif to TSS", y = "percentage") +
		ggtitle(args[1])
png(filename = paste(args[1], "promotermotif.hist.png", sep = "."), width = 1000, height = 800)
print(myplot)
dev.off()

