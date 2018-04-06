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
		geom_freqpoly(binwidth = 50, position = "identity", data = datax[types.sim == "HKG",], aes(y = ..count.. / 2551)) +
		geom_freqpoly(binwidth = 50, position = "identity", data = datax[types.sim == "mixTSG",], aes(y = ..count.. / 1574)) +
		geom_freqpoly(binwidth = 50, position = "identity", data = datax[types.sim == "singleTSG",], aes(y = ..count.. / 2735)) +
		facet_wrap( ~ motif, nrow = 2) + 
		scale_x_continuous(limits = c(-1000,5000)) +
		labs(x = "relative postion of motif to TSS", y = "percentage") +
		ggtitle(args[1])
png(filename = paste(args[1], "promotermotif.hist.png", sep = "."), width = 1000, height = 800)
print(myplot)
dev.off()

########
### ETS SP1 YY1
library(ggplot2)
library(data.table)

args <- c("hkg.vert.known.ETS.Sp1.YY1.dis.txt", 
	"mixTSG.vert.known.ETS.Sp1.YY1.dis.txt", 
	"singleTSG.vert.known.ETS.Sp1.YY1.dis.txt")
mymotifs <- c("Elk4(ETS)", "ETS(ETS)", "Elk1(ETS)", "ELF1(ETS)", 
	"YY1(Zf)", "Fli1(ETS)", "Sp1(Zf)", "GABPA(ETS)",
	"ETS1(ETS)", "EWS:FLI1-fusion(ETS)", "Etv2(ETS)", "PU.1(ETS)"
	)
input1 <- fread(args[1], sep = "\t", header = F, na.strings = "/")
input2 <- fread(args[2], sep = "\t", header = F, na.strings = "/")
input3 <- fread(args[3], sep = "\t", header = F, na.strings = "/")

types.sim <- c(rep("HKG",nrow(input1)), rep("mixTSG",nrow(input2)), rep("singleTSG",nrow(input3)))
datax <- data.frame(rbind.data.frame(input1, input2, input3), types.sim)
colnames(datax) <- c("motif", "dis2tss", "class")
datax <- datax[datax[,1] %in% mymotifs,]

myplot <- ggplot(data = datax, aes(x = dis2tss, colour = class)) + 
		geom_freqpoly(binwidth = 10, position = "identity", data = datax[datax[,3] == "HKG",], aes(y = ..count.. / 2551)) +
		geom_freqpoly(binwidth = 10, position = "identity", data = datax[datax[,3] == "mixTSG",], aes(y = ..count.. / 1574)) +
		geom_freqpoly(binwidth = 10, position = "identity", data = datax[datax[,3] == "singleTSG",], aes(y = ..count.. / 2735)) +
		facet_wrap( ~ motif, ncol = 4) + 
		scale_x_continuous(limits = c(-1000,1000)) +
		labs(x = "relative postion of motif to TSS", y = "percentage") +
		ggtitle(args[1])
png(filename = paste(args[1], "promoterETSmotif.hist.png", sep = "."), width = 800, height = 800)
print(myplot)
dev.off()

########
### neighboring TSS
library(ggplot2)
library(data.table)

args <- c("hkg.tsg.srtbyPCA.class.tss.neighbors.sim")
input <- fread(args[1], sep = "\t", header = T, na.strings = "/")
types.sim <- as.matrix(input[,2])[,1]
types.sim[grep("hkg",types.sim)] <- "HKG"
types.sim[types.sim!="HKG" & types.sim!="mixTSG"] <- "singleTSG"
strand.lab <- as.matrix(input[,3])[,1]
strand.lab[input[,3]==input[,5]] <- "same"
strand.lab[input[,3]!=input[,5]] <- "diff"

datax <- data.frame(input[,6], strand.lab, types.sim)
colnames(datax) <- c("tss.dis", "strand", "class")

myplot <- ggplot(data = datax, aes(x = tss.dis, colour = class)) + 
		geom_freqpoly(binwidth = 30, position = "identity", data = datax[types.sim == "HKG",], aes(y = ..count.. / 2551)) +
		geom_freqpoly(binwidth = 30, position = "identity", data = datax[types.sim == "mixTSG",], aes(y = ..count.. / 1574)) +
		geom_freqpoly(binwidth = 30, position = "identity", data = datax[types.sim == "singleTSG",], aes(y = ..count.. / 2735)) +
		scale_x_continuous(limits = c(-3000,3000)) +
		facet_wrap( ~ strand, nrow = 1) + 
		labs(x = "nearest neighboring TSS distance", y = "percentage") +
		ggtitle(args[1])
png(filename = paste(args[1], "neighboringTSS.hist.png", sep = "."), width = 1000, height = 500)
print(myplot)
dev.off()

########
### neighboring TSS annotation
library(ggplot2)
library(data.table)

args <- c("hkg.tsg.srtbyPCA.class.tss.neighbors.sim.annotation")
input <- fread(args[1], sep = "\t", header = T, na.strings = "/")
types.sim <- as.matrix(input[,2])[,1]
types.sim[grep("hkg",types.sim)] <- "HKG"
types.sim[types.sim!="HKG" & types.sim!="mixTSG"] <- "singleTSG"

anno.class <- as.matrix(input)[,7]
anno.class[anno.class!="protein_coding"] <- "others"
anno.class.1 <- anno.class
anno.class <- as.matrix(input)[,8]
anno.class[anno.class!="protein_coding"] <- "others"
# anno.class[anno.class!="protein_coding" & anno.class!="lincRNA" & anno.class!="antisense" & anno.class!="processed_transcript" & anno.class!="pseudogene" & anno.class!="snoRNA" & anno.class!="snRNA" & anno.class!="sense_overlapping" & anno.class!="sense_intronic"] <- "others"
anno.class.2 <- anno.class
annotation.type <- paste(anno.class.1, anno.class.2, sep = ": ")
# annotation.type <- paste(as.matrix(input)[,7], as.matrix(input)[,8], sep = ":")
x <- melt(sort(table(annotation.type[types.sim=="HKG" & abs(input[,6]) < 1000]), decreasing = F) / 2551)
y <- melt(sort(table(annotation.type[types.sim=="singleTSG" & abs(input[,6]) < 1000]), decreasing = F) / 2735)
z <- melt(sort(table(annotation.type[types.sim=="mixTSG" & abs(input[,6]) < 1000]), decreasing = F) / 1574)
datax <- data.frame(rbind.data.frame(x,y,z),c(rep("HKG", nrow(x)),rep("singleTSG", nrow(y)),rep("mixTSG", nrow(z))))
colnames(datax) <- c("type", "value", "class")
# sums <- c(2551, 1574, 2735)
# names(sums) <- c("HKG", "mixTSG", "singleTSG")


myplot <- ggplot(data = datax, aes(x = class, y = value, fill = type)) + 
	# geom_bar(position = position_dodge(1), stat = "identity") + 
	geom_bar(stat = "identity", width = 0.5) + 
	# geom_text(aes(y = label_ypos, label = value), vjust=1.6, color = "white", size = 3.5)+
	ggtitle(args[1])# + coord_flip()
	# theme(axis.text.x = element_text(angle = 60, hjust = 1))
# png(filename = paste(args[1], "annotype.barplot.png", sep = "."), width = 1200, height = 600)
png(filename = paste(args[1], "annotype.1k.barplot.png", sep = "."), width = 500, height = 800)
print(myplot)
dev.off()

anno.class <- as.matrix(input)[,7]
anno.class[anno.class!="protein_coding" & anno.class!="lincRNA" & anno.class!="antisense"] <- "others"
anno.class.1 <- anno.class
anno.class <- as.matrix(input)[,8]
anno.class[anno.class!="protein_coding" & anno.class!="lincRNA" & anno.class!="antisense"] <- "others"
anno.class.2 <- anno.class
annotation.type <- paste(anno.class.1, anno.class.2, sep = ": ")
x <- melt(sort(table(annotation.type[types.sim=="HKG" & abs(input[,6]) < 1000]), decreasing = T) / 2551)
y <- melt(sort(table(annotation.type[types.sim=="singleTSG" & abs(input[,6]) < 1000]), decreasing = T) / 2735)
z <- melt(sort(table(annotation.type[types.sim=="mixTSG" & abs(input[,6]) < 1000]), decreasing = T) / 1574)
datax <- data.frame(rbind.data.frame(x,y,z),c(rep("HKG", nrow(x)),rep("singleTSG", nrow(y)),rep("mixTSG", nrow(z))))
colnames(datax) <- c("type", "value", "class")
myplot <- ggplot(data = datax, aes(x = type, y = value, fill = class)) + 
	geom_bar(stat = "identity", width = 0.5) + 
	ggtitle(args[1]) + coord_flip()
png(filename = paste(args[1], "annotype.1k.t.barplot.png", sep = "."), width = 800, height = 600)
print(myplot)
dev.off()

myplot <- ggplot(data = datax, aes(x = class, y = type)) + 
		geom_point(aes(colour = value, size = value)) +
		ggtitle(args[1])
png(filename = paste(args[1], "annotype.1k.t.dotplot.png", sep = "."), width = 400, height = 600)
print(myplot)
dev.off()

