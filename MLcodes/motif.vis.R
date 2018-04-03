library(ggplot2)
library(pheatmap)
library(data.table)

# args <- c("known.motif.mat")
args <- c("sim.known.motif.mat")
input <- fread(args[1], sep = "\t", header = T, na.strings = "/")

motifs <- as.matrix(input[,1])
motifs.seq <- as.matrix(input[,2])
# qv <- data.matrix(input[,c(3,5,7,9,11,13,15)])
# rv <- data.matrix(input[,c(4,6,8,10,12,14,16)])
qv <- data.matrix(input[,c(3,5,7)])
rv <- data.matrix(input[,c(4,6,8)])
class <- gsub(".q", "", colnames(qv))
colnames(qv) <- class

# a <- data.frame(motifs[!is.na(qv[,5]),], qv[!is.na(qv[,5]),])
# b <- data.frame(motifs[!is.na(qv[,5]),], rv[!is.na(qv[,5]),])
labels <- paste(motifs, motifs.seq, sep = ":")
labels <- gsub("/Homer", "", labels)
labels <- factor(labels, levels = labels)
# a <- data.frame(labels[rv[,1]-rv[,2]>15], qv[rv[,1]-rv[,2]>15,])
# b <- data.frame(labels[rv[,1]-rv[,2]>15], rv[rv[,1]-rv[,2]>15,])
a <- data.frame(labels[-grep("DAP-Seq", labels)], qv[-grep("DAP-Seq", labels),])
b <- data.frame(labels[-grep("DAP-Seq", labels)], rv[-grep("DAP-Seq", labels),])
data <- data.frame(melt(a), melt(b))
colnames(data) <- c("labels", "class", "qvalue", "labels.1", "class.1", "percentage")

png(filename = paste(args[1], "enrichment.plot.png", sep = "."), width = 700, height = 500)
# corrplot(data, method = "circle", is.corr = FALSE)
myplot <- ggplot(data = data, aes(x = class, y = labels)) + 
		geom_point(aes(colour = qvalue, size = percentage)) +
		scale_colour_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0.5, space = "Lab", na.value = "grey50", guide = "colourbar") + 
		ggtitle(args[1])
print(myplot)
dev.off()



################################
library(ggplot2)
library(pheatmap)
library(data.table)

args <- c("hkg.vert.known.ETS.Sp1.YY1.bin.txt", 
	"mixTSG.vert.known.ETS.Sp1.YY1.bin.txt", 
	"singleTSG.vert.known.ETS.Sp1.YY1.bin.txt")
mymotifs <- c("Elk4(ETS)", "ETS(ETS)", "Elk1(ETS)", "ELF1(ETS)", 
	"YY1(Zf)", "Fli1(ETS)", "Sp1(Zf)", "GABPA(ETS)",
	"ETS1(ETS)", "EWS:FLI1-fusion(ETS)", "Etv2(ETS)" 
	)
input1 <- fread(args[1], sep = "\t", header = F, na.strings = "/")
input2 <- fread(args[2], sep = "\t", header = F, na.strings = "/")
input3 <- fread(args[3], sep = "\t", header = F, na.strings = "/")
percentage <- c(data.matrix(input1[,3])[,1]/2551, 
	data.matrix(input2[,3])[,1]/1574, 
	data.matrix(input3[,3])[,1]/2735)
class <- c(rep("hkg",nrow(input1)), rep("mixTSG",nrow(input2)), rep("singleTSG",nrow(input3)))
datax <- data.frame(rbind.data.frame(input1[,1:2], input2[,1:2], input3[,1:2]), class, percentage)
colnames(datax)[1:2] <- c("motif", "bin")

datax <- datax[abs(datax[,2])<11 & datax[,1] %in% mymotifs,]
myplot <- ggplot(data = datax, aes(x = bin, y = percentage, fill = class)) + 
		geom_bar(stat = "identity", position = position_dodge()) +
		facet_grid(motif ~ . + class) +
		scale_x_continuous(limits = c(-10,10)) +
		ggtitle(args[1])
png(filename = paste(args[1], "motif.bar.png", sep = "."), width = 900, height = 800)
print(myplot)
dev.off()


################################
library(ggplot2)
library(pheatmap)
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

class <- c(rep("hkg",nrow(input1)), rep("mixTSG",nrow(input2)), rep("singleTSG",nrow(input3)))
datax <- data.frame(rbind.data.frame(input1, input2, input3), class)
colnames(datax)[1:2] <- c("motif", "dis2tss")

datax <- datax[datax[,1] %in% mymotifs,]
myplot <- ggplot(data = datax, aes(x = dis2tss, colour = class)) + 
		#geom_histogram(fill = class, binwidth = 10)
		geom_freqpoly(binwidth = 10) +
		facet_wrap( ~ motif, nrow = 3) + 
		# facet_grid(motif ~ .) + 
		# facet_grid(motif ~ . + class) +
		scale_x_continuous(limits = c(-1000,1000)) +
		ggtitle(args[1])
png(filename = paste(args[1], "motif.hist.png", sep = "."), width = 800, height = 800)
print(myplot)
dev.off()


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
		geom_freqpoly(binwidth = 50, position = "identity", data = datax[class == 'HKG',], aes(y = ..count.. / 2551)) +
		geom_freqpoly(binwidth = 50, position = "identity", data = datax[class == 'mixTSG',], aes(y = ..count.. / 1574)) +
		geom_freqpoly(binwidth = 50, position = "identity", data = datax[class == 'singleTSG',], aes(y = ..count.. / 2735)) +
		scale_x_continuous(limits = c(-100,20000)) +
		labs(x = "relative postion of CAGE to TSS", y = "percentage") +
		ggtitle(args[1])
png(filename = paste(args[1], "cage.hist.png", sep = "."), width = 1500, height = 500)
print(myplot)
dev.off()

