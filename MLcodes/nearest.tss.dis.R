library(data.table)
library(ggplot2)
args <- commandArgs(TRUE)
args <- c("hkg.tsg.srtbyPCA.class.tss.neighbors.sim")
input <- fread(args[1], sep = "\t", header = T)

genes <- as.matrix(input[,1])
class <- as.matrix(input[,2])
types <- factor(class)
types.sim <- as.vector(types)
types.sim[1:2561] <- "HKG"
types.sim[2562:5303] <- "singleTSG"


divergent <- input[Gene.strand!=neighbors.strand & distance <= 0,]
convergent <- input[Gene.strand!=neighbors.strand & distance > 0,]
uphkg <- input[Gene.strand==neighbors.strand & distance <= 0,]
downhkg <- input[Gene.strand==neighbors.strand & distance > 0,]

strands <- types.sim
strands[input[,3]!=input[,5] & input[,6] <= 0] <- "divergent"
strands[input[,3]!=input[,5] & input[,6] > 0] <- "convergent"
strands[input[,3]==input[,5] & input[,6] <= 0] <- "neighbors.up"
strands[input[,3]==input[,5] & input[,6] > 0] <- "neighbors.down"

datax <- data.frame(types.sim, strands, abs(input[,6]))
colnames(datax) <- c("type", "strands", "value")
ymax <- 10000
myplot <- ggplot(data = datax, aes(x = type, y = value, fill = type)) + 
	geom_boxplot(position = position_dodge(1), outlier.alpha = 0.1, outlier.size = 0.1) + 
	scale_y_continuous(limits = c(0, ymax)) +
	ggtitle(args[1]) + 
	facet_grid(. ~ strands) +
	theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1))
png(filename = paste(args[1], "dis.strand.boxplot.png", sep = "."), width = 1200, height = 600)
print(myplot)
dev.off()

datax <- data.frame(types, abs(input[,6]))
colnames(datax) <- c("type", "value")
# ymax <- quantile(datax$value,probs = 0.95)
ymax <- 20000
myplot <- ggplot(data = datax, aes(x = type, y = value, fill = type)) + 
	geom_boxplot(position = position_dodge(1), outlier.alpha = 0.1, outlier.size = 0.1) + 
	scale_y_continuous(limits = c(0, ymax)) +
	ggtitle(args[1]) + 
	theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1))
png(filename = paste(args[1], "dis.boxplot.png", sep = "."), width = 1000, height = 600)
print(myplot)
dev.off()


datax <- data.frame(types.sim, abs(input[,6]))
colnames(datax) <- c("type", "value")
# ymax <- quantile(datax$value,probs = 0.9)
ymax <- 10000
myplot <- ggplot(data = datax, aes(x = type, y = value, fill = type)) + 
	geom_boxplot(position = position_dodge(1), outlier.alpha = 0.1, outlier.size = 0.1) + 
	scale_y_continuous(limits = c(0, ymax)) +
	ggtitle(args[1]) + 
	theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1))
png(filename = paste(args[1], "dis.boxplot.sim.png", sep = "."), width = 400, height = 600)
print(myplot)
dev.off()

myplot <- ggplot(data = datax, aes(x = type, y = value, fill = type)) + 
	geom_violin(position = position_dodge(1)) + 
	scale_y_continuous(limits = c(0, ymax)) +
	ggtitle(args[1]) + 
	theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1))
png(filename = paste(args[1], "dis.vioplot.sim.png", sep = "."), width = 400, height = 600)
print(myplot)
dev.off()
