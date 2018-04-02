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

