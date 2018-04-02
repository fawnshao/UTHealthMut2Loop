library(ggplot2)
library(pheatmap)
library(data.table)

args <- c("hkg.tsg.srtbyPCA.class.neighbors.sim")
input <- fread(args[1], sep = "\t", header = T)


genes <- as.matrix(input[,1])
class <- as.matrix(input[,2])
types <- factor(class)
types.sim <- as.vector(types)
types.sim[1:3958] <- "HKG"
types.sim[3959:7748] <- "singleTSG"

gene.strand <- as.matrix(input[,3])
neighbors.strand <- as.matrix(input[,5])
relations <- gene.strand==neighbors.strand
strand <- relations
strand[relations==T] <- 1
strand[relations==F] <- -1

# convergent.dis <- abs(as.matrix(input[,6]))
# divergent.dis <- abs(as.matrix(input[,6]))
same.dis <- as.matrix(input[,6])
diff.dis <- as.matrix(input[,6])
same.dis[relations==T] <- NA
diff.dis[relations==F] <- NA

data <- data.frame(paste(genes,as.matrix(input[,4])),same.dis, diff.dis)
rownames(data) <- paste(class, rownames(data))
colnames(data) <- c("gene.pairs", "same", "diff")
data.m <- melt(data)
annos <- data.frame(class = types, class.sim = types.sim)
rownames(annos) <- rownames(data)

colors <- colorRampPalette(c("white", "blue"))(100)
data.x <- abs(as.matrix((input[,6])))
data.x[data.x > 1000] <- 1000
data.x <- data.x * strand
png(filename = paste(args[1], "dis.pheatmap.png", sep = "."), width = 1000, height = 1500)
pheatmap(as.matrix(data.x), scale = "none", # annotation_row = annos,
	show_rownames = F, show_colnames = T, color = colors, # gaps_row = c(2552, 5287),
	cluster_cols = F, cluster_rows = F)
dev.off()

png(filename = paste(args[1], "displot.png", sep = "."), width = 1000, height = 2500)
# corrplot(data, method = "circle", is.corr = FALSE)
myplot <- ggplot(data = data.m, aes(x = variable, y = gene.pairs)) + 
		geom_point(aes(colour = value)) +
		ggtitle(args[1])
print(myplot)
dev.off()

