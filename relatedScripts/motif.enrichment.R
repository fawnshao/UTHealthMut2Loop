library(data.table)
library(ggplot2)
# library(reshape2)
library(pheatmap)
args <- c("housekeepinggene.known.motifs.mat", 
          "tissuespecificgene.known.motifs.mat")
# hkg <- read.table(args[1], sep = "\t", header = T, row.names = 1)
hkg <- fread(args[1], sep = "\t", header = T)
tsg <- fread(args[2], sep = "\t", header = T)
hkg.score <- hkg[,-1]
tsg.score <- tsg[,-1]
rownames(hkg.score) <- as.matrix(hkg[,1])
rownames(tsg.score) <- as.matrix(tsg[,1])

colorn <- 100
colors <- colorRampPalette(c("blue", "white", "red"))(colorn)

data <- log2(data.matrix(hkg.score) + 1)
png(filename = "hkg.homer.motif.png", width = 1500, height = 1200)
myplot <- pheatmap(data, scale = "column", 
	show_rownames = F, show_colnames = F, color = colors, 
	cluster_cols = F, cluster_rows = T)
dev.off()
cluster <- cutree(myplot$tree_row, k = 5)
write.table(data.frame(cluster[myplot$tree_row$order], 
	rownames(data)[myplot$tree_row$order],
	data[myplot$tree_row$order, ]), 
	file = "hkg.homer.motif.tsv", 
	sep = "\t", row.names = FALSE, quote = FALSE)

data <- log2(data.matrix(tsg.score) + 1)
png(filename = "tsg.homer.motif.png", width = 1500, height = 1200)
myplot <- pheatmap(data, scale = "column", 
	show_rownames = F, show_colnames = F, color = colors, 
	cluster_cols = F, cluster_rows = T)
dev.off()
cluster <- cutree(myplot$tree_row, k = 5)
write.table(data.frame(cluster[myplot$tree_row$order], 
	rownames(data)[myplot$tree_row$order],
	data[myplot$tree_row$order, ]), 
	file = "tsg.homer.motif.tsv", 
	sep = "\t", row.names = FALSE, quote = FALSE)
save.image("homer.motif.RData")

##################
args <- c("housekeepinggene.known.motifs.txt", 
          "tissuespecificgene.known.motifs.txt") #, 
          # "GTEx.100way.tab")
hkg <- fread(args[1], sep = "\t", header = T)
tsg <- fread(args[2], sep = "\t", header = T)
# all <- fread(args[3], sep = "\t", header = F)

hkg.score <- hkg[,-c(1:21)]
tsg.score <- tsg[,-c(1:21)]
# hkg.score[!is.na(hkg.score)] <- 1
# hkg.score[is.na(hkg.score)] <- 0
# tsg.score[!is.na(tsg.score)] <- 1
# tsg.score[is.na(tsg.score)] <- 0
# fruit <- c("apple", "banana", "pear", "pinapple")
# length(fruit[str_detect(fruit, "b")])
# length(strsplit(as.matrix(hkg.score[250,25]),"),")[[1]])
for(i in 1:nrow(hkg.score)){
	for(j in 1:ncol(hkg.score)){
		if(is.na(hkg.score[i,j,with=F])){
			hkg.score[i,j] <- 0
		}
		else{
			hkg.score[i,j] <- length(strsplit(as.character(hkg.score[i,j,with=F]),"),")[[1]])
		}
	}
}
for(i in 1:nrow(tsg.score)){
	for(j in 1:ncol(tsg.score)){
		if(is.na(tsg.score[i,j,with=F])){
			tsg.score[i,j] <- 0
		}
		else{
			tsg.score[i,j] <- length(strsplit(as.character(tsg.score[i,j,with=F]),"),")[[1]])
		}
	}
}
rownames(hkg.score) <- as.matrix(hkg[,1])
rownames(tsg.score) <- as.matrix(tsg[,1])
write.csv(hkg.score, "hkg.homer.motifs.occurence.csv")
write.csv(tsg.score, "tsg.homer.motifs.occurence.csv")
# all.dis <- unique(all[,c(1,6)])

colorn <- 10
colors <- colorRampPalette(c("white", "blue"))(colorn)

data <- data.matrix(hkg.score)
png(filename = "hkg.homer.motif.png", width = 1500, height = 1200)
myplot <- pheatmap(data, scale = "none", 
	show_rownames = F, show_colnames = F, color = colors, 
	cluster_cols = F, cluster_rows = T)
dev.off()
cluster <- cutree(myplot$tree_row, k = 5)
write.table(data.frame(cluster[myplot$tree_row$order], 
	rownames(data)[myplot$tree_row$order],
	data[myplot$tree_row$order, ]), 
	file = "hkg.homer.motif.tsv", 
	sep = "\t", row.names = FALSE, quote = FALSE)

data <- data.matrix(tsg.score)
png(filename = "tsg.homer.motif.png", width = 1500, height = 1200)
myplot <- pheatmap(data, scale = "none", 
	show_rownames = F, show_colnames = F, color = colors, 
	cluster_cols = F, cluster_rows = T)
dev.off()
cluster <- cutree(myplot$tree_row, k = 5)
write.table(data.frame(cluster[myplot$tree_row$order], 
	rownames(data)[myplot$tree_row$order],
	data[myplot$tree_row$order, ]), 
	file = "tsg.homer.motif.tsv", 
	sep = "\t", row.names = FALSE, quote = FALSE)
save.image("homer.motif.RData")

#############################################################
args <- c("housekeepinggene.homer.motifs/knownResults.txt", 
          "tissuespecificgene.homer.motifs/knownResults.txt") #, 
          # "GTEx.100way.tab")
hkg <- read.table(args[1], sep = "\t", header = F, skip = 1)
tsg <- read.table(args[2], sep = "\t", header = F, skip = 1)

# hkg.sig <- as.matrix(hkg[hkg[,5] < 0.05, c(7,9)])
# tsg.sig <- as.matrix(tsg[tsg[,5] < 0.05, c(7,9)])
# rownames(hkg.sig) <- hkg[hkg[,5] < 0.05, 1]
# rownames(tsg.sig) <- tsg[tsg[,5] < 0.05, 1]
hkg.sig <- as.matrix(hkg[1:30, c(7,9)])
tsg.sig <- as.matrix(tsg[1:30, c(7,9)])
rownames(hkg.sig) <- hkg[1:30, 1]
rownames(tsg.sig) <- tsg[1:30, 1]
hkg.sig[,1] <- as.numeric(unlist(strsplit(hkg.sig[,1],"%")))
hkg.sig[,2] <- as.numeric(unlist(strsplit(hkg.sig[,2],"%")))
tsg.sig[,1] <- as.numeric(unlist(strsplit(tsg.sig[,1],"%")))
tsg.sig[,2] <- as.numeric(unlist(strsplit(tsg.sig[,2],"%")))

# data <- data.frame(c(rownames(hkg.sig),rownames(tsg.sig)),
data <- data.frame(c(t(as.data.frame(strsplit(rownames(hkg.sig),"\\/")))[,1],
					t(as.data.frame(strsplit(rownames(tsg.sig),"\\/")))[,1]),
				c(as.numeric(unlist(strsplit(hkg.sig[,1],"%"))), 
					as.numeric(unlist(strsplit(tsg.sig[,1],"%")))), 
				c(rep("hkg",nrow(hkg.sig)),rep("tsg",nrow(tsg.sig))))
colnames(data) <- c("TF", "gene", "class")
# png(filename = "enriched.homer.motifs.png", width = 1800, height = 600)
data$TF <- factor(data$TF, levels=unique(data$TF))
png(filename = "top30enriched.homer.motifs.png", width = 1500, height = 600)
ggplot (data = data, aes(x = TF, y = gene, fill = class)) + 
	geom_bar(stat="identity", position=position_dodge(), width = 0.6) +
	theme(axis.text.x = element_text(angle = 60, hjust = 1)) 
dev.off()

# strsplit(rownames(hkg.sig),"\\(")



####

total <- as.matrix(read.table("hkg.motif.id.tsv", sep = "\t", header = T, row.names = 1))
values <- data.frame(t(as.data.frame(strsplit(rownames(total),"\\/")))[,1],
	as.numeric(unlist(strsplit(total[,1],"%"))),
	as.numeric(unlist(strsplit(total[,2],"%"))),
	as.numeric(unlist(strsplit(total[,3],"%")))
	)
colnames(values) <- c("TF", "GTEx", "hkg","tsg")
x <- values[1:30,]
# png(filename = "enriched.homer.motifs.png", width = 1800, height = 600)
data <- melt(x)
data$TF <- factor(data$TF, levels=unique(x$TF))
png(filename = "HKGtop30enriched.homer.motifs.png", width = 1500, height = 600)
ggplot (data = data, aes(x = TF, y = value, fill = variable)) + 
	geom_bar(stat="identity", position=position_dodge(), width = 0.6) +
	theme(axis.text.x = element_text(angle = 60, hjust = 1)) 
dev.off()
