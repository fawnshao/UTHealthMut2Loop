library(ggplot2)
library(pheatmap)

args <- c("housekeepinggene.GTRD.count.txt", 
          "tissuespecificgene.GTRD.count.txt") #, 
          # "GTEx.100way.tab")
hkg <- read.table(args[1], sep = "\t", header = F)
tsg <- read.table(args[2], sep = "\t", header = F)
# all <- fread(args[3], sep = "\t", header = F)

makematrix <- function(x) {
	col.lists <- unique(x[,1])
	row.lists <- unique(x[,2])
	out.matrix <- matrix(nrow = length(row.lists), ncol = length(col.lists))
	for(i in 1:nrow(out.matrix)){
		for(j in 1:ncol(out.matrix)){
			a <- x[x[,1]==col.lists[j] & x[,2]==row.lists[i],3]
			if(length(a) == 1){
				out.matrix[i,j] <- a
			}
		}
	}
	colnames(out.matrix) <- col.lists
	rownames(out.matrix) <- row.lists
	return(out.matrix)
}
hkg.score <- matrix(hkg)
hkg.score <- matrix(hkg)

# all.dis <- unique(all[,c(1,6)])

colorn <- 10
colors <- colorRampPalette(c("white", "blue"))(colorn)

data <- hkg.score
png(filename = "hkg.GTRD.peaks.png", width = 1500, height = 1200)
myplot <- pheatmap(data, scale = "none", 
	show_rownames = F, show_colnames = F, color = colors, 
	cluster_cols = F, cluster_rows = T)
dev.off()
cluster <- cutree(myplot$tree_row, k = 5)
write.table(data.frame(cluster[myplot$tree_row$order], 
	rownames(data)[myplot$tree_row$order],
	data[myplot$tree_row$order, ]), 
	file = "hkg.GTRD.peaks.tsv", 
	sep = "\t", row.names = FALSE, quote = FALSE)

data <- tsg.score
png(filename = "tsg.GTRD.peaks.png", width = 1500, height = 1200)
myplot <- pheatmap(data, scale = "none", 
	show_rownames = F, show_colnames = F, color = colors, 
	cluster_cols = F, cluster_rows = T)
dev.off()
cluster <- cutree(myplot$tree_row, k = 5)
write.table(data.frame(cluster[myplot$tree_row$order], 
	rownames(data)[myplot$tree_row$order],
	data[myplot$tree_row$order, ]), 
	file = "tsg.GTRD.peaks.tsv", 
	sep = "\t", row.names = FALSE, quote = FALSE)
save.image("GTRD.RData")
