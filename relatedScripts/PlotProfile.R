args <- commandArgs(TRUE)

input <- args[1]
library(ggplot2)

a <- read.table(input, sep = "\t")
a[,3] <- log10(a[,3])
coords <- a[order(a[,2] - a[,1], decreasing = TRUE),]
colnames(coords) <- c("x1", "x2", "y1")
p0 <- ggplot(coords, aes(xmin = x1, xmax = x2, ymin = 0, ymax = y1, fill = "red")) + 
	geom_rect(colour = "black")
png(filename = paste(input, "rect.png", sep = "."), width = 2000, height = 2000)
print(p0)
dev.off()
