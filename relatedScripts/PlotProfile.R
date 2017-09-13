args <- commandArgs(TRUE)

input <- args[1]
library(ggplot2)

coords <- read.table(input, sep = "\t")
colnames(coords) <- c("x1", "x2", "y1")
p0 <- ggplot(coords, aes(xmin = x1, xmax = x2, ymin = 0, ymax = y1, fill = "red")) + 
	geom_rect(colour = "red")
png(filename = paste(input, "rect.png", sep = "."), width = 2000, height = 2000)
print(p0)
dev.off()