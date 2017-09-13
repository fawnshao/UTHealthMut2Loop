args <- commandArgs(TRUE)

input <- args[1]
library(ggplot2)

a <- read.table(input, sep = "\t")
coords <- cbind.data.frame(a, a[,3]);
colnames(coords) <- c("x1", "x2", "y1", "y2");
p0 <- ggplot(coords, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2) + 
    geom_rect()
png(filename = paste(input, "rect.png", sep = "."), width = 2000, height = 2000)
print(p0)
dev.off()