library(data.table)
library(ggplot2)

args <- c("housekeepinggene.100way.tab", 
          "tissuespecificgene.100way.tab", 
          "GTEx.100way.tab")
hkg <- fread(args[1], sep = "\t", header = F)
tsg <- fread(args[2], sep = "\t", header = F)
all <- fread(args[3], sep = "\t", header = F)

hkg.dis <- unique(hkg[,c(1,6)])
tsg.dis <- unique(tsg[,c(1,6)])
all.dis <- unique(all[,c(1,6)])

colnames(hkg.dis) <- c("gene", "score")
colnames(tsg.dis) <- c("gene", "score")
colnames(all.dis) <- c("gene", "score")

dis.abs <- data.frame(rbind(hkg.dis[,2], tsg.dis[,2], all.dis[,2]),
                      c(rep("hkg", nrow(hkg.dis)),rep("tsg", nrow(tsg.dis)), rep("all", nrow(all.dis))))
colnames(dis.abs) <- c("score", "group")
cdf <- ggplot (data = dis.abs, aes(x = score, group = group, color = group)) + 
    stat_ecdf() + scale_x_continuous(limits = c(0, 1))
png(filename = "proteincoding.noMT.100way.phastCons.ecdf.1.png", width = 600, height = 600)
cdf
dev.off()
cdf <- ggplot (data = dis.abs, aes(x = score, group = group, color = group)) + 
    stat_ecdf() + scale_x_continuous(limits = c(0.25, 1))
png(filename = "proteincoding.noMT.100way.phastCons.ecdf.2.png", width = 600, height = 600)
cdf
dev.off()

cdf <- ggplot (data = dis.abs, aes(x = -1 * score, group = group, color = group)) + 
    stat_ecdf() + scale_x_continuous(limits = c(-1, -0.25)) # + scale_x_reverse()
png(filename = "proteincoding.noMT.100way.phastCons.revecdf.png", width = 600, height = 600)
cdf
dev.off()

ks.test(as.matrix(hkg.dis[,2]),as.matrix(all.dis[,2]))
png(filename = "proteincoding.noMT.100way.phastCons.violin.png", width = 600, height = 600)
ggplot (data = dis.abs, aes(x = group, y = score, fill = group)) + 
    geom_violin() + stat_summary(fun.data = "mean_sdl", geom = "pointrange", color = "black")
dev.off()
