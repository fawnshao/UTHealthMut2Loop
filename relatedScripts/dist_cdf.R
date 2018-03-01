library(data.table)
library(ggplot2)

args <- c("housekeepinggene.byGTEx.neighbors.txt", 
          "tissuespecificgene.byGTEx.neighbors.txt", 
          "GTEx.neighbors.txt")
          # "gencode.neighbors.txt")
hkg <- fread(args[1], sep = "\t", header = F)
tsg <- fread(args[2], sep = "\t", header = F)
all <- fread(args[3], sep = "\t", header = F)

hkg.dis <- unique(hkg[,c(4,13)])
tsg.dis <- unique(tsg[,c(4,13)])
all.dis <- unique(all[,c(4,13)])

colnames(hkg.dis) <- c("gene", "dis")
colnames(tsg.dis) <- c("gene", "dis")
colnames(all.dis) <- c("gene", "dis")

dis.abs <- data.frame(rbind(abs(hkg.dis[,2]), abs(tsg.dis[,2]), abs(all.dis[,2])),
                      c(rep("hkg", nrow(hkg.dis)),rep("tsg", nrow(tsg.dis)), rep("all", nrow(all.dis))))
colnames(dis.abs) <- c("distance", "group")
cdf <- ggplot (data = dis.abs, aes(x = distance, group = group, color = group)) + 
    stat_ecdf() + scale_x_continuous(limits = c(0, 100000))
png(filename = "proteincoding.noMT.dis.ecdf.png", width = 1200, height = 600)
cdf
dev.off()

dis.abs.log2 <- data.frame(log2(dis.abs[,1] + 1), dis.abs[,2])
colnames(dis.abs.log2) <- c("distance", "group")
cdf <- ggplot (data = dis.abs.log2, aes(x = distance, group = group, color = group)) + 
    stat_ecdf()
png(filename = "proteincoding.noMT.log2dis.ecdf.png", width = 1200, height = 600)
cdf
dev.off()

dis.abs <- data.frame(rbind(abs(hkg.dis[,2]), abs(tsg.dis[,2]), abs(all.dis[sample(1:nrow(all.dis), 2000),2])),
                      c(rep("hkg", nrow(hkg.dis)),rep("tsg", nrow(tsg.dis)), rep("rand", 2000)))
colnames(dis.abs) <- c("distance", "group")
cdf <- ggplot (data = dis.abs[dis.abs[,1] > 0,], aes(x = distance, group = group, color = group)) + 
    stat_ecdf() + scale_x_continuous(limits = c(0, 100000))
png(filename = "proteincoding.noMT.nooverlap.dis.ecdf.png", width = 1200, height = 600)
cdf
dev.off()

a <- ecdf(as.matrix(abs(hkg.dis[,2])))
b <- ecdf(as.matrix(abs(all.dis[,2])))
ks.test(as.matrix(abs(hkg.dis[,2])), as.matrix(abs(all.dis[,2])))


#####
args <- c("housekeepinggene.byGTEx.neighbors.txt", 
          "tissuespecificgene.byGTEx.neighbors.txt", 
          "GTEx.neighbors.txt",
          "gencode.neighbors.txt")
hkg <- fread(args[1], sep = "\t", header = F)
tsg <- fread(args[2], sep = "\t", header = F)
pcg <- fread(args[3], sep = "\t", header = F)
all <- fread(args[4], sep = "\t", header = F)

hkg.dis <- unique(hkg[,c(4,13)])
tsg.dis <- unique(tsg[,c(4,13)])
pcg.dis <- unique(pcg[,c(4,13)])
all.dis <- unique(all[,c(4,13)])

colnames(hkg.dis) <- c("gene", "dis")
colnames(tsg.dis) <- c("gene", "dis")
colnames(pcg.dis) <- c("gene", "dis")
colnames(all.dis) <- c("gene", "dis")

dis.abs <- data.frame(rbind(abs(hkg.dis[,2]), abs(tsg.dis[,2]), abs(pcg.dis[,2]), abs(all.dis[,2])),
                      c(rep("hkg", nrow(hkg.dis)),rep("tsg", nrow(tsg.dis)), rep("pcg", nrow(pcg.dis)),rep("all", nrow(all.dis))))
colnames(dis.abs) <- c("distance", "group")
cdf <- ggplot (data = dis.abs, aes(x = distance, group = group, color = group)) + 
    stat_ecdf() + scale_x_continuous(limits = c(0, 100000))
png(filename = "4groups.dis.ecdf.png", width = 1200, height = 600)
cdf
dev.off()

cdf <- ggplot (data = dis.abs[dis.abs[,1] > 0,], aes(x = distance, group = group, color = group)) + 
    stat_ecdf() + scale_x_continuous(limits = c(0, 100000))
png(filename = "4groups.nooverlap.dis.ecdf.png", width = 1200, height = 600)
cdf
dev.off()

######
args <- c("housekeepinggene.25k.cluster.count", 
          "tissuespecificgene.25k.cluster.count", 
          "GTEx.25k.cluster.count")
          # "gencode.25k.cluster.count")
hkg <- fread(args[1], header = F)
tsg <- fread(args[2], header = F)
all <- fread(args[3], header = F)

hkg.dis <- hkg[,c(2,1)]
tsg.dis <- tsg[,c(2,1)]
all.dis <- all[,c(2,1)]

colnames(hkg.dis) <- c("cluster", "dis")
colnames(tsg.dis) <- c("cluster", "dis")
colnames(all.dis) <- c("cluster", "dis")

# dis.abs <- data.frame(rbind(-1 * hkg.dis[,2], -1 * tsg.dis[,2], -1 * all.dis[,2]),
#                       c(rep("hkg", nrow(hkg.dis)),rep("tsg", nrow(tsg.dis)), rep("all", nrow(all.dis))))
# dis.abs <- data.frame(rbind(1/hkg.dis[,2], 1/tsg.dis[,2], 1/all.dis[,2]),
#                       c(rep("hkg", nrow(hkg.dis)),rep("tsg", nrow(tsg.dis)), rep("all", nrow(all.dis))))
dis.abs <- data.frame(rbind(hkg.dis[,2], tsg.dis[,2], all.dis[,2]),
                      c(rep("hkg", nrow(hkg.dis)),rep("tsg", nrow(tsg.dis)), rep("all", nrow(all.dis))))
colnames(dis.abs) <- c("cluster.count", "group")
cdf <- ggplot (data = dis.abs, aes(x = cluster.count, group = group, color = group)) + 
    stat_ecdf(geom = "step", pad = FALSE) + scale_x_reverse()
png(filename = "25k.cluster.count.ecdf.png", width = 1200, height = 600)
cdf
dev.off()

cdf <- ggplot (data = dis.abs, aes(x = cluster.count, group = group, color = group)) + 
    stat_ecdf(geom = "step", pad = FALSE) + scale_x_continuous(limits = c(0,100)) 
png(filename = "25k.cluster.count.ecdf.1.png", width = 1200, height = 600)
cdf
dev.off()
ks.test(as.matrix(abs(hkg.dis[,2])), as.matrix(abs(all.dis[,2])))


histg <- ggplot (data = dis.abs, aes(x = cluster.count, group = group, fill = group)) + 
	geom_histogram(binwidth = 1, position = "dodge") + xlab("gene count in a cluster") +
		facet_grid(~group) + scale_x_continuous(limits = c(0,25))
png(filename = "25k.cluster.count.histg.png", width = 1200, height = 400)
histg
dev.off()

######
args <- c("housekeepinggene.10k.cluster.count", 
          "tissuespecificgene.10k.cluster.count", 
          "GTEx.10k.cluster.count")
          # "gencode.10k.cluster.count")
hkg <- fread(args[1], header = F)
tsg <- fread(args[2], header = F)
all <- fread(args[3], header = F)

hkg.dis <- hkg[,c(2,1)]
tsg.dis <- tsg[,c(2,1)]
all.dis <- all[,c(2,1)]

colnames(hkg.dis) <- c("cluster", "dis")
colnames(tsg.dis) <- c("cluster", "dis")
colnames(all.dis) <- c("cluster", "dis")

dis.abs <- data.frame(rbind(hkg.dis[,2], tsg.dis[,2], all.dis[,2]),
                      c(rep("hkg", nrow(hkg.dis)),rep("tsg", nrow(tsg.dis)), rep("all", nrow(all.dis))))
colnames(dis.abs) <- c("cluster.count", "group")
cdf <- ggplot (data = dis.abs, aes(x = cluster.count, group = group, color = group)) + 
    stat_ecdf(geom = "step", pad = FALSE) + scale_x_reverse()
png(filename = "10k.cluster.count.ecdf.png", width = 1200, height = 600)
cdf
dev.off()

cdf <- ggplot (data = dis.abs, aes(x = cluster.count, group = group, color = group)) + 
    stat_ecdf(geom = "step", pad = FALSE) + scale_x_continuous(limits = c(0,30)) 
png(filename = "10k.cluster.count.ecdf.1.png", width = 1200, height = 600)
cdf
dev.off()
ks.test(as.matrix(abs(hkg.dis[,2])), as.matrix(abs(all.dis[,2])))


histg <- ggplot (data = dis.abs, aes(x = cluster.count, group = group, fill = group)) + 
	geom_histogram(binwidth = 1, position = "dodge") + xlab("gene count in a cluster") +
		facet_grid(~group) + scale_x_continuous(limits = c(0,25))
png(filename = "10k.cluster.count.histg.png", width = 1200, height = 400)
histg
dev.off()
