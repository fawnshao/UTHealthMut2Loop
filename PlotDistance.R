args <- commandArgs(TRUE)

mut2tss_by_disease <- args[1]
mut2motif_by_disease <- args[2]
motif2tss_by_motif <- args[3]
motifmut2tss_by_disease <- args[4]
outpre <- args[5]
libpath <- args[6]
source(libpath)

library(ggplot2)

a <- read.table(mut2tss_by_disease, sep = "\t")
b <- read.table(mut2motif_by_disease, sep = "\t")
c <- read.table(motif2tss_by_motif, sep = "\t")
d <- read.table(motifmut2tss_by_disease, sep = "\t")
colnames(a) <- c("group", "distance")
colnames(b) <- c("group", "distance")
colnames(c) <- c("group", "distance")
colnames(d) <- c("group", "distance")

p0 <- ggplot(a, aes(x = factor(group), y = distance, fill = factor(group))) + 
    geom_violin() + ylim(c(-1e4, 1e4)) +
    # xlab("") + ylab("distance") +
    labs(title = "Mutation To TSS", x = "", y = "distance") +
    annotate("text", label = unique(a$group), x = unique(a$group), y = 0)
p1 <- ggplot(b, aes(x = factor(group), y = distance, fill = factor(group))) + 
    geom_violin() + ylim(c(-1e4, 1e4)) +
    labs(title = "Mutation To motif", x = "", y = "distance") +
    annotate("text", label = unique(b$group), x = unique(b$group), y = 0)
p2 <- ggplot(c, aes(x = factor(group), y = distance, fill = factor(group))) + 
    geom_violin() + ylim(c(-1e4, 1e4)) +
    labs(title = "Motif To TSS", x = "", y = "distance") +
    annotate("text", label = unique(c$group), x = unique(c$group), y = 0)
p3 <- ggplot(d, aes(x = factor(group), y = distance, fill = factor(group))) + 
    geom_violin() + ylim(c(-1e4, 1e4)) +
    labs(title = "motif-Mutation To TSS", x = "", y = "distance") +
    annotate("text", label = unique(d$group), x = unique(d$group), y = 0)
png(filename = paste(outpre, "violin.png", sep = "."), width = 2000, height = 2000)
multiplot(p0, p1, p2, p3, cols = 1)
dev.off()
pdf(file = paste(outpre, "violin.pdf", sep = "."), width = 30, height = 20)
multiplot(p0, p1, p2, p3, cols = 1)
dev.off()