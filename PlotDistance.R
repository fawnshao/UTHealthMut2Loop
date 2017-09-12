args <- commandArgs(TRUE)

mut2tss_by_disease <- args[1]
mut2motif_by_disease <- args[2]
motif2tss_by_motif <- args[3]
motifmut2tss_by_disease <- args[4]
outpre <- args[5]
libpath <- args[6]
source(libpath)

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
    annotate("text", label = labels, x = unique(a$group), y = 0)
p1 <- ggplot(a, aes(x = factor(group), y = distance, fill = factor(group))) + 
    geom_violin() + ylim(c(-1e4, 1e4)) +
    labs(title = "Mutation To TSS", x = "", y = "distance") +
    annotate("text", label = labels, x = unique(a$group), y = 0)
p2 <- ggplot(a, aes(x = factor(group), y = distance, fill = factor(group))) + 
    geom_violin() + ylim(c(-1e4, 1e4)) +
    labs(title = "Mutation To TSS", x = "", y = "distance") +
    annotate("text", label = labels, x = unique(a$group), y = 0)
p2 <- ggplot(a, aes(x = factor(group), y = distance, fill = factor(group))) + 
    geom_violin() + ylim(c(-1e4, 1e4)) +
    labs(title = "Mutation To TSS", x = "", y = "distance") +
    annotate("text", label = labels, x = unique(a$group), y = 0)
png(filename = paste(outpre, "violin.png", sep = "."), width = 2000, height = 2000)
multiplot(p, p1, p2, p3, cols = 1)
dev.off()
