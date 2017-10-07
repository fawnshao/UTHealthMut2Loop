args <- commandArgs(TRUE)
expr <- args[1]
list <- args[2]
fc <- as.numeric(args[3])

data <- read.table(expr, sep = "\t")
gene.name <- as.vector(data[-1, 1])
sample.name <- data[1, -1]
sample.name.sim <- apply(sample.name, 2, function(x){substr(x, start = 1, stop = 15)})
expr.value <- matrix(as.numeric(as.matrix(data[-1, -1])), byrow = F, nrow = length(gene.name))

toprocess <- read.table(list, header = T)
ifLoop <- rep(" ", nrow(toprocess))
for(i in 1:nrow(toprocess)){
	x <- unlist(strsplit(as.vector(toprocess[i,1]), ","))
	y <- unlist(strsplit(as.vector(toprocess[i,2]), ","))
	z <- unlist(strsplit(as.vector(toprocess[i,3]), ","))
	u <- unlist(strsplit(as.vector(toprocess[i,4]), ","))
	p <- as.vector(toprocess[i,9])
	un <- length(intersect(gene.name, u))
	if(un > 1){
		ctr <- expr.value[is.element(gene.name, u), !is.element(sample.name.sim, z)]
		mut <- matrix(expr.value[is.element(gene.name, u), is.element(sample.name.sim, y)], 
			byrow = F, nrow = un)
		ctr.mean <- apply(ctr, 1, mean)
		mut.mean <- apply(mut, 1, mean)
		v <- wilcox.test(ctr.mean, mut.mean)$p.value
		w <- mut.mean/ctr.mean
		t <- data.frame(w, mut.mean, ctr.mean)
		rownames(t) <- gene.name[is.element(gene.name, u)]
		colnames(t) <- c("Fold-Change", paste(y, v, sep = ": "), "Average-nonmut-Tumor-Sample")

		if(nrow(t[t[,1] > fc & t[,2] > 10, ]) > 0 && nrow(t[t[,1] < 1/fc & t[,3] > 10, ]) > 0){
			out <- rbind.data.frame(t[t[,1] > fc & t[,2] > 10, ], t[t[,1] < 1/fc & t[,3] > 10, ])
			flag <- 0
			for(j in 1:nrow(out)){
				if(length(grep(rownames(out)[j],p)) > 0){
					flag <- 1
				}
			}
			if(flag == 1){
				write.csv(t, file = paste(x, y, "csv", sep = "."))
				ifLoop[i] <- "LoopBroken"
			}
		}
	}
}
write.csv(data.frame(toprocess, ifLoop), file = paste(list, "labeled", "csv", sep = "."))
