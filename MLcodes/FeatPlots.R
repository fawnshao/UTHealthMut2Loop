library(data.table)
library(ggplot2)
args <- c("hkg.tsg.srtbyPCA.class.sequenceFeatures", 
	"hkg.tsg.srtbyPCA.class.Homer", 
	"hkg.tsg.srtbyPCA.class.GTRD",
	"hkg.tsg.srtbyPCA.class.roadmap.DNase",
	"hkg.tsg.srtbyPCA.class.roadmap.histone",
	"hkg.tsg.srtbyPCA.class.roadmap.meth",
	"hkg.tsg.srtbyPCA.class.HiC",
	"hkg.tsg.srtbyPCA.class.HiChIP",
	"hkg.tsg.srtbyPCA.class.PCHiC",
	"hkg.tsg.srtbyPCA.class.PCHiC.oe",
	"hkg.tsg.srtbyPCA.class.enhanceratlas"
	)

simFeatBox <- function(x = input, pre = "hkg.tsg.srtbyPCA"){
	require(ggplot2)
	x[1:2551,2] <- "HKG"
	x[2552:5286,2] <- "singleTSG"
	for (i in 3:ncol(x)){
		png(filename = paste(pre, colnames(x)[i], "sim.boxplot.png", sep = "."), width = 600, height = 600)
		datax <- x[,c(2,i),with = FALSE]
		colnames(datax) <- c("type", "value")
		ymax <- quantile(datax$value, probs = 0.95)
		myplot <- ggplot(data = datax, aes(x = type, y = value, fill = type)) + 
			geom_boxplot(position = position_dodge(1), outlier.alpha = 0.1, outlier.size = 0.1) + 
			scale_y_continuous(limits = c(0, ymax)) +
			ggtitle(colnames(x)[i]) + 
			theme(legend.position = "none")
		print(myplot)
		dev.off()
	}
}

simFeatBar <- function(x = input, pre = "hkg.tsg.srtbyPCA"){
	require(ggplot2)
	x[1:2551,2] <- "HKG"
	x[2552:5286,2] <- "singleTSG"
	count <- table(x[,2])
	class <- length(count)
	for (i in 3:ncol(x)){
		datax <- x[,c(2,i),with = FALSE]
		colnames(datax) <- c("type", "value")
		a <- table(datax[datax$value > 0,1])
		if(length(a) < class){
			a <- count - table(datax[datax$value == 0,1])
		}
		b <- a/count
		if(max(b) > 0.5 && min(b) < 0.5){
			datay <- as.data.frame(b)
			png(filename = paste(pre, colnames(x)[i], "sim.barplot.png", sep = "."), width = 600, height = 600)
			myplot <- ggplot(data = datay, aes(x = Var1, y = Freq, fill = Var1)) + 
				geom_bar(stat = "identity") +
				ggtitle(colnames(x)[i]) + 
				theme(legend.position = "none")
			print(myplot)
			dev.off()
		}
	}
}

fullFeatBox <- function(x = input, pre = "hkg.tsg.srtbyPCA"){
	require(ggplot2)
	for (i in 3:ncol(x)){
		png(filename = paste(pre, colnames(x)[i], "boxplot.png", sep = "."), width = 1200, height = 600)
		datax <- x[,c(2,i),with = FALSE]
		colnames(datax) <- c("type", "value")
		ymax <- quantile(datax$value, probs = 0.95)
		myplot <- ggplot(data = datax, aes(x = type, y = value, fill = type)) + 
			geom_boxplot(position = position_dodge(1), outlier.alpha = 0.1, outlier.size = 0.1) + 
			scale_y_continuous(limits = c(0, ymax)) +
			ggtitle(colnames(x)[i]) + 
			theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1))
		print(myplot)
		dev.off()
	}
}

fullFeatBar <- function(x = input, pre = "hkg.tsg.srtbyPCA"){
	require(ggplot2)
	count <- table(x[,2])
	class <- length(count)
	for (i in 3:ncol(x)){
		datax <- x[,c(2,i),with = FALSE]
		colnames(datax) <- c("type", "value")
		aa <- names(count)
		a <- count
		for(j in 1:length(aa)){
			a[j] <- length(na.omit(match(as.matrix(datax[datax$value > 0,])[,1], aa[j])))
		}
		# a <- table(datax[datax$value > 0,1])
		# if(length(a) < class){
		# 	a <- count - table(datax[datax$value == 0,1])
		# }
		b <- a/count
		# if(max(b) > 0){
		if(max(b) > 0.5 && min(b) < 0.5){
			datay <- as.data.frame(b)
			png(filename = paste(pre, colnames(x)[i], "barplot.png", sep = "."), width = 1200, height = 600)
			myplot <- ggplot(data = datay, aes(x = Var1, y = Freq, fill = Var1)) + 
				geom_bar(stat = "identity") +
				ggtitle(colnames(x)[i]) + 
				theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1))
			print(myplot)
			dev.off()
		}
	}
}

simFeatPerc <- function(x = input, pre = "hkg.tsg.srtbyPCA"){
	require(ggplot2)
	x[1:2551,2] <- "HKG"
	x[2552:5286,2] <- "singleTSG"
	count <- table(x[,2])
	class <- length(count)
	output <- data.frame()
	for (i in 3:ncol(x)){
		datax <- x[,c(2,i),with = FALSE]
		colnames(datax) <- c("type", "value")
		a <- table(datax[datax$value > 0,1])
		if(length(a) < class){
			a <- count - table(datax[datax$value == 0,1])
		}
		b <- a/count
		output <- rbind.data.frame(output, c(a,b))
	}
	colnames(output) <- c(paste("count",names(b)), paste("ratio",names(b)))
	rownames(output) <- colnames(x)[3:ncol(x)]
	write.table(data.frame(rownames(output), output), 
		file = paste(pre, "percentage.tsv", sep = "."), 
		append = T,
		row.names = F, sep = "\t")
}


input <- fread(args[1], sep = "\t", header = T)
# simFeatBox(input)
# simFeatBar(input)
# fullFeatBar(input)
# fullFeatBox(input)
simFeatPerc(input)

input <- fread(args[2], sep = "\t", header = T)
# colnames(input) <- gsub("/",".",colnames(input))
# colnames(input) <- gsub("%",".",colnames(input))
# simFeatBar(input)
# fullFeatBar(input)
simFeatPerc(input)


input <- fread(args[3], sep = "\t", header = T)
# colnames(input) <- gsub("/",".",colnames(input))
# colnames(input) <- gsub("%",".",colnames(input))
# simFeatBar(input)
# fullFeatBar(input)
simFeatPerc(input)

input <- fread(args[4], sep = "\t", header = T)
# colnames(input) <- gsub("/",".",colnames(input))
# colnames(input) <- gsub("%",".",colnames(input))
# simFeatBar(input)
# fullFeatBar(input)
simFeatPerc(input)

input <- fread(args[5], sep = "\t", header = T)
# colnames(input) <- gsub("/",".",colnames(input))
# colnames(input) <- gsub("%",".",colnames(input))
# simFeatBar(input)
# fullFeatBar(input)
simFeatPerc(input)

input <- fread(args[6], sep = "\t", header = T)
# colnames(input) <- gsub("/",".",colnames(input))
# colnames(input) <- gsub("%",".",colnames(input))
# simFeatBox(input)
simFeatPerc(input)

input <- fread(args[7], sep = "\t", header = T)
# colnames(input) <- gsub("/",".",colnames(input))
# colnames(input) <- gsub("%",".",colnames(input))
# simFeatBar(input)
# fullFeatBar(input)
simFeatPerc(input)

input <- fread(args[8], sep = "\t", header = T)
# colnames(input) <- gsub("/",".",colnames(input))
# colnames(input) <- gsub("%",".",colnames(input))
# simFeatBar(input)
# fullFeatBar(input)
simFeatPerc(input)

input <- fread(args[9], sep = "\t", header = T)
# colnames(input) <- gsub("/",".",colnames(input))
# colnames(input) <- gsub("%",".",colnames(input))
# simFeatBar(input)
# fullFeatBar(input)
simFeatPerc(input)

input <- fread(args[10], sep = "\t", header = T)
# colnames(input) <- gsub("/",".",colnames(input))
# colnames(input) <- gsub("%",".",colnames(input))
# simFeatBar(input)
# fullFeatBar(input)
simFeatPerc(input)

input <- fread(args[11], sep = "\t", header = T)
# colnames(input) <- gsub("/",".",colnames(input))
# colnames(input) <- gsub("%",".",colnames(input))
# simFeatBox(input)
# fullFeatBox(input)
# simFeatBar(input)
# fullFeatBar(input)
colnames(input) <- gsub("-","_",colnames(input))
simFeatPerc(input)
