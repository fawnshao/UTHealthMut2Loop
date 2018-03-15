library(data.table)
library(ggplot2)
args <- commandArgs(TRUE)
input <- fread(args[1], sep = "\t", header = T)
count <- table(input[,2])
class <- length(count)

for (i in 3:ncol(input)){
	datax <- input[,c(2,i),with = FALSE]
	colnames(datax) <- c("type", "value")
	a <- table(datax[datax$value > 0,1])
	if(length(a) < class){
		a <- count - table(datax[datax$value == 0,1])
	}
	b <- a/count
	datay <- as.data.frame(b)
	write.table(data.frame(colnames(input)[i],t(datay[,2])), 
		file = paste(args[1],"stats.tsv", sep = "."), 
		append = TRUE, quote = TRUE, sep = "\t",
		row.names = FALSE, col.names = FALSE)
}

