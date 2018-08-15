library(data.table)
library(ggplot2)
args <- commandArgs(TRUE) 
# args <-  c("SRR5831489.q0.01.PET.dis")
datax <- fread(args[1], sep = "\t", header = F)
distancebin <- ceiling(datax$V4/1000)
datax[, disbin := distancebin ]
bins <- unique(distancebin)
med <- bins
sdv <- bins
for(i in 1:length(bins)){
    med[i] <- median(datax[disbin==bins[i],V3])
    sdv[i] <- sd(datax[disbin==bins[i],V3])
}
datax[, medianofdisbin := sapply(distancebin, function(x){med[bins==x]}) ]
datax[, sdofdisbin := sapply(distancebin, function(x){sdv[bins==x]}) ]
colnames(datax)[1:4] <-  c("node1", "node2", "ObservedPairs", "distance")
datax[, NormPairs := ObservedPairs / medianofdisbin]
write.table(datax, file = paste(args[1], "normpairs.tsv", sep = "."), 
	sep = "\t", quote = F, row.names = F)