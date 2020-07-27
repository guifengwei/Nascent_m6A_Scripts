
library(ggplot2)
library(dplyr)

splitstring <- function(data_vector, sep=":", k = 1){
	x_vector <- c()
	n=length(data_vector)
	for(i in 1:n){
		x_vector = append(x_vector,   unlist(strsplit(toString(data_vector[i]), sep))[k],  after=length(x_vector))
		}
	return(x_vector)
	}

data <- read.table("leafcutter_ds_cluster_significance.txt", header=T, sep="\t")
data <- data[!is.na(data$p),]

data_sig1 <- data[data$p.adjust<0.05,]
data_sig2 <- data[data$p.adjust<0.2 & data$p.adjust>=0.05,]
data_sig3 <- data[data$p <=0.05 & data$p.adjust>=0.2, ]
data_sig4 <- data[data$p >0.05, ]

plot(1:dim(data)[1], -log10(sort(data$p.adjust)), ylim=c(5,20), xlim=c(1,8))
plot(1:dim(data)[1], -log10(sort(data$p.adjust)), ylim=c(0,5))
abline(h=c(-log10(0.05), -log10(0.2), -log10(0.72889)), lty=2, col="red")

LargeBed <- read.table("leafcutter_ds_effect_sizes_Tolarge.bed", header=F)

LargeBed_sig1 <- LargeBed[LargeBed$V4 %in% splitstring(data_sig1$cluster, sep=":", k=2),]
LargeBed_sig2 <- LargeBed[LargeBed$V4 %in% splitstring(data_sig2$cluster, sep=":", k=2),]
LargeBed_sig3 <- LargeBed[LargeBed$V4 %in% splitstring(data_sig3$cluster, sep=":", k=2),]
LargeBed_sig4 <- LargeBed[LargeBed$V4 %in% splitstring(data_sig4$cluster, sep=":", k=2),]

write.table(LargeBed_sig1, file="leafcutter_ds_effect_sizes_Tolarge.sig1.bed", sep="\t", quote=F, col.names=F, row.names=F)
write.table(LargeBed_sig2, file="leafcutter_ds_effect_sizes_Tolarge.sig2.bed", sep="\t", quote=F, col.names=F, row.names=F)
write.table(LargeBed_sig3, file="leafcutter_ds_effect_sizes_Tolarge.sig3.bed", sep="\t", quote=F, col.names=F, row.names=F)
write.table(LargeBed_sig4, file="leafcutter_ds_effect_sizes_Tolarge.sig4.bed", sep="\t", quote=F, col.names=F, row.names=F)
