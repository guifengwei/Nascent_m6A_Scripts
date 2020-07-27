


data_A5SSA3SS <- read.table("Splicing_A5SSA3SS", header=F)
data_ES <- read.table("Splicing_ES", header=F)
data_pfIR <- read.table("Splicing_pfIR", header=F)

Data <- data.frame("Delta_Psi"= c(data_A5SSA3SS$V4, data_ES$V4, data_pfIR$V7),
                  "m6A_type" = factor(c(as.vector(data_A5SSA3SS$V2), as.vector(data_ES$V2), as.vector(data_pfIR$V2)), levels=c("m6A", "nom6A")),
"Splicing_type" = factor(c(rep("ASS", dim(data_A5SSA3SS)[1]),  rep("ES", dim(data_ES)[1]), rep('pfIR', dim(data_pfIR)[1]) ), levels=c("ASS", "ES", "pfIR"))
			)

library(ggplot2)

pdf("Splicing_plot.pdf", width=7, height=4)
ggplot(Data, aes(y=Delta_Psi, x=m6A_type)) + facet_grid(. ~ Splicing_type) + geom_boxplot(outlier.shape = 3) + geom_point(position=position_jitter(width=0.2, height=0.001),shape=16, colour="purple", alpha=.55, size=2) + coord_cartesian(ylim=c(-0.3, 0.45))
dev.off()
