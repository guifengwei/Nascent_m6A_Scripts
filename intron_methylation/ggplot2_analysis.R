
library(ggplot2)
library(dplyr)

Data <- read.table("m6AIP_Cfg1_Cfg2.summits.Intersects_with_GENCODE_vM24_MaxORF_LongestNcRNA.intron.Table", header=T)

pdf("Pattern.pdf", height=4.5, width=6)
#ggplot(Data, aes(x=ExonORIntron, y=m6A_Intensity))+ geom_boxplot(notch=T) + facet_grid(. ~ Class)
ggplot(Data, aes(x=ExonORIntron, y=m6A_Intensity))+ geom_boxplot(notch=T) + facet_grid(. ~ Class) + coord_cartesian(ylim=c(-1, 6)) 

#########################
Data_Phastcon1 <- Data %>% select(Overlap_Phastcon, ExonORIntron, Class) %>% rename(Phastcon=Overlap_Phastcon)
Data_Phastcon2 <- Data %>% select(Control_Phastcon, ExonORIntron, Class) %>% rename(Phastcon=Control_Phastcon)
Data_Phastcon <- rbind(Data_Phastcon1, Data_Phastcon2)

Data_Phastcon$Region<- factor(rep(c("Peak", "Control"), each=dim(Data_Phastcon1)[1]), levels=c("Peak", "Control"))
ggplot(Data_Phastcon, aes(x=Region, y=Phastcon) ) + geom_boxplot(notch=T) + facet_grid(. ~ Class) 

########################
Data_GC1 <- Data %>% select(Overlap_GC, ExonORIntron, Class) %>% rename(GC=Overlap_GC)
Data_GC2 <- Data %>% select(Control_GC, ExonORIntron, Class) %>% rename(GC=Control_GC)
Data_GC <- rbind(Data_GC1, Data_GC2)

Data_GC$Region<- factor(rep(c("Peak", "Control"), each=dim(Data_GC1)[1]), levels=c("Peak", "Control"))
#ggplot(Data_GC, aes(x=Region, y=GC) ) + geom_boxplot(notch=T) + facet_grid(ExonORIntron ~ Class) 
ggplot(Data_GC, aes(x=Region, y=GC) ) + geom_boxplot(notch=T) + facet_grid(. ~ Class) 

#########################
Data_L1 <- Data %>% select(Intron_length, ExonORIntron, Class) %>% rename(Length=Intron_length)
Data_L2 <- Data %>% select(Control_length, ExonORIntron, Class) %>% rename(Length=Control_length)
Data_Length <- rbind(Data_L1, Data_L2)

Data_Length$Region<- factor(rep(c("Peak", "Control"), each=dim(Data_L1)[1]), levels=c("Peak", "Control"))
#ggplot(Data_Length, aes(x=Region, y=log2(Length)) ) + geom_boxplot(notch=T) + facet_grid(. ~ Class)  
ggplot(Data_Length, aes(x=Region, y=log2(Length)) ) + geom_boxplot(notch=T) + facet_grid(. ~ Class) + coord_cartesian(ylim=c(6,22.5))
dev.off()

pdf("rPos_m6AIntensity_Phastcon.pdf", width=12, height=10)
#ggplot(Data, aes(x=rPos, y=m6A_Intensity, shape=ExonORIntron)) + geom_point(aes(size=Overlap_Phastcon)) + scale_shape_manual(values=c(19,1)) + facet_grid(Class ~ .) 
ggplot(Data, aes(x=rPos, y=m6A_Intensity, shape=ExonORIntron)) + geom_point(aes(size=Overlap_Phastcon, col=ExonORIntron), alpha=0.75) + scale_shape_manual(values=c(19,19)) + scale_colour_manual(values=c("#F8766D", "slateblue4")) + facet_grid(Class ~ .) 
