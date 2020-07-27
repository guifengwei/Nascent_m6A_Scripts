	
splitstring <- function(data_vector, sep="|", k = 1){
	x_vector <- c()
	n=length(data_vector)
	for(i in 1:n){
	    x_vector = append(x_vector,  unlist(strsplit(toString(data_vector[i]), "[|]"))[k],  after=length(x_vector))
	    }
	return(x_vector)
}


Cfg1 <- read.table("../m6AIP_ConfGroup1.narrowPeak.refined.bed", header=F)
Cfg2 <- read.table("../m6AIP_ConfGroup2.narrowPeak.refined.bed_abcam", header=F)

Cfg1_peak = Cfg1$V4
Cfg2_abcam_peak = Cfg2[grep("abcam_", Cfg2$V4),]$V4
Cfg2_sysy_peak = Cfg2[grep("SySy_", Cfg2$V4),]$V4


data_GC_overlap <- read.table("m6AIP_Cfg1_Cfg2.summits.Intersects_with_GENCODE_vM24_MaxORF_LongestNcRNA.intron.Overlap.GC", header=F)
data_GC_control <- read.table("m6AIP_Cfg1_Cfg2.summits.Intersects_with_GENCODE_vM24_MaxORF_LongestNcRNA.intron.matchControl.GC", header=F)
data_GC_overlap$name = splitstring(data_GC_overlap$V1)
data_GC_control$name = splitstring(data_GC_control$V1)

data_Phastcon_overlap <- read.table("m6AIP_Cfg1_Cfg2.summits.Intersects_with_GENCODE_vM24_MaxORF_LongestNcRNA.intron.Overlap.bed.Phastcon.bed", header=F)
data_Phastcon_control <- read.table("m6AIP_Cfg1_Cfg2.summits.Intersects_with_GENCODE_vM24_MaxORF_LongestNcRNA.intron.matchControl.bed.Phastcon.bed", header=F)
data_Phastcon_overlap$name = splitstring(data_Phastcon_overlap$V4)
data_Phastcon_control$name = splitstring(data_Phastcon_control$V4)

data_intron_length <- read.table("m6AIP_Cfg1_Cfg2.summits.Intersects_with_GENCODE_vM24_MaxORF_LongestNcRNA.intron.Overlap_and_matchControl2", header=F)


### Phastcon
boxplot(data_Phastcon_overlap[data_Phastcon_overlap$name %in% Cfg1_peak,]$V7, data_Phastcon_control[data_Phastcon_control$name %in% Cfg1_peak,]$V7,                             data_Phastcon_overlap[data_Phastcon_overlap$name %in% Cfg2_sysy_peak,]$V7, data_Phastcon_control[data_Phastcon_control$name %in% Cfg2_sysy_peak,]$V7,                   data_Phastcon_overlap[data_Phastcon_overlap$name %in% Cfg2_abcam_peak,]$V7, data_Phastcon_control[data_Phastcon_control$name %in% Cfg2_abcam_peak,]$V7, notch=T)


boxplot(data_GC_overlap[data_GC_overlap$name %in% Cfg1_peak,]$V2, data_GC_control[data_GC_control$name %in% Cfg1_peak,]$V2,                                                   data_GC_overlap[data_GC_overlap$name %in% Cfg2_sysy_peak,]$V2, data_GC_control[data_GC_control$name %in% Cfg2_sysy_peak,]$V2,                                           data_GC_overlap[data_GC_overlap$name %in% Cfg2_abcam_peak,]$V2, data_GC_control[data_GC_control$name %in% Cfg2_abcam_peak,]$V2, notch=T)

boxplot(log2(data_intron_length[data_intron_length$V4 %in% Cfg1_peak,]$V18), log2(data_intron_length[data_intron_length$V4 %in% Cfg1_peak,]$V19),                               log2(data_intron_length[data_intron_length$V4 %in% Cfg2_sysy_peak,]$V18),log2(data_intron_length[data_intron_length$V4 %in% Cfg2_sysy_peak,]$V19),                     log2(data_intron_length[data_intron_length$V4 %in% Cfg2_abcam_peak,]$V18),log2(data_intron_length[data_intron_length$V4 %in% Cfg2_abcam_peak,]$V19), notch=TRUE)


