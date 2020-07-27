
#### positive strand is 1 and negative strand is 2
m6AIP16_p="m6A_IP16_STAR/m6A_IP16_p.Norm.bw"
m6AIP16_n="m6A_IP16_STAR/m6A_IP16_n.Norm.bw"
m6AIP17_p="m6A_IP17_STAR/m6A_IP17_p.Norm.bw"
m6AIP17_n="m6A_IP17_STAR/m6A_IP17_n.Norm.bw"

m6AInput16_p="m6A_Input16_STAR/m6A_Input16_p.Norm.bw"
m6AInput16_n="m6A_Input16_STAR/m6A_Input16_n.Norm.bw"
m6AInput17_p="m6A_Input17_STAR/m6A_Input17_p.Norm.bw"
m6AInput17_n="m6A_Input17_STAR/m6A_Input17_n.Norm.bw"

#m6AIP16_p="Abcam_chrRNA_m6A-seq/17_E14_chrRNA12_m6AIP_STAR/17_E14_chrRNA12_m6AIP_p.Norm.bw"
#m6AIP16_n="Abcam_chrRNA_m6A-seq/17_E14_chrRNA12_m6AIP_STAR/17_E14_chrRNA12_m6AIP_n.Norm.bw"
#m6AIP17_p="Abcam_chrRNA_m6A-seq/18_E14_chrRNA13_m6AIP_STAR/18_E14_chrRNA13_m6AIP_p.Norm.bw"
#m6AIP17_n="Abcam_chrRNA_m6A-seq/18_E14_chrRNA13_m6AIP_STAR/18_E14_chrRNA13_m6AIP_n.Norm.bw"
#
#m6AInput16_p="Abcam_chrRNA_m6A-seq/15_E14_chrRNA12_input_STAR/15_E14_chrRNA12_input_p.Norm.bw"
#m6AInput16_n="Abcam_chrRNA_m6A-seq/15_E14_chrRNA12_input_STAR/15_E14_chrRNA12_input_n.Norm.bw"
#m6AInput17_p="Abcam_chrRNA_m6A-seq/16_E14_chrRNA13_input_STAR/16_E14_chrRNA13_input_p.Norm.bw"
#m6AInput17_n="Abcam_chrRNA_m6A-seq/16_E14_chrRNA13_input_STAR/16_E14_chrRNA13_input_n.Norm.bw"

#m6A_peak="/usr/people/bioc1387/Project/METTL3dTag/METTL3dTag_m6A/m6A_DataSet/mESC_Geula_m6A_peak.m6A.bed"
m6A_peak=$1

name="$(echo $m6A_peak | awk -F "/" '{print $NF}')"
echo $name

for file in $m6AIP16_p $m6AIP16_n $m6AIP17_p $m6AIP17_n $m6AInput16_p $m6AInput16_n $m6AInput17_p $m6AInput17_n
do
	echo $file
	~/Project/UCSCTools/bigWigAverageOverBed $file $m6A_peak ${file}_m6A.txt -bedOut=${file}_m6A.bed
	rm ${file}_m6A.txt
done

####
paste ${m6AIP16_p}_m6A.bed ${m6AIP16_n}_m6A.bed ${m6AIP17_p}_m6A.bed ${m6AIP17_n}_m6A.bed | cut -f1-7,14,21,28 > ${m6A_peak}_IP.bed

paste ${m6AInput16_p}_m6A.bed ${m6AInput16_n}_m6A.bed ${m6AInput17_p}_m6A.bed ${m6AInput17_n}_m6A.bed | cut -f1-7,14,21,28 > ${m6A_peak}_input.bed

paste ${m6A_peak}_IP.bed ${m6A_peak}_input.bed | cut -f1-10,17-20 > ${m6A_peak}_output.bed

###
rm ${m6AIP16_p}_m6A.bed ${m6AIP16_n}_m6A.bed ${m6AIP17_p}_m6A.bed ${m6AIP17_n}_m6A.bed
rm  ${m6AInput16_p}_m6A.bed ${m6AInput16_n}_m6A.bed ${m6AInput17_p}_m6A.bed ${m6AInput17_n}_m6A.bed

rm ${m6A_peak}_IP.bed ${m6A_peak}_input.bed

#more ${m6A_peak}_output.bed | awk '{FS=OFS="\t"}{if($6=="+"){print $1,$2,$3,$4,$5,$6,$7,$9,$11,$13,$15,$17,$19,$21}else{print $1,$2,$3,$4,$5,$6,$8,$10,$12,$14,$16,$18,$20,$22}}' > $name.result

#rm ${m6A_peak}_output.bed

