
sort -k1,1 -k2,2n leafcutter_ds_effect_sizes_Tolarge.sig1.bed > leafcutter_ds_effect_sizes_Tolarge.sig1.bed2 
sort -k1,1 -k2,2n leafcutter_ds_effect_sizes_Tolarge.sig2.bed > leafcutter_ds_effect_sizes_Tolarge.sig2.bed2 
sort -k1,1 -k2,2n leafcutter_ds_effect_sizes_Tolarge.sig3.bed > leafcutter_ds_effect_sizes_Tolarge.sig3.bed2 
sort -k1,1 -k2,2n leafcutter_ds_effect_sizes_Tolarge.sig4.bed > leafcutter_ds_effect_sizes_Tolarge.sig4.bed2

mv leafcutter_ds_effect_sizes_Tolarge.sig1.bed2 leafcutter_ds_effect_sizes_Tolarge.sig1.bed
mv leafcutter_ds_effect_sizes_Tolarge.sig2.bed2 leafcutter_ds_effect_sizes_Tolarge.sig2.bed
mv leafcutter_ds_effect_sizes_Tolarge.sig3.bed2 leafcutter_ds_effect_sizes_Tolarge.sig3.bed
mv leafcutter_ds_effect_sizes_Tolarge.sig4.bed2 leafcutter_ds_effect_sizes_Tolarge.sig4.bed

m6A_peak="/usr/people/bioc1387/Project/ChrM6A-seq/Peaks/m6AIP_Cfg1_Cfg2.narrowPeak.refined.bed"

prefix="leafcutter_ds_effect_sizes_Tolarge.sig"

for file in 1 2 3 4
do
	file=$prefix$file
	echo $file.bed
	echo "number of lines"
	more $file.bed | wc -l
	echo "## Overlapped m6A Cfg1 Cfg2 peaks: "
	closestBed -b $m6A_peak -a $file.bed -d | awk '$13<=500 && $13>=0' | cut -f1-6 | sort | uniq |wc -l
	closestBed -b $m6A_peak -a $file.bed -d | awk '$13<=500 && $13>=0' | cut -f4 | sort | uniq | awk '{print $0"\t"}'> $file.Cluster
	closestBed -b $m6A_peak -a $file.bed -d | awk '$13>500||$13<0'     | cut -f1-6 | sort | uniq | awk '{print $0}'> $file.nom6A.bed
	closestBed -b $m6A_peak -a $file.bed -d | awk '$13<=500 && $13>=0' | cut -f1-6 | sort | uniq | awk '{print $0}'> $file.m6A.bed
	closestBed -b $m6A_peak -a $file.bed -d | awk '$13<=500 && $13>=0' | cut -f7-12 | sort | uniq > $file.m6APeaks.bed
done

