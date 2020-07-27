


namefile="m6AIP_Cfg1_Cfg2.summits.Intersects_with_GENCODE_vM24_MaxORF_LongestNcRNA.intron"

python FindingMatchControl.py $namefile > $namefile.Overlap_and_matchControl2

genome="/usr/people/bioc1387/Project/mm10/Sequences/WholeGenome/mm10.fa"
Phastcon="/usr/people/bioc1387/Data/mm10_PhastCons/mm10.phastCon60way.bw"


### for the Phastcon analysis and GC analysis

## Overlap.bed
more $namefile.Overlap_and_matchControl2 | awk '{FS=OFS="\t"}{print $1,$14,$15,$4"|"$10,0,$6}'  > $namefile.Overlap.bed

## Control.bed
more $namefile.Overlap_and_matchControl2 | awk '{FS=OFS="\t"}{print $1,$16,$17,$4"|"$10,0,$6}'  > $namefile.matchControl.bed



for bed in $namefile.Overlap.bed $namefile.matchControl.bed
do
    /usr/people/bioc1387/Project/UCSCTools/bigWigAverageOverBed $Phastcon $bed $bed.Phastcon.tab -bedOut=$bed.Phastcon.bed
    fastaFromBed -fi $genome -bed $bed -name -s -fo $bed.fa 
    python /usr/people/bioc1387/Scripts/My_Scripts/xCalculateGC.py $bed.fa > $bed.GC
done


#python Aggregate.ExonIntron_anno.Phastcon.GC.m6A_intensity.py
