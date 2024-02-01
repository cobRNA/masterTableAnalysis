###############
# Merging TMs from RNA catalogues using anchored merging.
###############

# -------------------
# extract ends
# -------------------

for file in `ls ./Data/Catalogues/*.hg38.bed12`
do
    echo $file
    lab=`basename $file| awk -F "." '{print $1}'`
    while read end dist
    do
        echo $end
        cat $file | ./Utils/extractTranscriptEndsFromBed12.pl $end | ./Utils/sortbed | ./Utils/bedtools merge -s -d $dist -c 4 -o collapse -i stdin | awk '{print $1"\t"$2"\t"$3"\t"$5"\t0\t"$4}' > ./Data/Processed/$lab.$end.bed
    done < ./Data/Others/ends.dist.tsv
done


# -------------------
# polyA
# -------------------

for file in `ls ./Data/Catalogues/*.hg38.bed12`
do
    echo $file
    lab=`basename $file| awk -F "." '{print $1}'`
    while read end dist
    do
        echo $end
        cat ./Data/Processed/$lab.$end.bed | ./Utils/sortbed | ./Utils/bedtools slop -s -l 5 -r -5 -i stdin -g ./Data/Others/chromInfo_hg38.txt | ./Utils/bedtools intersect -u -s -a stdin -b ./Data/Others/hg38.polyAsignals.bed > ./Data/Processed/$lab.$end.bed.vspolyAsignals.bedtsv
    done < ./Data/Others/ends.dist.tsv
done




 

