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
        cat $file | ./Utils/extractTranscriptEndsFromBed12.pl $end | ./Utils/sortbed  > ./Data/Processed/Extracted_ends/$lab.$end.bed
    done < ./Data/Others/ends.dist.tsv
done

# bez tego:
# | ./Utils/bedtools merge -s -d $dist -c 4 -o collapse -i stdin
# | awk '{print $1"\t"$2"\t"$3"\t"$5"\t0\t"$4}'


# -------------------
# polyA
# -------------------

while read sampl
do
    name=`echo $sampl`
    echo $name
    while read end dist
    do
       cat ./Data/Processed/Extracted_ends/$name.$end.bed | ./Utils/sortbed | ./Utils/jlagarde/bin.bkp/bedtools2/bin/bedtools slop -s -l 5 -r 5 -i stdin -g ./Data/Others/chromInfo_hg38.txt | ./Utils/jlagarde/bin.bkp/bedtools2/bin/bedtools intersect -u -s -a stdin -b ./Data/Others/hg38.polyAsignals.bed > ./Data/Processed/PolyA/$name.$end.bed.vspolyAsignals.bedtsv
    done < ./Data/Others/ends.dist.tsv
done < ./Data/Others/samples.tsv



# -----------------------
# CAGE FANTOM pahse 1+2
# -----------------------

while read sampl
do
    name=`echo $sampl`
    echo $name
    while read end dist
    do
        cat ./Data/Processed/Extracted_ends/$name.$end.bed | ./Utils/sortbed | ./Utils/jlagarde/bin.bkp/bedtools2/bin/bedtools slop -s -l 50 -r 50 -i stdin -g ./Data/Others/chromInfo_hg38.txt | ./Utils/jlagarde/bin.bkp/bedtools2/bin/bedtools intersect -u -s -a stdin -b ./Data/Others/hg38_fair+new_CAGE_peaks_phase1and2.bed > ./Data/Processed/Cage/$name.$end.bed.vsCage.fantom.bedtsv
    done < ./Data/Others/ends.dist.tsv
done < ./Data/Others/samples.tsv

