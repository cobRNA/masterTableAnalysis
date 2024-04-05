#!/usr/bin/env bash

while read catalogue || [ -n "$catalogue" ]  # solution to no terminal newline problem
do  
    echo $catalogue
done < ./Data/Catalogues/catalogues_list








# Convert catalogue to gff
cat cls.hg38.bed12 | ../../../Utils/bed12togff > cls.hg38.bed12.gff

# Sort catalogue
cat cls.hg38.bed12.gff | sort -k1,1 -k12,12 -k4,4n -k5,5n "$@" > cls.hg38_sorted.gff

# Create clusters
# cat cls.5p.ends | sort -k1,1 -k2,2n | bedtools merge -c 4 -o collapse -s -d 5 -i - > cls_5p_clusters.bed
# cat cls.3p.ends | sort -k1,1 -k2,2n | bedtools merge -c 4 -o collapse -s -d 5 -i - > cls_3p_clusters.bed
## With proper BED6 output formatting
#cat cls.5p.ends | sort -k1,1 -k2,2n | bedtools merge -c 4,5,6 -o collapse,distinct,distinct -s -d 5 -i - > cls_5p_clusters.bed
#cat cls.3p.ends | sort -k1,1 -k2,2n | bedtools merge -c 4,5,6 -o collapse,distinct,distinct -s -d 5 -i - > cls_3p_clusters.bed
## Using Julien magic from article
#cat cls.5p.ends | ../../../Utils/sortbed | bedtools merge -c 4,5,6 -o collapse,count,distinct -s -d 5 -i stdin | ../../../Utils/sortbed > cls_5p_clusters.bed
#cat cls.3p.ends | ../../../Utils/sortbed | bedtools merge -c 4,5,6 -o collapse,count,distinct -s -d 5 -i stdin | perl -F"\t" -lane 'if($F[5] eq "+"){$F[1]=$F[2]-1}elsif($F[5] eq "-"){$F[2]=$F[1]+1}else{die} print join("\t",@F);' | ../../../Utils/sortbed > cls_3p_clusters.bed
# | awk '$5>1'
## Using only polyA and CAGE supported ends - stupid idea
#cat cls.5.bed.vsCage.fantom.bedtsv |../../../Utils/sortbed | bedtools merge -c 4,5,6 -o collapse,count,distinct -s -d 5 -i stdin | ../../../Utils/sortbed > cls_5p_clusters.bed
#cat cls.3.bed.vspolyAsignals.bedtsv | ../../../Utils/sortbed | bedtools merge -c 4,5,6 -o collapse,count,distinct -s -d 5 -i stdin | ../../../Utils/sortbed > cls_3p_clusters.bed
## Using 0 distance between ends (the default for bedtools merge) - still fails
#cat cls.5p.ends | sort -k1,1 -k2,2n | bedtools merge -c 4,5,6 -o collapse,distinct,distinct -s -d 0 -i - | ../../../Utils/sortbed > cls_5p_clusters.bed
#cat cls.3p.ends | sort -k1,1 -k2,2n | bedtools merge -c 4,5,6 -o collapse,distinct,distinct -s -d 0 -i - | ../../../Utils/sortbed > cls_3p_clusters.bed
## Using 3p and 5p ends as clusters - not merging anything at all
cat cls.5p.ends | ../../../Utils/sortbed > cls_5p_clusters.bed
cat cls.3p.ends | ../../../Utils/sortbed > cls_3p_clusters.bed

# Sort transcript files
cat cls.tx.cage | sort -k1,1 > cls.tx.cage_sorted
cat cls.tx.polyA | sort -k1,1 > cls.tx.polyA_sorted


# Run anchoring
../../../Utils/anchorTranscriptsEnds.pl cls.hg38_sorted.gff cls.tx.cage_sorted cls.tx.polyA_sorted cls_5p_clusters.bed cls_3p_clusters.bed > anchored_cls.gff                     

# Look for erroneously modified exones (corrupted exones i.e. $4>$5 were removed to run bedtools)
cat anchored_cls.gff | awk '$4<=$5' | bedtools intersect -a - -b cls.hg38_sorted.gff -s -f 1.00 -r -loj > loj_anchored
cat anchored_cls.gff | awk '$4<=$5' | bedtools intersect -a cls.hg38_sorted.gff -b - -s -f 1.00 -r -loj > loj_cls.hg38

# Check resulting file
echo '##############################'
echo 'Looking for corrupted entries...'
cat anchored_cls.gff | awk '$4>$5'
echo '##############################'
echo 'Cheking ENST00000583400.3 transcript...'
cat anchored_cls.gff | grep ENST00000583400.3
echo '##############################'
echo 'Resulting file preview:'
cat anchored_cls.gff | head -n 10
echo '##############################'
echo 'Looking for Erroneously modified exones in loj_anchored file...'
echo '------------------------------'
cat loj_anchored | grep '\-1' | awk -F';'  '$3!=" fakeExon \"yes\""' | head
printf 'Erroneously modified exons: '
cat loj_anchored | grep '\-1' | awk -F';'  '$3!=" fakeExon \"yes\""' | wc -l
echo '##############################'
echo 'Looking for erroneously modified exones using loj_cls.hg38 file...'
echo '------------------------------'
cat loj_cls.hg38 | grep '\-1' | head
printf 'Erroneously modified exons: '
cat loj_cls.hg38 | grep '\-1' | wc -l
