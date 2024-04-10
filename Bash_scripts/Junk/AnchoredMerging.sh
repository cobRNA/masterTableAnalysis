#!/usr/bin/env bash


#############################################
# Clean-up generated data
#############################################

rm ./Data/.temp/*
rm ./Data/Processed/Catalogues_gff/*
rm ./Data/Processed/Anchored_catalogues/*
rm ./Data/Processed/Merged_anchored_catalogues/*
# Clean-up report file
> ./Data/Reports/anchoring


#############################################
# Prepare FantomCat catalogue
#############################################

echo '### Preparing FantomCat catalogue...'
echo '-----------------------------------------'
# Filter monoexonic tx and patches from FantomCat catalogue
cat ./Data/Catalogues/fantomCat.hg38.bed12 | grep -v -E 'chr7_KI270803v1_alt|chr14_GL000009v2_random|chr4_GL000008v2_random' | awk '$10 > 1' > ./Data/Catalogues/fantomCat.hg38_polyexonic.bed12
echo '############# DONE ####################'

#############################################
# Prepare catalogues
#############################################

echo '### Concatenating catalogues...'
echo '-----------------------------------------'
# Concatenate catalogues
while read catalogue || [ -n "$catalogue" ]  # solution to no terminal newline problem
do  
    # Extract base catalogue name
    cat_name=$(echo $catalogue | awk -F'.' '{print $1}')

    # Report progress
    echo $cat_name

    # # Convert catalogue to required gff format
    # cat ./Data/Catalogues/$catalogue.bed12 | ./Utils/bed12togff > ./Data/Processed/Catalogues_gff/$catalogue.gff

    # Concatenate all catalogues so they can be anchored all at once at later step
    cat ./Data/Catalogues/$catalogue.bed12 >> ./Data/.temp/ALL_catalogues.bed12

done < ./Data/Catalogues/catalogues_list

# Remove duplicated entries
duplicated=$(cat ./Data/.temp/ALL_catalogues.gff | wc -l)
cat ./Data/.temp/ALL_catalogues.gff | sort | uniq | ./Utils/sortgff > ./Data/Processed/Catalogues_gff/ALL_catalogues.gff
((duplicated=duplicated-$(cat ./Data/Processed/Catalogues_gff/ALL_catalogues.gff | wc -l)))
echo "Found and removed $duplicated entries."

echo '################## DONE ##################'

#############################################
# Extract concatenated transcripts ends
#############################################






#############################################
# Anchor full-length transcript models
#############################################

echo '#####################'
echo '### Anchoring concatenated catalogues...'
echo '#####################'

## Concatenate tx.cage for all catalogues
rm ./Data/Processed/Cage/ALL.tx.cage
cat ./Data/Processed/Cage/*.tx.cage | sort | uniq > ./Data/Processed/Cage/ALL.tx.cage

## Concatenate tx.polyA for all catalogues
rm ./Data/Processed/PolyA/ALL.tx.polyA
cat ./Data/Processed/PolyA/*.tx.polyA | sort | uniq > ./Data/Processed/PolyA/ALL.tx.polyA

## Concatenate 5p ends
rm ./Data/Processed/Extracted_ends/ALL.5.bed
cat ./Data/Processed/Extracted_ends/*.5.bed | ./Utils/sortbed | uniq > ./Data/Processed/Extracted_ends/ALL.5.bed

## Concatenate 3p ends
rm ./Data/Processed/Extracted_ends/ALL.3.bed
cat ./Data/Processed/Extracted_ends/*.3.bed | ./Utils/sortbed | uniq > ./Data/Processed/Extracted_ends/ALL.3.bed

## Sort input gff
cat ./Data/Processed/Catalogues_gff/ALL_catalogues.gff | sort -k1,1 -k12,12 -k4,4n -k5,5n "$@" | uniq > ./Data/Processed/Catalogues_gff/sorted_ALL_catalogues.gff

## Anchor TM
./Utils/anchorTranscriptsEnds.pl ./Data/Processed/Catalogues_gff/sorted_ALL_catalogues.gff ./Data/Processed/Cage/ALL.tx.cage ./Data/Processed/PolyA/ALL.tx.polyA ./Data/Processed/Extracted_ends/ALL.5.bed ./Data/Processed/Extracted_ends/ALL.3.bed > ./Data/Processed/Anchored_catalogues/anchored_ALL.gff



#############################################
# Run tmerge
#############################################

# Change dir for cleaner script
cd ./Data/Processed/Anchored_catalogues/

# Run Tmerge for each catalogue
echo '#####################'
echo '### Merging catalogues...'
echo '#####################'

cat anchored_fantomCat.hg38_polyexonic.gff anchored_cls.hg38.gff anchored_bigtrans.hg38.gff anchored_clsFL.hg38.gff anchored_gen.hg38.gff anchored_mitrans.hg38.gff anchored_noncode.hg38.gff anchored_pcConf.hg38.gff anchored_refseq.hg38.gff | ../../../Utils/sortgff | ../../../Utils/tmerge - --tmPrefix anch > ../../../Data/Processed/Merged_anchored_catalogues/anchored_tmerged_catalogues.gtf

# Run Tmerge for concatenated catalogues
echo '#####################'
echo '### Merging concatenatd catalogues...'
echo '#####################'

cat anchored_ALL.gff | ../../../Utils/sortgff | ../../../Utils/tmerge - --tmPrefix anch > ../../../Data/Processed/Merged_anchored_catalogues/anchored_tmerged_ALL_catalogues.gtf

# Change dir to previous state
cd ../../..




# Remove fake-exons
#cat ./Data/Processed/Merged_anchored_catalogues/anchored_tmerged_catalogues.gtf | grep -Fv 'fa'
