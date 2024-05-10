#!/usr/bin/env bash


##########################################################################################
# Clean-up generated data
##########################################################################################

rm ./Data/.temp/*
# rm ./Data/Processed/Catalogues_gff/*
# rm ./Data/Processed/Anchored_catalogues/*
# rm ./Data/Processed/Merged_anchored_catalogues/*
rm ./Data/Processed/mT+Catalogues/*
# Clean-up report file
> ./Data/Reports/anchoring


##########################################################################################
# Prepare FantomCat catalogue
##########################################################################################

echo '### Preparing FantomCat catalogue...'
# Filter monoexonic tx and patches from FantomCat catalogue
cat ./Data/Catalogues/fantomCat.hg38.bed12 | grep -v -E 'chr7_KI270803v1_alt|chr14_GL000009v2_random|chr4_GL000008v2_random' | awk '$10 > 1' > ./Data/Catalogues/fantomCat.hg38_polyexonic.bed12
echo 'Monoexonic transcripts and patches sucessfully removed from FantomCat catalogue.'
echo '############# DONE ####################'

##########################################################################################
# Prepare catalogues
##########################################################################################

echo '### Concatenating catalogues with the spliced-mT lncRNA slice...'
echo '@ Concatenating catalogues...'
# Concatenate catalogues
while read catalogue || [ -n "$catalogue" ]  # solution to no terminal newline problem
do  
    # Extract base catalogue name
    cat_name=$(echo $catalogue | awk -F'.' '{print $1}')

    # Concatenate all catalogues so they can be anchored all at once at later step
    cat ./Data/Catalogues/$catalogue.bed12 >> ./Data/.temp/ALL_catalogues.bed12

    # Report progress
    echo "Added $cat_name to ./Data/.temp/ALL_catalogues.bed12 file"

done < ./Data/Catalogues/catalogues_list
echo '~~NOTE: skipped gen.hg38 catalogue.~~'
echo '++ Done'

# Concatenate further with lncRNA spliced_mT slice.
echo '@ Concatenating catalogues with the lncRNA spliced_mT slice...'
cat ./Data/Processed/mT/Hv3_splicedmasterTable_refined_lncRNA.gtf | ./Utils/gff2bed_full.pl - | cat - ./Data/.temp/ALL_catalogues.bed12 > ./Data/.temp/ALL_catalogues+mT.bed12
echo '++ Done'

# Remove duplicated entries ($7 and $8 column from bed12 file was substituted with $2 and $3 respectively to avoid duplicated entires)
echo '@ Removing duplicated entries from resulting ALL_catalogues+mT.bed12 file...'
duplicated=$(cat ./Data/.temp/ALL_catalogues+mT.bed12 | wc -l)
cat ./Data/.temp/ALL_catalogues+mT.bed12 | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$2"\t"$3"\t"$9"\t"$10"\t"$11"\t"$12}' | sort | uniq | ./Utils/sortbed > ./Data/Processed/mT+Catalogues/ALL_catalogues+mT.bed12
((duplicated=duplicated-$(cat ./Data/Processed/mT+Catalogues/ALL_catalogues+mT.bed12 | wc -l)))
echo "Found and removed $duplicated entries."
echo '++ Done'

# Converting to gff file format
echo '@ Converting from bed12 to gff...'
cat ./Data/Processed/mT+Catalogues/ALL_catalogues+mT.bed12 | ./Utils/bed12togff > ./Data/Processed/mT+Catalogues/ALL_catalogues+mT.gff
echo '++ Done'
echo '################## DONE ##################'

##########################################################################################
# Concatenate CAGE and polyA supported transcript ids
##########################################################################################
echo '### Concatenating CAGE and polyA supported transcript ids...'
echo '@ Concatenating transcript ids from catalogues...'
# Concatenate catalogues
while read catalogue || [ -n "$catalogue" ]  # solution to no terminal newline problem
do  
    # Extract base catalogue name
    cat_name=$(echo $catalogue | awk -F'.' '{print $1}')

    # Concatenate CAGE supported 5p transcript ids
    cat ./Data/Processed/Cage/$cat_name.tx.cage >> ./Data/.temp/ALL_catalogues.tx.cage

    # Report progress
    echo "Added $cat_name to ./Data/.temp/ALL_catalogues.tx.cage file"

    # Concatenate polyA supported 3p transcript ids
    cat ./Data/Processed/PolyA/$cat_name.tx.polyA >> ./Data/.temp/ALL_catalogues.tx.polyA

    # Report progress
    echo "Added $cat_name to ./Data/.temp/ALL_catalogues.tx.polyA file"

done < ./Data/Catalogues/catalogues_list
echo '~~NOTE: skipped gen.hg38 catalogue.~~'
echo '++ Done'

# Concatenate further with transcript ids originating from lncRNA spliced_mT slice.
## CAGE supported
echo '@ Concatenating with CAGE supported transcript ids from lncRNA spliced_mT slice...'
cat ./Data/.temp/ALL_catalogues.tx.cage ./Data/Processed/Cage/mT.tx.cage | sort | uniq > ./Data/Processed/Cage/ALL_catalogues+mT.tx.cage
echo '++ Done'
## polyA supported
echo '@ Concatenating with polyA supported transcript ids from lncRNA spliced_mT slice...'
cat ./Data/.temp/ALL_catalogues.tx.polyA ./Data/Processed/PolyA/mT.tx.polyA | sort | uniq > ./Data/Processed/PolyA/ALL_catalogues+mT.tx.polyA
echo '++ Done'
echo '################## DONE ##################'

##########################################################################################
# Extract transcripts ends from concatenated transcripts
##########################################################################################
# Extract 5p and 3p tx ends
echo '### Extracting transcripts ends...'
echo '@ Extracting 5p ends...'
cat ./Data/Processed/mT+Catalogues/ALL_catalogues+mT.bed12 | ./Utils/extractTranscriptEndsFromBed12.pl 5 > ./Data/.temp/ALL_catalogues+mT.5.bed
echo '++ Done'
echo '@ Extracting 3p ends...'
cat ./Data/Processed/mT+Catalogues/ALL_catalogues+mT.bed12 | ./Utils/extractTranscriptEndsFromBed12.pl 3 > ./Data/.temp/ALL_catalogues+mT.3.bed
echo '++ Done'
# Filter out ends that won't be used for anchoring and remove duplicated entries
## 5p ends
cat ./Data/.temp/ALL_catalogues+mT.5.bed | grep -Ff ./Data/Processed/Cage/ALL_catalogues+mT.tx.cage > ./Data/Processed/Extracted_ends/ALL_catalogues+mT.5.bed
## 3p ends
cat ./Data/.temp/ALL_catalogues+mT.3.bed | grep -Ff ./Data/Processed/PolyA/ALL_catalogues+mT.tx.polyA > ./Data/Processed/Extracted_ends/ALL_catalogues+mT.3.bed

###### WYWALIC
echo 'DIAGNOSTIC PART'
echo 'Raw 5p:'
cat ./Data/Processed/Extracted_ends/ALL_catalogues+mT.5.bed | wc -l
echo 'After deduplication 5p:'
cat ./Data/Processed/Extracted_ends/ALL_catalogues+mT.5.bed | awk '{print $4}' | sort | uniq | wc -l
echo '#########################'
echo 'Raw 3p:'
cat ./Data/Processed/Extracted_ends/ALL_catalogues+mT.3.bed | wc -l
echo 'After deduplication 3p:'
cat ./Data/Processed/Extracted_ends/ALL_catalogues+mT.3.bed | awk '{print $4}' | sort | uniq | wc -l

#cat ./Data/Processed/Extracted_ends/ALL_catalogues+mT.5.bed | awk '{if (count[$4] > 1) print $0; else if (count[$4] == 1) { print save[$4]; print $0; } else save[$4] = $0; count[$4]++; }'

# # Problem:
# 7786 transkryptów kompletnych na koncu 5p oraz
# 1954 transkryptów kompletnych na koncu 3p posiada duplikaty
# tzn. wystepuje inny transkrypt o tym samym tx.id, który ma inne współrzędne odpowiednio 5p i 3p.
# W części przypadków różnice te są bardzo małe ~ kilka zasad.
# Może rozwiązaniem byłby merge za z wykorzystaniem bedtoold merge.
# Duplicated and 5p supported only: 6718
# Duplicated and 3p supported only: 886
# Duplicated and supported on both ends: 1068


cat ./Data/Processed/Extracted_ends/ALL_catalogues+mT.5.bed | awk '{if (count[$4] > 1) print $0; else if (count[$4] == 1) { print save[$4]; print $0; } else save[$4] = $0; count[$4]++; }' | awk '{print $4}' | sort | uniq > ./Data/.temp/duplicated.5p.tx
cat ./Data/Processed/Extracted_ends/ALL_catalogues+mT.3.bed | awk '{if (count[$4] > 1) print $0; else if (count[$4] == 1) { print save[$4]; print $0; } else save[$4] = $0; count[$4]++; }' | awk '{print $4}' | sort | uniq > ./Data/.temp/duplicated.3p.tx

echo "Duplicated and 5p supported only: $(cat ./Data/.temp/duplicated.5p.tx | grep -vFf ./Data/.temp/duplicated.3p.tx | wc -l)"
echo "Duplicated and 3p supported only: $(cat ./Data/.temp/duplicated.3p.tx | grep -vFf ./Data/.temp/duplicated.5p.tx | wc -l)"
echo "Duplicated and supported on both ends: $(cat ./Data/.temp/duplicated.5p.tx | grep -Ff ./Data/.temp/duplicated.3p.tx | wc -l)"
###### WYWALIC

#echo '################## DONE ##################'





# ##########################################################################################
# # Anchor full-length transcript models
# ##########################################################################################

# echo '#####################'
# echo '### Anchoring concatenated catalogues...'
# echo '#####################'

# ## Concatenate tx.cage for all catalogues
# rm ./Data/Processed/Cage/ALL.tx.cage
# cat ./Data/Processed/Cage/*.tx.cage | sort | uniq > ./Data/Processed/Cage/ALL.tx.cage

# ## Concatenate tx.polyA for all catalogues
# rm ./Data/Processed/PolyA/ALL.tx.polyA
# cat ./Data/Processed/PolyA/*.tx.polyA | sort | uniq > ./Data/Processed/PolyA/ALL.tx.polyA

# ## Concatenate 5p ends
# rm ./Data/Processed/Extracted_ends/ALL.5.bed
# cat ./Data/Processed/Extracted_ends/*.5.bed | ./Utils/sortbed | uniq > ./Data/Processed/Extracted_ends/ALL.5.bed

# ## Concatenate 3p ends
# rm ./Data/Processed/Extracted_ends/ALL.3.bed
# cat ./Data/Processed/Extracted_ends/*.3.bed | ./Utils/sortbed | uniq > ./Data/Processed/Extracted_ends/ALL.3.bed

# ## Sort input gff
# cat ./Data/Processed/Catalogues_gff/ALL_catalogues.gff | sort -k1,1 -k12,12 -k4,4n -k5,5n "$@" | uniq > ./Data/Processed/Catalogues_gff/sorted_ALL_catalogues.gff

# ## Anchor TM
# ./Utils/anchorTranscriptsEnds.pl ./Data/Processed/Catalogues_gff/sorted_ALL_catalogues.gff ./Data/Processed/Cage/ALL.tx.cage ./Data/Processed/PolyA/ALL.tx.polyA ./Data/Processed/Extracted_ends/ALL.5.bed ./Data/Processed/Extracted_ends/ALL.3.bed > ./Data/Processed/Anchored_catalogues/anchored_ALL.gff



# #############################################
# # Run tmerge
# #############################################

# # Change dir for cleaner script
# cd ./Data/Processed/Anchored_catalogues/

# # Run Tmerge for each catalogue
# echo '#####################'
# echo '### Merging catalogues...'
# echo '#####################'

# cat anchored_fantomCat.hg38_polyexonic.gff anchored_cls.hg38.gff anchored_bigtrans.hg38.gff anchored_clsFL.hg38.gff anchored_gen.hg38.gff anchored_mitrans.hg38.gff anchored_noncode.hg38.gff anchored_pcConf.hg38.gff anchored_refseq.hg38.gff | ../../../Utils/sortgff | ../../../Utils/tmerge - --tmPrefix anch > ../../../Data/Processed/Merged_anchored_catalogues/anchored_tmerged_catalogues.gtf

# # Run Tmerge for concatenated catalogues
# echo '#####################'
# echo '### Merging concatenatd catalogues...'
# echo '#####################'

# cat anchored_ALL.gff | ../../../Utils/sortgff | ../../../Utils/tmerge - --tmPrefix anch > ../../../Data/Processed/Merged_anchored_catalogues/anchored_tmerged_ALL_catalogues.gtf

# # Change dir to previous state
# cd ../../..




# # Remove fake-exons
# #cat ./Data/Processed/Merged_anchored_catalogues/anchored_tmerged_catalogues.gtf | grep -Fv 'fa'
