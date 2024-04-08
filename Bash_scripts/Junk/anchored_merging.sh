#!/usr/bin/env bash


#############################################
# Clean-up generated data
#############################################

rm ./Data/Processed/Catalogues_gff/*
rm ./Data/Processed/Anchored_catalogues/*
# Clean-up report file
> ./Data/Reports/anchoring

#############################################
# Anchor full-length transcript models
#############################################

while read catalogue || [ -n "$catalogue" ]  # solution to no terminal newline problem
do  
    # Create base catalogue name
    cat_name=$(echo $catalogue | awk -F'.' '{print $1}')

    # Report progress
    echo $cat_name

    # Convert catalogue to required gff format
    cat ./Data/Catalogues/$catalogue.bed12 | ./Utils/bed12togff > ./Data/Processed/Catalogues_gff/$catalogue.gff
       
    # Sort catalogue before anchoring
    cat ./Data/Processed/Catalogues_gff/$catalogue.gff | sort -k1,1 -k12,12 -k4,4n -k5,5n "$@" > ./Data/Processed/Catalogues_gff/${catalogue}_sorted.gff

    # Run anchoring on sorted content
    ## Use 3p and 5p ends as clusters - not merging anything at all
    ./Utils/anchorTranscriptsEnds.pl ./Data/Processed/Catalogues_gff/${catalogue}_sorted.gff ./Data/Processed/Cage/$cat_name.tx.cage ./Data/Processed/PolyA/$cat_name.tx.polyA ./Data/Processed/Extracted_ends/${cat_name}.5.bed ./Data/Processed/Extracted_ends/${cat_name}.3.bed > ./Data/Processed/Anchored_catalogues/anchored_${catalogue}.gff                     

    # Look for corrupted entries and create report file
    echo '##############################' >> ./Data/Reports/anchoring
    echo 'Looking for corrupted entries...' >> ./Data/Reports/anchoring
    cat ./Data/Processed/Anchored_catalogues/anchored_${catalogue}.gff  | awk '$4>$5' >> ./Data/Reports/anchoring
    echo '##############################' >> ./Data/Reports/anchoring

done < ./Data/Catalogues/catalogues_list

# Preview report file
cat ./Data/Reports/anchoring



#############################################
# Run tmerge
#############################################


