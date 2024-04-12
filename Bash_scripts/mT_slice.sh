#!/usr/bin/env bash

#############################################
# Clean-up generated data
#############################################

echo 'Cleaning-up data from previous run..'
rm ./Data/Processed/mT/*
echo 'Done.'

#############################################
# Prepare lncRNA masterTable slice
#############################################

# Extract target_ids for "longNonCoding" target_category
echo 'Extracting target_ids...'
zcat ./Data/masterTable/Hv3_CLS3_targetDesign.gtf.gz | grep 'longNonCoding' | ./Utils/extract.gtf.tags.sh - target_id > ./Data/Processed/mT/longNonCoding.target.ids
echo 'Done.'

# Aquire transcript_ids for longNonCoding target_ids
## This step is needed, cuz exons do not contain info about targets
echo 'Extracting lncRNA related transcript_ids from mT...'
zcat ./Data/masterTable/Hv3_splicedmasterTable_refined.gtf.gz | grep -Ff ./Data/Processed/mT/longNonCoding.target.ids | ./Utils/extract.gtf.tags.sh - transcript_id > ./Data/Processed/mT/longNonCoding.transcript.ids
echo 'Done.'

# Use aquired longNonCoding transcript_ids to prepare mT slice (with exons)
echo 'Slicing mT...'
zcat ./Data/masterTable/Hv3_splicedmasterTable_refined.gtf.gz | grep -Ff ./Data/Processed/mT/longNonCoding.transcript.ids > ./Data/Processed/mT/Hv3_splicedmasterTable_refined_lncRNA.gtf
echo 'Done.'



