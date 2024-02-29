#!/usr/bin/env bash

###############
# Anchor full-length transcript models.
###############

for catalogue in ./Data/Catalogues/*bed12
    do
        # Get filename
        filename=$(basename -- "$catalogue")
        extension="${filename##*.}"
        filename="${filename%.*}"
        echo $filename

        # Convert catalogue to required gff format
        cat $catalogue | ./Utils/bed12togff > ./Data/Processed/Catalogues_gff/$filename.gff

        # Create clusters for each end
        cat_name=$(echo $filename | awk -F'.' '{print $1}')
        ## 5' (cage)
        bedtools merge -c 4 -o collapse -s -d 5 -i ./Data/Processed/Cage/$cat_name.5.bed.vsCage.fantom.bedtsv > ./Data/Processed/Clusters/${cat_name}_5p_clusters.bed
        ## 3' (polyA)
        bedtools merge -c 4 -o collapse -s -d 5 -i ./Data/Processed/PolyA/$cat_name.3.bed.vspolyAsignals.bedtsv > ./Data/Processed/Clusters/${cat_name}_3p_clusters.bed

        # Run anchoring
        ./Utils/anchorTranscriptsEnds.pl ./Data/Processed/Catalogues_gff/$filename.gff ./Data/Processed/Cage/$cat_name.tx.cage ./Data/Processed/PolyA/$cat_name.tx.polyA ./Data/Processed/Clusters/${cat_name}_5p_clusters.bed ./Data/Processed/Clusters/${cat_name}_3p_clusters.bed > ./Data/Processed/Anchored_catalogues/anchored_${filename}.gff                     
        done;







###############
# Merge catalogues with anchored FL trancript models.
# Remove fake exons after merging.
###############

cat ./Data/Processed/Anchored_catalogues/*.gff | sort -k1,1 -k4,4n -k5,5n | ./Utils/tmerge - > ./Data/Processed/Merged_anchored_catalogues/lncrna.anchor.merged.gtf

# | grep -v 'fakeExon "yes";'





# for catalogue in ./Data/Catalogues/*bed12
#     do
#         # Get filename
#         filename=$(basename -- "$catalogue")
#         extension="${filename##*.}"
#         filename="${filename%.*}"
#         echo $filename

#         # Convert catalogue to required gff format
#         cat $catalogue | ./Utils/bed12togff > ./Data/Processed/Catalogues_gff/$filename.gff

#         # Create clusters for each end
#         cat_name=$(echo $filename | awk -F'.' '{print $1}')
#         ## 5' (cage)
#         bedtools merge -d 5 -s -c 4,6 -o distinct -i ./Data/Processed/Cage/$cat_name.5.bed.vsCage.fantom.bedtsv > ./Data/Processed/Clusters/${cat_name}_5p_clusters.bed
#         ## 3' (polyA)
#         bedtools merge -d 5 -s -c 4,6 -o distinct -i ./Data/Processed/PolyA/$cat_name.3.bed.vspolyAsignals.bedtsv > ./Data/Processed/Clusters/${cat_name}_3p_clusters.bed

#         # Run anchoring
#         ./Utils/anchorTranscriptsEnds.pl ./Data/Processed/Catalogues_gff/$filename.gff ./Data/Processed/Cage/$cat_name.tx.cage ./Data/Processed/PolyA/$cat_name.tx.polyA ./Data/Processed/Clusters/${cat_name}_5p_clusters.bed ./Data/Processed/Clusters/${cat_name}_3p_clusters.bed > ./Data/Processed/Anchored_catalogues/anchored_${filename}.gff                     
#         done;