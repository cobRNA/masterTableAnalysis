cat ./Data/Processed/Anchored_catalogues/*.gff | sort -k1,1 -k4,4n -k5,5n "$@" | ./Utils/tmerge - > ./Data/Processed/Merged_anchored_catalogues/lncrna.anchor.merged.gtf

# | grep -v 'fakeExon "yes";'

# Stary fantom
cat ./Data/Catalogues/cls.hg38.bed12 ./Data/Catalogues/fantomCat.hg38.bed12 | ./Utils/bed12togff | sort -k1,1 -k4,4n -k5,5n "$@" | ./Utils/tmerge - > ./Random/Remove/tmerge_old_fantom.gtf
## Only polyexonic fantom
cat ./Data/Processed/Catalogues_gff/fantomCat.hg38.gff | awk '{if (count[$12] > 1) print $0; else if (count[$12] == 1) { print save[$12]; print $0; } else save[$12] = $0; count[$12]++; }' | cat - ./Data/Processed/Catalogues_gff/cls.hg38.gff | sort -k1,1 -k4,4n -k5,5n "$@" | ./Utils/tmerge - > ./Random/Remove/tmerge_old_fantom.gtf
## Only polyexonic fantom and excluding chr7_KI270803v1_alt chr14_GL000009v2_random chr4_GL000008v2_random
cat ./Data/Processed/Catalogues_gff/fantomCat.hg38.gff | grep -v -E 'chr7_KI270803v1_alt|chr14_GL000009v2_random|chr4_GL000008v2_random' | awk '{if (count[$12] > 1) print $0; else if (count[$12] == 1) { print save[$12]; print $0; } else save[$12] = $0; count[$12]++; }' | cat - ./Data/Processed/Catalogues_gff/cls.hg38.gff | sort -k1,1 -k4,4n -k5,5n "$@" | ./Utils/tmerge - > ./Random/Remove/tmerge_old_fantom.gtf



## Only polyexonic fantom and cls
cat ./Data/Processed/Catalogues_gff/fantomCat.hg38.gff ./Data/Processed/Catalogues_gff/cls.hg38.gff | awk '{if (count[$12] > 1) print $0; else if (count[$12] == 1) { print save[$12]; print $0; } else save[$12] = $0; count[$12]++; }' | sort -k1,1 -k4,4n -k5,5n "$@" | ./Utils/tmerge - > ./Random/Remove/tmerge_old_fantom.gtf

## With gencode
cat ./Data/Catalogues/gen.hg38.bed12 ./Data/Catalogues/fantomCat.hg38.bed12 | ./Utils/bed12togff | sort -k1,1 -k4,4n -k5,5n "$@" | ./Utils/tmerge - > ./Random/Remove/tmerge_old_fantom_gencode.gtf
## Only polyexonic fantom
cat ./Data/Processed/Catalogues_gff/fantomCat.hg38.gff | awk '{if (count[$12] > 1) print $0; else if (count[$12] == 1) { print save[$12]; print $0; } else save[$12] = $0; count[$12]++; }' | cat - ./Data/Catalogues/gen.hg38.bed12 | sort -k1,1 -k4,4n -k5,5n "$@" | ./Utils/tmerge - > ./Random/Remove/tmerge_old_fantom.gtf

## Gencode with cls
cat ./Data/Catalogues/gen.hg38.bed12 ./Data/Catalogues/cls.hg38.bed12 | ./Utils/bed12togff | sort -k1,1 -k4,4n -k5,5n "$@" | ./Utils/tmerge - > ./Random/Remove/tmerge_old_fantom_gencode.gtf



# Z nowym Fantomem
cat ./Data/Catalogues/cls.hg38.bed12 | ./Utils/bed12togff | cat - ./Data/Catalogues/FANTOM_CAT.lv4_stringent.only_lncRNA.gtf | sort -k1,1 -k4,4n -k5,5n "$@" | ./Utils/tmerge - > ./Random/Remove/tmerge_new_fantom.gtf



############
# Finalny Tmerge
############
cat ./Data/Processed/Catalogues_gff/fantomCat.hg38.gff | grep -v -E 'chr7_KI270803v1_alt|chr14_GL000009v2_random|chr4_GL000008v2_random' | awk '{if (count[$12] > 1) print $0; else if (count[$12] == 1) { print save[$12]; print $0; } else save[$12] = $0; count[$12]++; }' | cat - ./Data/Processed/Catalogues_gff/cls.hg38.gff ./Data/Processed/Catalogues_gff/bigtrans.hg38.gff ./Data/Processed/Catalogues_gff/clsFL.hg38.gff ./Data/Processed/Catalogues_gff/gen.hg38.gff ./Data/Processed/Catalogues_gff/mitrans.hg38.gff ./Data/Processed/Catalogues_gff/noncode.hg38.gff ./Data/Processed/Catalogues_gff/pcConf.hg38.gff ./Data/Processed/Catalogues_gff/refseq.hg38.gff | sort -k1,1 -k4,4n -k5,5n "$@" | ./Utils/tmerge - > ./Random/Remove/tmerge_old_fantom.gtf






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