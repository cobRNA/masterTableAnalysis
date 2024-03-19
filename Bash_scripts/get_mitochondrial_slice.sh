#!/usr/bin/env bash


#---------------------
# Aquire IDs
#---------------------

#----- Gencode.v45
# Removed version no. from gene_id and transcript_id ::: ENSG00000290825.1 -> ENSG00000290825
# Removed _PAR_Y from the end of gene_id and transcript_id
zcat ./Data/Gencode_annotations/gencode.v45.annotation.gtf.gz | awk 'NR>5{print$0}' | awk '$3=="transcript"' | ./Utils/extract.gtf.tags.sh - gene_id,transcript_id | awk -F'\t' 'gsub(/\.[0-9]+/,""){print $1"\t"$2}' | sed 's/_PAR_Y//g' > ./Data/Gencode_annotations/gencode.v45.annotation.ids


#----- Gencode.v44
# Removed version no. from gene_id and transcript_id ::: ENSG00000290825.1 -> ENSG00000290825
# Removed _PAR_Y from the end of gene_id and transcript_id
zcat ./Data/Gencode_annotations/gencode.v44.annotation.gtf.gz | awk 'NR>5{print$0}' | awk '$3=="transcript"' | ./Utils/extract.gtf.tags.sh - gene_id,transcript_id | awk -F'\t' 'gsub(/\.[0-9]+/,""){print $1"\t"$2}' | sed 's/_PAR_Y//g' > ./Data/Gencode_annotations/gencode.v44.annotation.ids



#----- Gencode.v27
# Removed version no. from gene_id and transcript_id ::: ENSG00000290825.1 -> ENSG00000290825
# Removed _PAR_Y from the end of gene_id and transcript_id
zcat ./Data/Gencode_annotations/gencode.v27.annotation.gtf.gz | awk 'NR>5{print$0}' | awk '$3=="transcript"' | ./Utils/extract.gtf.tags.sh - gene_id,transcript_id | awk -F'\t' 'gsub(/\.[0-9]+/,""){print $1"\t"$2}' | sed 's/_PAR_Y//g' > ./Data/Gencode_annotations/gencode.v27.annotation.ids



#----- Hv3_enhancedCLS3_refined_+gencodev43.loci.refmerged
# Removed version no. from gene_id and transcript_id ::: ENSG00000290825.1 -> ENSG00000290825
# Removed _PAR_Y from the end of gene_id and transcript_id
zcat ./Data/masterTable/Hv3_enhancedCLS3_refined_+gencodev43.loci.refmerged.gff.gz | awk '$3=="transcript"' | ./Utils/extract.gtf.tags.sh - gene_id,transcript_id | awk -F'\t' 'gsub(/\.[0-9]+/,""){print $1"\t"$2}' | sed 's/_PAR_Y//g' > ./Data/masterTable/Hv3_enhancedCLS3_refined_+gencodev43.loci.refmerged.ids


#----- IMPI
# Downloaded latest version as IMPI 2021 version (IMPI-2021-Q4pre) xlsx file from:
# https://www.mrc-mbu.cam.ac.uk/research-resources-and-facilities/impi
# Saved only ENSEMBL column with Verifid Mitochondrial IMPI Class as impi-2021-q4pre-20211001-dist_0.ids


#----- MitoCarta3.0
# Downloaded latest version as Human.MitoCarta3.0.xls file from:
# https://www.broadinstitute.org/mitocarta/mitocarta30-inventory-mammalian-mitochondrial-proteins-and-pathways
# Saved only EnsemblGeneID_mapping_version_20200130 column from A Human MitoCarta3.0 sheet as csv file
# Some lines contained multiple IDs separated by |
# I splitted them using:
cat ./Data/Mitochondrial_catalogues/Human.MitoCarta3.0.csv | sed -z 's/|/\n/g' > ./Data/Mitochondrial_catalogues/Human.MitoCarta3.0.ids


#----- IMPI and MitoCarta3.0
# ID files from both catalogues where concatenated, duplicates were removed using:
cat ./Data/Mitochondrial_catalogues/impi-2021-q4pre-20211001-dist_0.ids ./Data/Mitochondrial_catalogues/Human.MitoCarta3.0.ids | sort | uniq > ./Data/Mitochondrial_catalogues/impi_and_MitoCarta.ids


#---------------------
# Prepare gencode.v45/v44/v27.annotations slices
#---------------------

#----- Gencode.v45/v44/v27.annotations
# Comprehensive gene annotation gtf file was downloaded from https://www.gencodegenes.org/human/

for id_file in ./Data/Mitochondrial_catalogues/*ids
    do
        filename=$(basename -- "$id_file")
        extension="${filename##*.}"
        filename="${filename%.*}"
        echo $filename
        # Prepare slice
        for gen_annot in ./Data/Gencode_annotations/*gtf.gz
        do
            annot_filename=$(basename -- "$gen_annot")
            annot_extension="${annot_filename##*.}"
            annot_filename="${annot_filename%.*}"
            echo $annot_filename
            zcat ./Data/Gencode_annotations/${annot_filename} | grep -f $id_file > ./Data/Processed/Mitochondrial/${filename}.${annot_filename}_slice.gtf
            # Extract slice ids
            cat ./Data/Processed/Mitochondrial/${filename}.${annot_filename}_slice.gtf | awk 'NR>5{print$0}' | awk '$3=="transcript"' | ./Utils/extract.gtf.tags.sh - gene_id,transcript_id | awk -F'\t' 'gsub(/\.[0-9]+/,""){print $1"\t"$2}' | sed 's/_PAR_Y//g' > ./Data/Processed/Mitochondrial/${filename}.${annot_filename}_slice.ids
            # Extract missing ids
            cat ./Data/Processed/Mitochondrial/${filename}.${annot_filename}_slice.gtf | awk '{print $10}' | sed 's/"//g ; s/;//g'| awk -F'.' '{print $1}' | sed 's/_PAR_Y//g' | sort | uniq | grep -vf - $id_file > ./Data/Processed/Mitochondrial/${filename}.${annot_filename}_slice.missing.ids      
        done;
    done;


#---------------------
# Prepare masterTable slices
#---------------------

#----- masterTable annotation
# Hv3_enhancedCLS3_refined_+gencodev43.loci.refmerged gtf file was downloaded from https://github.com/guigolab/gencode-cls-master-table/releases


for id_file in ./Data/Mitochondrial_catalogues/*ids
    do
        filename=$(basename -- "$id_file")
        extension="${filename##*.}"
        filename="${filename%.*}"
        # Prepare slice
        zcat ./Data/masterTable/Hv3_enhancedCLS3_refined_+gencodev43.loci.refmerged.gff.gz | grep -f $id_file > ./Data/Processed/Mitochondrial/${filename}_mT_Hv3_enhancedCLS3_slice.gtf
        # Extract slice ids
        cat ./Data/Processed/Mitochondrial/${filename}_mT_Hv3_enhancedCLS3_slice.gtf | awk '$3=="transcript"' | ./Utils/extract.gtf.tags.sh - gene_id,transcript_id | awk -F'\t' 'gsub(/\.[0-9]+/,""){print $1"\t"$2}' | sed 's/_PAR_Y//g' > ./Data/Processed/Mitochondrial/${filename}_mT_Hv3_enhancedCLS3_slice.ids
        # Extract missing ids
        cat ./Data/Processed/Mitochondrial/${filename}_mT_Hv3_enhancedCLS3_slice.gtf | awk '{print $10}' | sed 's/"//g ; s/;//g'| awk -F'.' '{print $1}' | sed 's/_PAR_Y//g' | sort | uniq | grep -vf - $id_file > ./Data/Processed/Mitochondrial/${filename}_mT_Hv3_enhancedCLS3_slice.missing.ids
    done;


#---------------------
# Calculate stats for Gencode.v45.annotation and masterTable annotation slices
#---------------------

# printf 'Stats calculated for Gencode.v45 and mT_Hv3 slices\n' > ./Data/Processed/Mitochondrial/mito_stats # Remove file content
# for annotation_slice in ./Data/Processed/Mitochondrial/*.gtf
#     do
#         filename=$(basename -- "$annotation_slice")
#         extension="${filename##*.}"
#         filename="${filename%.*}"
#         echo $filename
#         missing_gene_ids=$(cat ./Data/Processed/Mitochondrial/${filename}.missing.ids | sort | uniq | wc -l)
#         slice_gene_ids=$(cat $annotation_slice | awk '{print $10}' | sort | uniq | wc -l)
#         transcript_ids=$(cat $annotation_slice | awk '$3=="transcript"' | awk '{print $12}' | sort | uniq | wc -l)

#         printf '%s:\n' "$filename" >> ./Data/Processed/Mitochondrial/mito_stats
#         printf 'Missing gene_ids: %s\n' "$missing_gene_ids" >> ./Data/Processed/Mitochondrial/mito_stats
#         printf 'Common gene_ids: %s\n' "$slice_gene_ids" >> ./Data/Processed/Mitochondrial/mito_stats
#         printf 'Common transcript_ids: %s\n' "$transcript_ids" >> ./Data/Processed/Mitochondrial/mito_stats
#         printf '%s\n' "-----------------------" >> ./Data/Processed/Mitochondrial/mito_stats
#     done

# "ENST00000482405.7_PAR_Y"; <- id sztucznie zwiększające liczbę genów

# sed "s/_PAR_Y//g"   
#  awk '$1 !~ /_PAR_Y/' | awk -F'.' '{print$1}' |

# echo $filename
#         gene_ids=$(cat $annotation_slice | awk '$3=="transcript"' | awk '{print $10}' | sort | uniq | wc -l)
#         missing_gene_ids=$(cat ${filename}.missing.ids | awk '$3=="transcript"' | awk '{print $10}' | sort | uniq | wc -l)
#         transcript_ids=$(cat $annotation_slice | awk '$3=="transcript"' | awk '{print $12}' | sort | uniq | wc -l)



#---------------------
# Prepare data for Venn diagram in R
#---------------------








