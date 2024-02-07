###############
# Merging TMs from RNA catalogues using anchored merging.
###############

# -----------------------------
# raw reads
# -----------------------------

# Moimi raw readsami w bed, są chyba katalogi

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


# --------------------
# stats and plot <- ze starego skryptu
# --------------------

while read lab
do
    for end in 3
    do 
        close=`cat ./Data/Processed/PolyA/$lab.$end.bed.vspolyAsignals.bedtsv | cut -f4 | sort |uniq | wc -l`
        total=`cat ./Data/Processed/Extracted_ends/$lab.$end.bed | cut -f4 | sort | uniq | wc -l`
        echo -e "$lab\t$total\t$close\t$end" | awk '{print $0"\t"$3/$2}'    
    done 
done < ./Data/Others/samples.tsv | sed 's/gen/GENCODE/g'| sed 's/refseq/RefSeq/g'|sed 's/fantomCat/FANTOM CAT/g'| sed 's/mitrans/MiTranscriptome/g'| sed 's/bigtrans/BIGTranscriptome/g'| sed 's/noncode/NONCODE/g'| sed 's/cls/CLS/g' | sed 's/CLS+FL/CLS FL/g' > ./Data/Processed/PolyA/annots.polyAsignals.stats.tsv 

#| sed 's/pcConf/Protein coding/g'| sed 's/GENCODE+FL/CLS FL/g'

echo "
library(ggplot2)
library(scales)
cbPalette <- c(\"#a6761d\", \"#e6ab02\", \"#FF7F00\", \"#984EA3\", \"#4DAF4A\", \"#377EB8\", \"#E41A1C\") #, \"#e7298a\",\"#999999\"
plot <- read.table(\"./Data/Processed/PolyA/annots.polyAsignals.stats.tsv\", header=F, as.is=T, sep=\"\t\")
colnames(plot)<-c(\"gene\", \"total\", \"count\", \"end\",\"prop\")
plot\$gene=factor(plot\$gene, levels=c(\"CLS\",\"GENCODE\", \"BIGTranscriptome\", \"RefSeq\", \"FANTOM CAT\", \"MiTranscriptome\", \"NONCODE\")) #, \"CLS FL\", \"Protein coding\"
pdf(\"./Plots/human.annots.polyAsignals.stats.pdf\", width=8, height=6)
#setEPS()
#postscript(\"human.annots.polyAsignals.stats.eps\", family=\"serif\", width=8, height=10)
ggplot(data=plot, aes(x=gene, y=prop)) +
geom_bar(stat=\"identity\", fill=cbPalette) +
ylab(\"% PAS(+) 3\' ends\") +
theme_bw(base_size = 28) +
xlab(\"\") +
scale_y_continuous(labels=percent, limits=c(0, 1)) +
coord_flip() +
geom_text(position = \"stack\", aes(x = gene, y = prop, label = comma(count), hjust = -0.1, vjust = 0.5), size=7) +
theme(axis.text.x  = element_text(angle=45, vjust=0.5),
axis.ticks.y = element_blank(),
axis.text.y  = element_text(angle=0, hjust=0.5)) +
theme(axis.line.x = element_line(colour = \"black\"),
axis.line.y = element_line(colour = \"black\"),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
panel.border = element_blank(),
panel.background = element_blank(),
strip.background = element_rect(colour=\"black\",fill=\"white\")) 
#+facet_wrap(~ end) 
dev.off()
" | R --slave


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


# --------------------
# stats and plot <- ze starego skryptu
# --------------------

while read lab
do
    for end in 5
    do 
        close=`cat ./Data/Processed/Cage/$lab.$end.bed.vsCage.fantom.bedtsv | cut -f4 | sort | uniq | wc -l`
        total=`cat ./Data/Processed/Extracted_ends/$lab.$end.bed | cut -f4 | sort | uniq | wc -l`
        echo -e "$lab\t$total\t$close\t$end" | awk '{print $0"\t"$3/$2}'    
    done 
done < ./Data/Others/samples.tsv | sed 's/gen/GENCODE/g'| sed 's/refseq/RefSeq/g'|sed 's/fantomCat/FANTOM CAT/g'| sed 's/mitrans/MiTranscriptome/g'| sed 's/bigtrans/BIGTranscriptome/g'| sed 's/noncode/NONCODE/g'| sed 's/cls/CLS/g' | sed 's/CLS+FL/CLS FL/g' > ./Data/Processed/Cage/annots.vsCage.fantom.stats.tsv 

echo "
library(ggplot2)
library(scales)
cbPalette <- c(\"#a6761d\", \"#e6ab02\", \"#FF7F00\", \"#984EA3\", \"#4DAF4A\", \"#377EB8\", \"#E41A1C\") #,\"#999999\", \"#e7298a\"
plot <- read.table(\"./Data/Processed/Cage/annots.vsCage.fantom.stats.tsv\", header=F, as.is=T, sep=\"\t\")
colnames(plot)<-c(\"gene\", \"total\", \"count\", \"end\",\"prop\")
plot\$gene=factor(plot\$gene, levels=c(\"CLS\",\"GENCODE\",\"BIGTranscriptome\",\"RefSeq\",\"FANTOM CAT\",\"MiTranscriptome\",\"NONCODE\")) #, \"Protein coding\", \"CLS FL\"
pdf(\"./Plots/human.annots.vsCage.fantom.stats.pdf\", width=7, height=6)
#setEPS()
#postscript(\"human.annots.vsCage.fantom.stats.eps\", family=\"serif\", width=6, height=10)
ggplot(data=plot, aes(x=gene, y=prop)) +
geom_bar(stat=\"identity\", fill=cbPalette) +
ylab(\"% CAGE(+) 5\' ends\") +
theme_bw(base_size = 28) +
xlab(\"\") +
scale_y_continuous(labels=percent, lim=c(0,1)) +
coord_flip() +
geom_text(position = \"stack\", aes(x = gene, y = prop, label = comma(count), hjust = -0.1, vjust = 0.5), size=7) +
theme(axis.text.x  = element_text(angle=45, vjust=0.5)) +
theme(axis.line.x = element_line(colour = \"black\"),
axis.line.y = element_line(colour = \"black\"),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
panel.border = element_blank(),
panel.background = element_blank(),
strip.background = element_rect(colour=\"black\",fill=\"white\"))
#+facet_wrap(~ end) 
dev.off()
" | R --slave

#, axis.ticks.y = element_blank(), axis.text.y  = element_blank()



# ------------------------------
# get the read ids for each category
# ------------------------------


# Convert catalogues in bed12 to gff
for file in `ls ./Data/Catalogues/*.hg38.bed12`
    do
    name=`basename $file | awk -F "." '{print $1"."$2}'`
    ./Utils/bed12togff ./Data/Catalogues/${name}.bed12 > ./Data/Catalogues/Gff/${name}.gff
    echo "$name converted."
    done

# Get ids
# UWAGA!! ZMIENIŁEM
# cat \$name.bed| cut -f4| sort| uniq ..........
# NA
# cat ./Data/Catalogues/$name.hg38.bed12 | cut -f4 | sort | uniq ..........

while read sampl
do
    name=$sampl
    echo $name
    cat ./Data/Catalogues/$name.hg38.bed12 | cut -f4 | sort | uniq > ./Data/Processed/Read_ids/${name}.readID
    cat ./Data/Processed/Cage/$name.5.bed.vsCage.fantom.bedtsv | cut -f4 | sed 's/,/\n/g' | sort | uniq > ./Data/Processed/Read_ids/${name}.cage.readID
    cat ./Data/Processed/PolyA/$name.3.bed.vspolyAsignals.bedtsv | cut -f4 | sed 's/,/\n/g' | sort | uniq > ./Data/Processed/Read_ids/${name}.polyA.readID
    cat ./Data/Catalogues/Gff/$name.hg38.gff | awk '{print $10}'| sed 's/\"//g'| sed 's/;//g'| sort | uniq -c | awk '$1==1{print $2}' > ./Data/Processed/Read_ids/${name}.mono.readID
    cat ./Data/Catalogues/Gff/$name.hg38.gff | awk '{print $10}'| sed 's/\"//g'| sed 's/;//g'| sort | uniq -c | awk '$1>1{print $2}' > ./Data/Processed/Read_ids/${name}.spliced.readID 
    cat ./Data/Catalogues/Gff/$name.hg38.gff | awk '{print $10}'| sed 's/\"//g'| sed 's/;//g'| sort | uniq > ./Data/Processed/Read_ids/${name}.all.readID 
done < ./Data/Others/samples.tsv

# --------------------------------
# stats polyA reads
# --------------------------------

# UWAGA
# ZAMIENIŁEM
# fgrep -f $name.polyAreads.list $name.$typ.readID > ...............
# NA
# fgrep -f ./Data/Others/hg38.polyAsignals.bed $name.$typ.readID > ...............

for typ in all spliced mono
do 
    echo $typ
    while read sampl
    do 
        name=`echo $sampl`
        total=`cat ./Data/Processed/Read_ids/$name.readID | wc -l`
        All=`cat ./Data/Processed/Read_ids/$name.$typ.readID | wc -l`
        fgrep -f ./Data/Others/hg38.polyAsignals.bed ./Data/Processed/Read_ids/$name.$typ.readID > ./Data/Processed/Read_ids/$name.$typ.polyAreads.readID
        fgrep -f ./Data/Processed/Read_ids/$name.cage.readID ./Data/Processed/Read_ids/$name.$typ.readID > ./Data/Processed/Read_ids/$name.$typ.cage.polyAreads.readID
        fgrep -f ./Data/Processed/Read_ids/$name.$typ.polyAreads.readID ./Data/Processed/Read_ids/$name.$typ.cage.polyAreads.readID > ./Data/Processed/Read_ids/$name.$typ.FL.polyAreads.readID

        PolyA=`fgrep -vf ./Data/Processed/Read_ids/$name.$typ.FL.polyAreads.readID ./Data/Processed/Read_ids/$name.$typ.polyAreads.readID | wc -l`
        CAGE=`fgrep -vf ./Data/Processed/Read_ids/$name.$typ.FL.polyAreads.readID ./Data/Processed/Read_ids/$name.$typ.cage.polyAreads.readID | wc -l`
        FL=`cat ./Data/Processed/Read_ids/$name.$typ.FL.polyAreads.readID | wc -l`
        supported=$(($PolyA+$CAGE+$FL))
        NA=`expr $All - $supported` 
        
        #PolyAprop=\`echo \$PolyA / \$All| bc -l\`
        #CAGEprop=\`echo \$CAGE / \$All| bc -l\`
        #FLprop=\`echo \$FL / \$All| bc -l\`
        #NAprop=\`echo \$NA / \$All| bc -l\`
        
        PolyAprop=`echo $PolyA / $total| bc -l`
        CAGEprop=`echo $CAGE / $total| bc -l`
        FLprop=`echo $FL / $total| bc -l`
        NAprop=`echo $NA / $total| bc -l`
        
        echo -e "$lab\tHpreCap\t0+\t$sampl\tcageAndPolyA\t$FL\t$FLprop\t$typ\n$lab\tHpreCap\t0+\t$sampl\tcageOnly\t$CAGE\t$CAGEprop\t$typ\n$lab\tHpreCap\t0+\t$sampl\tpolyAOnly\t$PolyA\t$PolyAprop\t$typ\n$lab\tHpreCap\t0+\t$sampl\tnoCageNoPolyA\t$NA\t$NAprop\t$typ"
    done < ./Data/Others/samples.tsv  > ./Data/Processed/PolyA/samples.end.completenes.raw.reads.stats.polyAreadslist.$typ.tsv
done  


