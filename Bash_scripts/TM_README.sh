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
    done < ./Data/Source/ends.dist.tsv
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
       cat ./Data/Processed/Extracted_ends/$name.$end.bed | ./Utils/sortbed | ./Utils/jlagarde/bin.bkp/bedtools2/bin/bedtools slop -s -l 5 -r 5 -i stdin -g ./Data/Source/chromInfo_hg38.txt | ./Utils/jlagarde/bin.bkp/bedtools2/bin/bedtools intersect -u -s -a stdin -b ./Data/Source/hg38.polyAsignals.bed > ./Data/Processed/PolyA/$name.$end.bed.vspolyAsignals.bedtsv
    done < ./Data/Source/ends.dist.tsv
done < ./Data/Source/samples.tsv


# --------------------
# stats and plot
# --------------------

while read lab
do
    for end in 3
    do 
        close=`cat ./Data/Processed/PolyA/$lab.$end.bed.vspolyAsignals.bedtsv | cut -f4 | sort | uniq | wc -l`
        total=`cat ./Data/Processed/Extracted_ends/$lab.$end.bed | cut -f4 | sort | uniq | wc -l`
        echo -e "$lab\t$total\t$close\t$end" | awk '{print $0"\t"$3/$2}'    
    done 
done < ./Data/Source/samples.tsv | sed 's/gen/GENCODE/g' | sed 's/refseq/RefSeq/g' | sed 's/fantomCat/FANTOM CAT/g' | sed 's/mitrans/MiTranscriptome/g' | sed 's/bigtrans/BIGTranscriptome/g' | sed 's/noncode/NONCODE/g' | sed 's/cls/GENCODE+/g' | sed 's/pcConf/Protein coding/g' | sed 's/GENCODE+FL/CLS FL/g' > ./Data/Processed/Plots_input/annots.polyAsignals.stats.tsv 

#| sed 's/pcConf/Protein coding/g'| sed 's/GENCODE+FL/CLS FL/g'

echo "
library(ggplot2)
library(scales)
cbPalette <- c(\"#a6761d\", \"#e6ab02\", \"#FF7F00\", \"#984EA3\", \"#4DAF4A\", \"#377EB8\", \"#E41A1C\", \"#e7298a\",\"#999999\") #
plot <- read.table(\"./Data/Processed/Plots_input/annots.polyAsignals.stats.tsv\", header=F, as.is=T, sep=\"\t\")
colnames(plot)<-c(\"gene\", \"total\", \"count\", \"end\",\"prop\")
plot\$gene=factor(plot\$gene, levels=c(\"GENCODE+\",\"GENCODE\", \"BIGTranscriptome\", \"RefSeq\", \"FANTOM CAT\", \"MiTranscriptome\", \"NONCODE\", \"CLS FL\", \"Protein coding\")) #
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
#axis.ticks.y = element_blank(), #axis.text.y  = element_text(angle=0, hjust=0.5)
) +
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
        cat ./Data/Processed/Extracted_ends/$name.$end.bed | ./Utils/sortbed | ./Utils/jlagarde/bin.bkp/bedtools2/bin/bedtools slop -s -l 50 -r 50 -i stdin -g ./Data/Source/chromInfo_hg38.txt | ./Utils/jlagarde/bin.bkp/bedtools2/bin/bedtools intersect -u -s -a stdin -b ./Data/Source/hg38_fair+new_CAGE_peaks_phase1and2.bed > ./Data/Processed/Cage/$name.$end.bed.vsCage.fantom.bedtsv
    done < ./Data/Source/ends.dist.tsv
done < ./Data/Source/samples.tsv


# --------------------
# stats and plot
# --------------------

while read lab
do
    for end in 5
    do 
        close=`cat ./Data/Processed/Cage/$lab.$end.bed.vsCage.fantom.bedtsv | cut -f4 | sort | uniq | wc -l`
        total=`cat ./Data/Processed/Extracted_ends/$lab.$end.bed | cut -f4 | sort | uniq | wc -l`
        echo -e "$lab\t$total\t$close\t$end" | awk '{print $0"\t"$3/$2}'    
    done 
done < ./Data/Source/samples.tsv | sed 's/gen/GENCODE/g' | sed 's/refseq/RefSeq/g' | sed 's/fantomCat/FANTOM CAT/g' | sed 's/mitrans/MiTranscriptome/g' | sed 's/bigtrans/BIGTranscriptome/g' | sed 's/noncode/NONCODE/g' | sed 's/cls/GENCODE+/g' | sed 's/pcConf/Protein coding/g' | sed 's/GENCODE+FL/CLS FL/g' > ./Data/Processed/Plots_input/annots.vsCage.fantom.stats.tsv 

echo "
library(ggplot2)
library(scales)
cbPalette <- c(\"#a6761d\", \"#e6ab02\", \"#FF7F00\", \"#984EA3\", \"#4DAF4A\", \"#377EB8\", \"#E41A1C\",\"#999999\", \"#e7298a\") #
plot <- read.table(\"./Data/Processed/Plots_input/annots.vsCage.fantom.stats.tsv\", header=F, as.is=T, sep=\"\t\")
colnames(plot)<-c(\"gene\", \"total\", \"count\", \"end\",\"prop\")
plot\$gene=factor(plot\$gene, levels=c(\"GENCODE+\",\"GENCODE\",\"BIGTranscriptome\",\"RefSeq\",\"FANTOM CAT\",\"MiTranscriptome\",\"NONCODE\", \"Protein coding\", \"CLS FL\")) #
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


# -----------------------
# CAGE strictsTSS  <------------# brakuje mi TSS.strict_hg38.bed
# -----------------------
#hg38_fair+new_CAGE_peaks_phase1and2 powinno byc wymienne ale sprawdzic

# Pliki wynikowe są znacząco inne dla
# TSS.strict_hg38.bed (OLD),
# jak i
# hg38_fair+new_CAGE_peaks_phase1and2.bed (NEW)
# Zrobiłem analizę dla obu, najwyżej się wywali


# OLD
while read lab
do
    echo $lab
    while read end dist
    do
        cat ./Data/Processed/Extracted_ends/$lab.$end.bed | ./Utils/sortbed | ./Utils/jlagarde/bin.bkp/bedtools2/bin/bedtools slop -s -l 50 -r 50 -i stdin -g ./Data/Source/chromInfo_hg38.txt | ./Utils/jlagarde/bin.bkp/bedtools2/bin/bedtools intersect -u -s -a stdin -b ./Data/Source/TSS.strict_hg38.bed > ./Data/Processed/Cage/$lab.$end.bed.vsCage.strictsTSS_OLD.bedtsv
    done < ./Data/Source/ends.dist.tsv
done < ./Data/Source/samples.tsv

# NEW
while read lab
do
    echo $lab
    while read end dist
    do
        cat ./Data/Processed/Extracted_ends/$lab.$end.bed | ./Utils/sortbed | ./Utils/jlagarde/bin.bkp/bedtools2/bin/bedtools slop -s -l 50 -r 50 -i stdin -g ./Data/Source/chromInfo_hg38.txt | ./Utils/jlagarde/bin.bkp/bedtools2/bin/bedtools intersect -u -s -a stdin -b ./Data/Source/hg38_fair+new_CAGE_peaks_phase1and2.bed > ./Data/Processed/Cage/$lab.$end.bed.vsCage.strictsTSS_NEW.bedtsv
    done < ./Data/Source/ends.dist.tsv
done < ./Data/Source/samples.tsv


# --------------------
# stats and plot <------- OLD
# --------------------

while read lab
do
    for end in 5
    do 
        close=`cat ./Data/Processed/Cage/$lab.$end.bed.vsCage.strictsTSS_OLD.bedtsv | cut -f4 | sort | uniq | wc -l`
        total=`cat ./Data/Processed/Extracted_ends/$lab.$end.bed | cut -f4 | sort | uniq | wc -l`
        echo -e "$lab\t$total\t$close\t$end" | awk '{print $0"\t"$3/$2}'    
    done 
done < ./Data/Source/samples.tsv | sed 's/gen/GENCODE/g' | sed 's/refseq/RefSeq/g' | sed 's/fantomCat/FANTOM CAT/g' | sed 's/mitrans/MiTranscriptome/g' | sed 's/bigtrans/BIGTranscriptome/g' | sed 's/noncode/NONCODE/g' | sed 's/cls/GENCODE+/g' | sed 's/pcConf/Protein coding/g' | sed 's/GENCODE+FL/CLS FL/g' > ./Data/Processed/Plots_input/annots.vsCage.strictsTSS.stats_OLD.tsv 

echo "
library(ggplot2)
library(scales)
cbPalette <- c(\"#a6761d\", \"#e6ab02\", \"#FF7F00\", \"#984EA3\", \"#4DAF4A\", \"#377EB8\", \"#E41A1C\",\"#999999\", \"#e7298a\") #
plot <- read.table(\"./Data/Processed/Plots_input/annots.vsCage.strictsTSS.stats_OLD.tsv\", header=F, as.is=T, sep=\"\t\")
colnames(plot)<-c(\"gene\", \"total\", \"count\", \"end\",\"prop\")
plot\$gene=factor(plot\$gene, levels=c(\"GENCODE+\",\"GENCODE\",\"BIGTranscriptome\",\"RefSeq\",\"FANTOM CAT\",\"MiTranscriptome\",\"NONCODE\", \"Protein coding\", \"CLS FL\")) #
pdf(\"./Plots/human.annots.vsCage.strictsTSS.stats_OLD.pdf\", width=9, height=6)
ggplot(data=plot, aes(x=gene, y=prop)) + geom_bar(stat=\"identity\", fill=cbPalette) +
ylab(\"% CAGE(+) 5\' ends\") +
theme_bw(base_size = 28) +
xlab(\"\") +
scale_y_continuous(labels=percent, lim=c(0,1)) +
coord_flip() +
geom_text(position = \"stack\", aes(x = gene, y = prop, label = comma(count), hjust = -0.1, vjust = 0.5), size=7) +
theme(axis.text.x = element_text(angle=45, vjust=0.5)) +
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


# --------------------
# stats and plot <------- NEW
# --------------------

while read lab
do
    for end in 5
    do 
        close=`cat ./Data/Processed/Cage/$lab.$end.bed.vsCage.strictsTSS_NEW.bedtsv | cut -f4 | sort | uniq | wc -l`
        total=`cat ./Data/Processed/Extracted_ends/$lab.$end.bed | cut -f4 | sort | uniq | wc -l`
        echo -e "$lab\t$total\t$close\t$end" | awk '{print $0"\t"$3/$2}'    
    done 
done < ./Data/Source/samples.tsv | sed 's/gen/GENCODE/g' | sed 's/refseq/RefSeq/g' | sed 's/fantomCat/FANTOM CAT/g' | sed 's/mitrans/MiTranscriptome/g' | sed 's/bigtrans/BIGTranscriptome/g' | sed 's/noncode/NONCODE/g' | sed 's/cls/GENCODE+/g' | sed 's/pcConf/Protein coding/g' | sed 's/GENCODE+FL/CLS FL/g' > ./Data/Processed/Plots_input/annots.vsCage.strictsTSS.stats_NEW.tsv 

echo "
library(ggplot2)
library(scales)
cbPalette <- c(\"#a6761d\", \"#e6ab02\", \"#FF7F00\", \"#984EA3\", \"#4DAF4A\", \"#377EB8\", \"#E41A1C\", \"#999999\", \"#e7298a\") #
plot <- read.table(\"./Data/Processed/Plots_input/annots.vsCage.strictsTSS.stats_NEW.tsv\", header=F, as.is=T, sep=\"\t\")
colnames(plot)<-c(\"gene\", \"total\", \"count\", \"end\",\"prop\")
plot\$gene=factor(plot\$gene, levels=c(\"GENCODE+\",\"GENCODE\",\"BIGTranscriptome\",\"RefSeq\",\"FANTOM CAT\",\"MiTranscriptome\",\"NONCODE\", \"Protein coding\", \"CLS FL\")) #
pdf(\"./Plots/human.annots.vsCage.strictsTSS.stats_NEW.pdf\", width=9, height=6)
ggplot(data=plot, aes(x=gene, y=prop)) + geom_bar(stat=\"identity\", fill=cbPalette) +
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


# ------------------------------------------------------------
# % Completness
# ------------------------------------------------------------

# I've downloaded lacking files ($lab.geneIds, $lab.txIds, $lab.buildLoci.geneIds, $lab.buildLoci.txIds)
# Using:
# wget -r -l1 -nd --no-parent -A "*.txIds" https://public-docs.crg.es/rguigo/lncRNA_review/TM/annot.comp/
# wget -r -l1 -nd --no-parent -A "*.geneIds" https://public-docs.crg.es/rguigo/lncRNA_review/TM/annot.comp/


while read lab
do
    cat ./Data/Processed/Cage/$lab.5.bed.vsCage.fantom.bedtsv | cut -f4 | sed 's/,/\n/g' | sort | uniq > ./Data/Processed/Cage/$lab.tx.cage
    cat ./Data/Processed/PolyA/$lab.3.bed.vspolyAsignals.bedtsv | cut -f4 | sed 's/,/\n/g' | sort | uniq > ./Data/Processed/PolyA/$lab.tx.polyA
    #fgrep -f $lab.tx.cage $lab.tx.polyA > $lab.fl.tx.tsv
    ./Utils/join.py -a ./Data/Processed/Cage/$lab.tx.cage -b ./Data/Processed/PolyA/$lab.tx.polyA -x 1 -y 1 > ./Data/Processed/Fl/$lab.fl.tx.tsv
done < ./Data/Source/samples.tsv


# Brakuje mi
# /users/rg/buszczynska/Projects/review/human/annot.comp/$lab.geneIds
# /users/rg/buszczynska/Projects/review/human/annot.comp/$lab.txIds

while read lab
do
    fl=`cat ./Data/Processed/Fl/$lab.fl.tx.tsv | wc -l`
    total=`cat ./Data/Catalogues/$lab.hg38.bed12 | cut -f4 | sort | uniq | wc -l`
    propFl=`echo $fl / $total | bc -l`  
    gene=`cat ./Data/Source/Completness/$lab.geneIds | sort | uniq | wc -l`
    tx=`cat ./Data/Source/Completness/$lab.txIds | sort | uniq | wc -l`
    nbTx=`echo $tx / $gene | bc -l`
    echo -e "$lab\t$propFl\t$gene\t$nbTx"  
done < ./Data/Source/samples.tsv | sed 's/gen/GENCODE/g' | sed 's/refseq/RefSeq/g' | sed 's/fantomCat/FANTOM CAT/g' | sed 's/mitrans/MiTranscriptome/g' | sed 's/bigtrans/BIGTranscriptome/g' | sed 's/noncode/NONCODE/g' | sed 's/cls/GENCODE+/g' | sed 's/pcConf/Protein coding/g' | sed 's/GENCODE+FL/CLS FL/g' > ./Data/Processed/Plots_input/annots.completeness.tsv


# Oraz:
# /users/rg/buszczynska/Projects/review/human/annot.comp/$lab.buildLoci.geneIds
# /users/rg/buszczynska/Projects/review/human/annot.comp/$lab.buildLoci.txIds


while read lab
do
    fl=`cat ./Data/Processed/Fl/$lab.fl.tx.tsv| wc -l`
    total=`cat ./Data/Catalogues/$lab.hg38.bed12| cut -f4| sort| uniq| wc -l`
    propFl=`echo $fl / $total| bc -l`  
    gene=`cat ./Data/Source/Completness/$lab.buildLoci.geneIds | sort | uniq | wc -l`
    tx=`cat ./Data/Source/Completness/$lab.buildLoci.txIds | sort | uniq | wc -l`
    nbTx=`echo $tx / $gene| bc -l`
    echo -e "$lab\t$propFl\t$gene\t$nbTx"  
done < ./Data/Source/samples.tsv | sed 's/gen/GENCODE/g' | sed 's/refseq/RefSeq/g' | sed 's/fantomCat/FANTOM CAT/g' | sed 's/mitrans/MiTranscriptome/g' | sed 's/bigtrans/BIGTranscriptome/g' | sed 's/noncode/NONCODE/g' | sed 's/cls/GENCODE+/g' | sed 's/pcConf/Protein coding/g' | sed 's/GENCODE+FL/CLS FL/g' > ./Data/Processed/Plots_input/annots.completeness.buildLoci.tsv


echo "
library(ggplot2)
library(scales)
cbPalette <- c(\"#a6761d\", \"#e6ab02\", \"#FF7F00\", \"#984EA3\", \"#4DAF4A\", \"#377EB8\", \"#E41A1C\", \"#999999\", \"#e7298a\") #
plot <- read.table(\"./Data/Processed/Plots_input/annots.completeness.buildLoci.tsv\", header=F, as.is=T, sep=\"\t\")
colnames(plot)<-c(\"Catalog\", \"fl\", \"gene\", \"Isoforms\")
plot\$Catalog=factor(plot\$Catalog, levels=c(\"GENCODE+\",\"GENCODE\",\"BIGTranscriptome\",\"RefSeq\",\"FANTOM CAT\",\"MiTranscriptome\",\"NONCODE\", \"Protein coding\", \"CLS FL\")) #
#png(\"annots.completeness.buildLoci.png\", bg = \"white\", units=\"in\", width=8, height=6, res=300)
pdf(\"./Plots/annots.completeness.buildLoci.pdf\", bg = \"white\", width=8, height=7)
#setEPS()
#postscript(\"annots.completeness.buildLoci.eps\", family=\"serif\", width=12, height=8) 
ggplot(data=plot, aes(x=gene, y=fl, size=Isoforms, fill=Catalog)) +
geom_point(shape = 21) +
ylab(\"% Completeness\") +
scale_fill_manual(values = cbPalette) +
theme_bw(base_size = 24) +
xlab(\"Number of loci\") +
scale_x_continuous(lim=c(0,100000), breaks=seq(0,100000, by=25000)) +
scale_y_continuous(labels=percent, lim=c(0,1)) +
theme(axis.text.x = element_text(vjust=0.5),legend.text = element_text(size=12.5)) +
theme(axis.line.x = element_line(colour = \"black\"),
axis.line.y = element_line(colour = \"black\"),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
panel.border = element_blank(),
panel.background = element_blank(),
strip.background = element_rect(colour=\"black\",fill=\"white\")) +
scale_size_area(max_size = 15) +
guides(fill = guide_legend(override.aes = list(size=6.5)))
#+facet_wrap(~ end) 
dev.off()
" | R --slave


# NO PC

cat ./Data/Processed/Plots_input/annots.completeness.buildLoci.tsv | grep -v Protein | grep -v CLS > ./Data/Processed/Plots_input/annots.completeness.buildLoci.noPC.tsv

echo "
library(ggplot2)
library(scales)
cbPalette <- c(\"#a6761d\", \"#e6ab02\", \"#FF7F00\", \"#984EA3\", \"#4DAF4A\", \"#377EB8\", \"#E41A1C\")
plot <- read.table(\"./Data/Processed/Plots_input/annots.completeness.buildLoci.noPC.tsv\", header=F, as.is=T, sep=\"\t\")
colnames(plot)<-c(\"Catalog\", \"fl\", \"gene\", \"Isoforms\")
plot\$Catalog=factor(plot\$Catalog, levels=c(\"GENCODE+\",\"GENCODE\",\"BIGTranscriptome\",\"RefSeq\",\"FANTOM CAT\",\"MiTranscriptome\",\"NONCODE\"))
#png(\"annots.completeness.buildLoci.noPC.png\", bg = \"white\", units=\"in\", width=8, height=6, res=300)
pdf(\"./Plots/annots.completeness.buildLoci.noPC.pdf\", bg = \"white\", width=8, height=7)
#setEPS()
#postscript(\"annots.completeness.buildLoci.noPC.eps\", family=\"serif\", width=12, height=8) 
ggplot(data=plot, aes(x=gene, y=fl, size=Isoforms, fill=Catalog)) +
geom_point(shape = 21) +
ylab(\"% Completeness\") +
scale_fill_manual(values = cbPalette)+ theme_bw(base_size = 24) +
xlab(\"Number of loci\") +
scale_x_continuous(lim=c(0,100000), breaks=seq(0,100000, by=25000)) +
scale_y_continuous(labels=percent, lim=c(0,0.3)) +
theme(axis.text.x = element_text(vjust=0.5),legend.text = element_text(size=12.5)) +
theme(axis.line.x = element_line(colour = \"black\"),
axis.line.y = element_line(colour = \"black\"),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
panel.border = element_blank(),
panel.background = element_blank(),
strip.background = element_rect(colour=\"black\",fill=\"white\")) +
scale_size_area(max_size = 17) +
guides(fill = guide_legend(override.aes = list(size=6.5)))
#+facet_wrap(~ end) 
dev.off()
" | R --slave



