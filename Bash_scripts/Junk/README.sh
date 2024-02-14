working=/users/rg/buszczynska/Projects/review/human/cage.tss

cd $working

# ---------------------------------------------------------------
# link BED12 files
# ---------------------------------------------------------------

cp /users/rg/buszczynska/Projects/review/human/annot.comp/samples.tsv .

# CAGE 

wget http://fantom.gsc.riken.jp/5/datafiles/latest/extra/CAGE_peaks/hg19.cage_peak_phase1and2combined_coord.bed.gz .
gunzip hg19.cage_peak_phase1and2combined_coord.bed.gz

liftOver hg19.cage_peak_phase1and2combined_coord.bed /users/rg/buszczynska/liftOver/hg19ToHg38.over.chain hg38.cage_peak_phase1and2combined_coord.bed hg38.cage_peak_phase1and2combined_coord.unMapped

# Get the full legth CLS

#wget https://public_docs.crg.es/rguigo/CLS/data/mergedTMs/non-anchored/gtf/hsAll_Cap1_all_noAnchor.compmerge.cage+polyASupported.gtf.gz . 
#gunzip hsAll_Cap1_all_noAnchor.compmerge.cage+polyASupported.gtf.gz 

#~abreschi/Documents/utils/extract.gtf.tags.sh hsAll_Cap1_all_noAnchor.compmerge.cage+polyASupported.gtf transcript_id | sort| uniq > clsFL.txIds
#fgrep -f clsFL.txIds cls.hg38.bed12 > clsFL.hg38.bed12 

#join.py -a clsFL.txIds -b cls.hg38.bed12 -x 1 -y 4 > clsFL.hg38.bed12

#cat /users/rg/jlagarde/projects/encode/scaling/whole_genome/lncRnaCapture/pacBioVsHiSeq/stringtieVsPacBio/transcripts/5pEnds_support/hs.CLS-noAnchor.onTargetLncRna.5pEndsClusters.cage+polyASupported.bed|cut -f4|sed 's/,/\n/g'| sort|uniq > clsFL.txIds

# -------------------
# link bed12
# -------------------
while read lab
do
    cat /users/rg/buszczynska/Projects/review/human/annot.comp/$lab.hg38.bed12 | awk 'BEGIN{FS=OFS="\t"}$6!="."' > $lab.hg38.bed12
done < samples.tsv

# -------------------
# extract ends
# -------------------

for file in `ls *.hg38.bed12`
do
    echo $file
    lab=`basename $file| awk -F "." '{print $1}'`
    while read end dist
    do
        echo $end
        cat $file | /users/rg/jlagarde/julien_utils/extractTranscriptEndsFromBed12.pl $end |sortbed | ~jlagarde/bin/bedtools2/bin/bedtools merge -s -d $dist -c 4 -o collapse -i stdin | awk '{print $1"\t"$2"\t"$3"\t"$5"\t0\t"$4}' > $lab.$end.bed
    done < ends.dist.tsv
done

# -------------------
# polyA
# -------------------

while read lab
do
    echo $lab
    while read end dist
    do
        cat $lab.$end.bed | sortbed | ~jlagarde/bin/bedtools2/bin/bedtools slop -s -l 50 -r -10 -i stdin -g ~jlagarde/hg38.genome | ~jlagarde/bin/bedtools2/bin/bedtools intersect -u -s -a stdin -b /users/rg/jlagarde/projects/polyAsignals/hg38.polyAsignals.bed > $lab.$end.bed.vspolyAsignals.bedtsv
    done < ends.dist.tsv
done < samples.tsv

# --------------------
# stats and plot
# --------------------

while read lab
do
    for end in 3
    do 
        close=`cat $lab.$end.bed.vspolyAsignals.bedtsv| cut -f4| sort|uniq| wc -l`
        total=`cat $lab.$end.bed| cut -f4| sort|uniq| wc -l`
        echo -e "$lab\t$total\t$close\t$end" | awk '{print $0"\t"$3/$2}'    
    done 
done < samples.tsv | sed 's/gen/GENCODE/g'| sed 's/refseq/RefSeq/g'|sed 's/fantomCat/FANTOM CAT/g'| sed 's/mitrans/MiTranscriptome/g'| sed 's/bigtrans/BIGTranscriptome/g'| sed 's/noncode/NONCODE/g'| sed 's/cls/GENCODE+/g'| sed 's/pcConf/Protein coding/g'| sed 's/GENCODE+FL/CLS FL/g' > annots.polyAsignals.stats.tsv 


echo "
library(ggplot2)
library(scales)
cbPalette <- c(\"#a6761d\", \"#e6ab02\", \"#FF7F00\", \"#984EA3\", \"#4DAF4A\", \"#377EB8\", \"#E41A1C\",\"#999999\", \"#e7298a\")
plot <- read.table(\"annots.polyAsignals.stats.tsv\", header=F, as.is=T, sep=\"\t\")
colnames(plot)<-c(\"gene\", \"total\", \"count\", \"end\",\"prop\")
plot\$gene=factor(plot\$gene, levels=c(\"GENCODE+\",\"GENCODE\",\"BIGTranscriptome\",\"RefSeq\",\"FANTOM CAT\",\"MiTranscriptome\",\"NONCODE\", \"Protein coding\", \"CLS FL\"))
png(\"human.annots.polyAsignals.stats.png\", bg = \"white\", units=\"in\", width=8, height=6, res=300)
#setEPS()
#postscript(\"human.annots.polyAsignals.stats.eps\", family=\"serif\", width=8, height=10)
ggplot(data=plot, aes(x=gene, y=prop)) + geom_bar(stat=\"identity\", fill=cbPalette) + ylab(\"% PAS(+) 3\' ends\")+ theme_bw(base_size = 28) + xlab(\"\") + scale_y_continuous(labels=percent, limits=c(0, 1)) + coord_flip()+ geom_text(position = \"stack\", aes(x = gene, y = prop, label = comma(count), hjust = 1, vjust = 0.5), size=7)+ theme(axis.text.x  = element_text(angle=45, vjust=0.5), axis.ticks.y = element_blank(), axis.text.y  = element_text(angle=0, hjust=0.5))+ theme(axis.line.x = element_line(colour = \"black\"), axis.line.y = element_line(colour = \"black\"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(), strip.background = element_rect(colour=\"black\",fill=\"white\")) 
#+facet_wrap(~ end) 
dev.off()
" | R --slave

cp human.annots.polyAsignals.stats.png /users/rg/buszczynska/public_html/gencode/tmp/B/review
#cp human.annots.polyAsignals.stats.eps /users/rg/buszczynska/public_html/gencode/tmp/B/review

# -----------------------
# CAGE FANTOM pahse 1+2
# -----------------------

while read lab
do
    echo $lab
    while read end dist
    do
        cat $lab.$end.bed | sortbed | ~jlagarde/bin/bedtools2/bin/bedtools slop -s -l 50 -r 50 -i stdin -g ~jlagarde/hg38.genome | ~jlagarde/bin/bedtools2/bin/bedtools intersect -u -s -a stdin -b hg38.cage_peak_phase1and2combined_coord.bed > $lab.$end.bed.vsCage.fantom.bedtsv
    done < ends.dist.tsv
done < samples.tsv

# --------------------
# stats and plot
# --------------------

while read lab
do
    for end in 5
    do 
        close=`cat $lab.$end.bed.vsCage.fantom.bedtsv| cut -f4| sort|uniq| wc -l`
        total=`cat $lab.$end.bed| cut -f4| sort|uniq| wc -l`
        echo -e "$lab\t$total\t$close\t$end" | awk '{print $0"\t"$3/$2}'    
    done 
done < samples.tsv | sed 's/gen/GENCODE/g'| sed 's/refseq/RefSeq/g'|sed 's/fantomCat/FANTOM CAT/g'| sed 's/mitrans/MiTranscriptome/g'| sed 's/bigtrans/BIGTranscriptome/g'| sed 's/noncode/NONCODE/g'| sed 's/cls/GENCODE+/g'| sed 's/pcConf/Protein coding/g'| sed 's/GENCODE+FL/CLS FL/g'  > annots.vsCage.fantom.stats.tsv 

echo "
library(ggplot2)
library(scales)
cbPalette <- c(\"#a6761d\", \"#e6ab02\", \"#FF7F00\", \"#984EA3\", \"#4DAF4A\", \"#377EB8\", \"#E41A1C\",\"#999999\", \"#e7298a\")
plot <- read.table(\"annots.vsCage.fantom.stats.tsv\", header=F, as.is=T, sep=\"\t\")
colnames(plot)<-c(\"gene\", \"total\", \"count\", \"end\",\"prop\")
plot\$gene=factor(plot\$gene, levels=c(\"GENCODE+\",\"GENCODE\",\"BIGTranscriptome\",\"RefSeq\",\"FANTOM CAT\",\"MiTranscriptome\",\"NONCODE\", \"Protein coding\", \"CLS FL\"))
png(\"human.annots.vsCage.fantom.stats.png\", bg = \"white\", units=\"in\", width=6, height=6, res=300)
#setEPS()
#postscript(\"human.annots.vsCage.fantom.stats.eps\", family=\"serif\", width=6, height=10)
ggplot(data=plot, aes(x=gene, y=prop)) + geom_bar(stat=\"identity\", fill=cbPalette) + ylab(\"% CAGE(+) 5\' ends\")+ theme_bw(base_size = 28) + xlab(\"\")+ scale_y_reverse(labels=percent, lim=c(1,0)) + coord_flip()+ geom_text(position = \"stack\", aes(x = gene, y = prop, label = comma(count), hjust = 1, vjust = 0.5), size=7) +theme(axis.text.x  = element_text(angle=45, vjust=0.5), axis.ticks.y = element_blank(), axis.text.y  = element_blank())+ theme(axis.line.x = element_line(colour = \"black\"), axis.line.y = element_line(colour = \"black\"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(), strip.background = element_rect(colour=\"black\",fill=\"white\"))
#+facet_wrap(~ end) 
dev.off()
" | R --slave

#cp human.annots.vsCage.fantom.stats.eps /users/rg/buszczynska/public_html/gencode/tmp/B/review
cp human.annots.vsCage.fantom.stats.png /users/rg/buszczynska/public_html/gencode/tmp/B/review

#montage -geometry 100% -tile 2x1 human.annots.vsCage.fantom.stats.png human.annots.polyAsignals.stats.png human.annots.cageFantom+polyAsignals.stats.png
#cp human.annots.cageFantom+polyAsignals.stats.png /users/rg/buszczynska/public_html/gencode/tmp/B/review

# -----------------------
# CAGE strictsTSS
# -----------------------

while read lab
do
    echo $lab
    while read end dist
    do
        cat $lab.$end.bed | sortbed | ~jlagarde/bin/bedtools2/bin/bedtools slop -s -l 50 -r 50 -i stdin -g ~jlagarde/hg38.genome | ~jlagarde/bin/bedtools2/bin/bedtools intersect -u -s -a stdin -b /no_backup/rg/jlagarde/projects/fantom5/TSS_classifier/TSS.strict_hg38.bed > $lab.$end.bed.vsCage.strictsTSS.bedtsv
    done < ends.dist.tsv
done < samples.tsv

# --------------------
# stats and plot
# --------------------

while read lab
do
    for end in 5
    do 
        close=`cat $lab.$end.bed.vsCage.strictsTSS.bedtsv| cut -f4| sort|uniq| wc -l`
        total=`cat $lab.$end.bed| cut -f4| sort|uniq| wc -l`
        echo -e "$lab\t$total\t$close\t$end" | awk '{print $0"\t"$3/$2}'    
    done 
done < samples.tsv | sed 's/gen/GENCODE/g'| sed 's/refseq/RefSeq/g'|sed 's/fantomCat/FANTOM CAT/g'| sed 's/mitrans/MiTranscriptome/g'| sed 's/bigtrans/BIGTranscriptome/g'| sed 's/noncode/NONCODE/g'| sed 's/cls/GENCODE+/g'| sed 's/pc/Protein coding/g'| sed 's/GENCODE+FL/CLS FL/g' > annots.vsCage.strictsTSS.stats.tsv 

echo "
library(ggplot2)
library(scales)
cbPalette <- c(\"#a6761d\", \"#e6ab02\", \"#FF7F00\", \"#984EA3\", \"#4DAF4A\", \"#377EB8\", \"#E41A1C\",\"#999999\", \"#e7298a\")
plot <- read.table(\"annots.vsCage.strictsTSS.stats.tsv\", header=F, as.is=T, sep=\"\t\")
colnames(plot)<-c(\"gene\", \"total\", \"count\", \"end\",\"prop\")
plot\$gene=factor(plot\$gene, levels=c(\"GENCODE+\",\"GENCODE\",\"BIGTranscriptome\",\"RefSeq\",\"FANTOM CAT\",\"MiTranscriptome\",\"NONCODE\", \"Protein coding\", \"CLS FL\"))
png(\"human.annots.vsCage.strictsTSS.stats.png\", bg = \"white\", units=\"in\", width=8, height=6, res=300)
ggplot(data=plot, aes(x=gene, y=prop)) + geom_bar(stat=\"identity\", fill=cbPalette) + ylab(\"% CAGE(+) 5\' ends\")+ theme_bw(base_size = 28) + xlab(\"\")+ scale_y_reverse(labels=percent, lim=c(1,0)) + coord_flip()+ geom_text(position = \"stack\", aes(x = gene, y = prop, label = comma(count), hjust = 1, vjust = 0.5), size=7) +theme(axis.text.x  = element_text(angle=45, vjust=0.5), axis.ticks.y = element_blank(), axis.text.y  = element_blank())+ theme(axis.line.x = element_line(colour = \"black\"), axis.line.y = element_line(colour = \"black\"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(), strip.background = element_rect(colour=\"black\",fill=\"white\"))
#+facet_wrap(~ end) 
dev.off()
" | R --slave

cp human.annots.vsCage.strictsTSS.stats.png /users/rg/buszczynska/public_html/gencode/tmp/B/review

montage -geometry 100% -tile 2x1 human.annots.vsCage.strictsTSS.stats.png human.annots.polyAsignals.stats.png human.annots.cageStrictsTSS+polyAsignals.stats.png
cp human.annots.cageStrictsTSS+polyAsignals.stats.png /users/rg/buszczynska/public_html/gencode/tmp/B/review

# ------------------------------------------------------------
# % Completness
# ------------------------------------------------------------

while read lab
do
    cat $lab.5.bed.vsCage.fantom.bedtsv | cut -f4| sed 's/,/\n/g'| sort|uniq > $lab.tx.cage
    cat $lab.3.bed.vspolyAsignals.bedtsv| cut -f4| sed 's/,/\n/g'| sort|uniq > $lab.tx.polyA
    #fgrep -f $lab.tx.cage $lab.tx.polyA > $lab.fl.tx.tsv
    join.py -a $lab.tx.cage -b $lab.tx.polyA -x 1 -y 1 > $lab.fl.tx.tsv
done < samples.tsv


while read lab
do
    fl=`cat $lab.fl.tx.tsv| wc -l`
    total=`cat $lab.hg38.bed12| cut -f4| sort| uniq| wc -l`
    propFl=`echo $fl / $total| bc -l`  
    gene=`cat /users/rg/buszczynska/Projects/review/human/annot.comp/$lab.geneIds|sort|uniq| wc -l`
    tx=`cat /users/rg/buszczynska/Projects/review/human/annot.comp/$lab.txIds|sort|uniq| wc -l`
    nbTx=`echo $tx / $gene| bc -l`
    echo -e "$lab\t$propFl\t$gene\t$nbTx"  
done < samples.tsv | sed 's/gen/GENCODE/g'| sed 's/refseq/RefSeq/g'|sed 's/fantomCat/FANTOM CAT/g'| sed 's/mitrans/MiTranscriptome/g'| sed 's/bigtrans/BIGTranscriptome/g'| sed 's/noncode/NONCODE/g'| sed 's/cls/GENCODE+/g'| sed 's/pcConf/Protein coding/g'| sed 's/GENCODE+FL/CLS FL/g' > annots.completeness.tsv

while read lab
do
    fl=`cat $lab.fl.tx.tsv| wc -l`
    total=`cat $lab.hg38.bed12| cut -f4| sort| uniq| wc -l`
    propFl=`echo $fl / $total| bc -l`  
    gene=`cat /users/rg/buszczynska/Projects/review/human/annot.comp/$lab.buildLoci.geneIds|sort|uniq| wc -l`
    tx=`cat /users/rg/buszczynska/Projects/review/human/annot.comp/$lab.buildLoci.txIds|sort|uniq| wc -l`
    nbTx=`echo $tx / $gene| bc -l`
    echo -e "$lab\t$propFl\t$gene\t$nbTx"  
done < samples.tsv | sed 's/gen/GENCODE/g'| sed 's/refseq/RefSeq/g'|sed 's/fantomCat/FANTOM CAT/g'| sed 's/mitrans/MiTranscriptome/g'| sed 's/bigtrans/BIGTranscriptome/g'| sed 's/noncode/NONCODE/g'| sed 's/cls/GENCODE+/g'| sed 's/pcConf/Protein coding/g'| sed 's/GENCODE+FL/CLS FL/g' > annots.completeness.buildLoci.tsv

# ------------------------------------------
# New with intropolis support
# ------------------------------------------

while read lab
do
    fl=`cat $lab.fl.tx.tsv| wc -l`
    total=`cat $lab.hg38.bed12| cut -f4| sort| uniq| wc -l`
    propFl=`echo $fl / $total| bc -l`  
    gene=`cat /users/rg/buszczynska/Projects/review/human/sjs/$lab.chrOnly.buildLoci.geneIds|sort|uniq| wc -l`
    total=`cat /users/rg/buszczynska/Projects/review/human/sjs/$lab.tx.supp.10smpl| cut -f1| sort| uniq| wc -l`
    validated=`cat /users/rg/buszczynska/Projects/review/human/sjs/$lab.tx.supp.10smpl| awk '$2=="VALIDATED"'| cut -f1| sort| uniq| wc -l`
    valP=`echo "$validated *100" / $total| bc  -l`
    echo -e "$lab\t$propFl\t$gene\t$valP"  
done < samples.tsv | sed 's/gen/GENCODE/g'| sed 's/refseq/RefSeq/g'|sed 's/fantomCat/FANTOM CAT/g'| sed 's/mitrans/MiTranscriptome/g'| sed 's/bigtrans/BIGTranscriptome/g'| sed 's/noncode/NONCODE/g'| sed 's/cls/GENCODE+/g'| sed 's/pcConf/Protein coding/g'| sed 's/GENCODE+FL/CLS FL/g' > annots.completeness.intropolis.buildLoci.tsv


echo "   <------tego nie robic
library(ggplot2)
library(scales)
cbPalette <- c(\"#a6761d\", \"#e6ab02\", \"#FF7F00\", \"#984EA3\", \"#4DAF4A\", \"#377EB8\", \"#E41A1C\",\"#999999\", \"#e7298a\")
plot <- read.table(\"annots.completeness.intropolis.buildLoci.tsv\", header=F, as.is=T, sep=\"\t\")
colnames(plot)<-c(\"Catalog\", \"fl\", \"gene\", \"Isoforms\")
plot\$Isoforms=round(plot\$Isoforms,2)
plot\$Catalog=factor(plot\$Catalog, levels=c(\"GENCODE+\",\"GENCODE\",\"BIGTranscriptome\",\"RefSeq\",\"FANTOM CAT\",\"MiTranscriptome\",\"NONCODE\", \"Protein coding\", \"CLS FL\"))
pdf(\"annots.completeness.intropolis.buildLoci.pdf\", bg = \"white\", width=8, height=6.2)
ggplot(data=plot, aes(x=gene, y=fl, size=Isoforms, fill=Catalog)) + geom_point(shape = 21) + ylab(\"% Completeness\")+ scale_fill_manual(values = cbPalette)+ theme_bw(base_size = 24) + xlab(\"Number of loci\")+ scale_x_continuous(lim=c(0,100000), breaks=seq(0,100000, by=25000))+ scale_y_continuous(labels=percent, lim=c(0,1)) +theme(axis.text.x  = element_text(vjust=0.5),legend.text = element_text(size=12.5))+ theme(axis.line.x = element_line(colour = \"black\"), axis.line.y = element_line(colour = \"black\"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(), strip.background = element_rect(colour=\"black\",fill=\"white\"))+ scale_size_area(max_size = 15)+ guides(fill = guide_legend(override.aes = list(size=6.5)))
#+facet_wrap(~ end) 
dev.off()
" | R --slave

cp annots.completeness.intropolis.buildLoci.pdf /users/rg/buszczynska/public_html/gencode/tmp/B/review

echo " <------ten robic
library(ggplot2)
library(scales)
cbPalette <- c(\"#a6761d\", \"#e6ab02\", \"#FF7F00\", \"#984EA3\", \"#4DAF4A\", \"#377EB8\", \"#E41A1C\",\"#999999\", \"#e7298a\")
plot <- read.table(\"annots.completeness.buildLoci.tsv\", header=F, as.is=T, sep=\"\t\")
colnames(plot)<-c(\"Catalog\", \"fl\", \"gene\", \"Isoforms\")
plot\$Catalog=factor(plot\$Catalog, levels=c(\"GENCODE+\",\"GENCODE\",\"BIGTranscriptome\",\"RefSeq\",\"FANTOM CAT\",\"MiTranscriptome\",\"NONCODE\", \"Protein coding\", \"CLS FL\"))
#png(\"annots.completeness.buildLoci.png\", bg = \"white\", units=\"in\", width=8, height=6, res=300)
pdf(\"annots.completeness.buildLoci.pdf\", bg = \"white\", width=8, height=6.2)
#setEPS()
#postscript(\"annots.completeness.buildLoci.eps\", family=\"serif\", width=12, height=8) 
ggplot(data=plot, aes(x=gene, y=fl, size=Isoforms, fill=Catalog)) + geom_point(shape = 21) + ylab(\"% Completeness\")+ scale_fill_manual(values = cbPalette)+ theme_bw(base_size = 24) + xlab(\"Number of loci\")+ scale_x_continuous(lim=c(0,100000), breaks=seq(0,100000, by=25000))+ scale_y_continuous(labels=percent, lim=c(0,1)) +theme(axis.text.x  = element_text(vjust=0.5),legend.text = element_text(size=12.5))+ theme(axis.line.x = element_line(colour = \"black\"), axis.line.y = element_line(colour = \"black\"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(), strip.background = element_rect(colour=\"black\",fill=\"white\"))+ scale_size_area(max_size = 15)+ guides(fill = guide_legend(override.aes = list(size=6.5)))
#+facet_wrap(~ end) 
dev.off()
" | R --slave

#cp annots.completeness.buildLoci.png /users/rg/buszczynska/public_html/gencode/tmp/B/review
cp annots.completeness.buildLoci.pdf /users/rg/buszczynska/public_html/gencode/tmp/B/review
#cp annots.completeness.buildLoci.eps /users/rg/buszczynska/public_html/gencode/tmp/B/review

# NO PC

cat annots.completeness.buildLoci.tsv| grep -v Protein| grep -v CLS > annots.completeness.buildLoci.noPC.tsv

echo "
library(ggplot2)
library(scales)
cbPalette <- c(\"#a6761d\", \"#e6ab02\", \"#FF7F00\", \"#984EA3\", \"#4DAF4A\", \"#377EB8\", \"#E41A1C\")
plot <- read.table(\"annots.completeness.buildLoci.noPC.tsv\", header=F, as.is=T, sep=\"\t\")
colnames(plot)<-c(\"Catalog\", \"fl\", \"gene\", \"Isoforms\")
plot\$Catalog=factor(plot\$Catalog, levels=c(\"GENCODE+\",\"GENCODE\",\"BIGTranscriptome\",\"RefSeq\",\"FANTOM CAT\",\"MiTranscriptome\",\"NONCODE\"))
#png(\"annots.completeness.buildLoci.noPC.png\", bg = \"white\", units=\"in\", width=8, height=6, res=300)
pdf(\"annots.completeness.buildLoci.noPC.pdf\", bg = \"white\", width=8, height=6.2)
#setEPS()
#postscript(\"annots.completeness.buildLoci.noPC.eps\", family=\"serif\", width=12, height=8) 
ggplot(data=plot, aes(x=gene, y=fl, size=Isoforms, fill=Catalog)) + geom_point(shape = 21) + ylab(\"% Completeness\")+ scale_fill_manual(values = cbPalette)+ theme_bw(base_size = 24) + xlab(\"Number of loci\")+ scale_x_continuous(lim=c(0,100000), breaks=seq(0,100000, by=25000))+ scale_y_continuous(labels=percent, lim=c(0,0.3)) +theme(axis.text.x  = element_text(vjust=0.5),legend.text = element_text(size=12.5))+ theme(axis.line.x = element_line(colour = \"black\"), axis.line.y = element_line(colour = \"black\"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(), strip.background = element_rect(colour=\"black\",fill=\"white\"))+ scale_size_area(max_size = 17)+ guides(fill = guide_legend(override.aes = list(size=6.5)))
#+facet_wrap(~ end) 
dev.off()
" | R --slave

#cp annots.completeness.buildLoci.noPC.png /users/rg/buszczynska/public_html/gencode/tmp/B/review
cp annots.completeness.buildLoci.noPC.pdf /users/rg/buszczynska/public_html/gencode/tmp/B/review
#cp annots.completeness.buildLoci.noPC.eps /users/rg/buszczynska/public_html/gencode/tmp/B/review

# ---------------------------------------------------
# QC for the CAGE/polyA support
# ---------------------------------------------------

while read lab
do
    gene=`cat /users/rg/buszczynska/Projects/review/human/annot.comp/$lab.geneIds|sort|uniq| wc -l`
    tx=`cat /users/rg/buszczynska/Projects/review/human/annot.comp/$lab.txIds|sort|uniq| wc -l`
    end5=`cat $lab.5.bed|sort|uniq| wc -l`
    end3=`cat $lab.3.bed|sort|uniq| wc -l`
    echo -e "$lab\t$gene\tLoci\n$lab\t$tx\tIsoforms\n$lab\t$end5\t5clusters\n$lab\t$end3\t3clusters" 
done < samples.tsv | sed 's/gen/GENCODE/g'| sed 's/refseq/RefSeq/g'|sed 's/fantomCat/FANTOM_CAT/g'| sed 's/mitrans/MiTranscriptome/g'| sed 's/bigtrans/BIGTranscriptome/g'| sed 's/noncode/NONCODE/g'| sed 's/cls/GENCODE+/g'| sed 's/pc/Protein_coding/g'| sed 's/GENCODE+FL/CLS_FL/g' > annots.qc.stats.tsv


echo "
library(ggplot2)
library(RColorBrewer)
library(scales)
df <- read.table('annots.qc.stats.tsv', header=F, sep=\"\t\")
colnames(df)<-c(\"lab\",\"count\",\"typ\")
df\$lab=factor(df\$lab, levels=c(\"NONCODE\",\"MiTranscriptome\",\"FANTOM_CAT\",\"RefSeq\",\"GENCODE\",\"BIGTranscriptome\",\"GENCODE+\", \"Protein_coding\", \"CLS_FL\"))
df\$typ=factor(df\$typ, levels=c(\"Loci\",\"Isoforms\",\"5clusters\",\"3clusters\"))
png(\"human.annots.qc.stats.png\", bg = \"white\", width = 10, height = 6, res=300, units=\"in\")
ggplot(df, aes(x = lab, y = count, fill = typ))+ geom_bar(stat = \"identity\", color=\"black\", position=position_dodge(), width=0.9)+ theme_bw()+ scale_fill_brewer(palette=\"Dark2\")+ theme(legend.position=\"top\")+ xlab(\"\") + ylab (\"# Features\")+theme(axis.text.x = element_text(size = 16, colour = \"black\", angle=90), axis.text.y = element_text(size = 16,  colour = \"black\"), legend.text = element_text(size=18), plot.title = element_text(size = 22), legend.title =element_text(size=0, color=\"white\"), axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18), legend.text.align=0.5, strip.text.x=element_text(size=18), strip.text.y=element_text(size=18)) + theme(legend.position=\"top\") + ggtitle('')+ theme(plot.title = element_text(hjust = 0.5))+ geom_text(aes(x = lab, y=count, label = comma(round(count,1))),position = position_dodge(width = 0.9), size=3.2, angle=90,hjust = 1)
dev.off()
" | R --slave

cp human.annots.qc.stats.png /users/rg/buszczynska/public_html/gencode/tmp/B/review

# --------------------------------------------------
# FANTOM CAT TIEScore
# --------------------------------------------------

# It's to see what's the distribution of TIEScore for CAGE/polyA supported fantomCat txs

# FL TX
#fgrep -f fantomCat.fl.tx.tsv /users/rg/buszczynska/Projects/review/human/annot.comp/fantomCat.lncRNAs.gtf | awk '$3=="transcript"' > fantomCat.lncRNAs.fl.tx.gtf
fgrep -f fantomCat.tx.cage /users/rg/buszczynska/Projects/review/human/annot.comp/fantomCat.lncRNAs.gtf | awk '$3=="transcript"' > fantomCat.lncRNAs.cage.tx.gtf

# NOT FL TX
#fgrep -vf fantomCat.fl.tx.tsv /users/rg/buszczynska/Projects/review/human/annot.comp/fantomCat.lncRNAs.gtf | awk '$3=="transcript"' > fantomCat.lncRNAs.nofl.tx.gtf
fgrep -vf fantomCat.tx.cage /users/rg/buszczynska/Projects/review/human/annot.comp/fantomCat.lncRNAs.gtf | awk '$3=="transcript"' > fantomCat.lncRNAs.nocage.tx.gtf
cat fantomCat.hg38.bed12| cut -f4| fgrep -f - fantomCat.lncRNAs.nocage.tx.gtf > tmp && mv tmp fantomCat.lncRNAs.nocage.tx.gtf

# get the tx and tiescore
#~abreschi/Documents/utils/extract.gtf.tags.sh fantomCat.lncRNAs.fl.tx.gtf transcript_id,TIEScore| awk 'BEGIN{FS=OFS="\t"}{print $1,$2,"Supported"}' > fantomCat.lncRNAs.fl.txID.tiescore.tsv
#~abreschi/Documents/utils/extract.gtf.tags.sh fantomCat.lncRNAs.nofl.tx.gtf transcript_id,TIEScore| awk 'BEGIN{FS=OFS="\t"}{print $1,$2,"Unsupported"}' > fantomCat.lncRNAs.nofl.txID.tiescore.tsv

~abreschi/Documents/utils/extract.gtf.tags.sh fantomCat.lncRNAs.cage.tx.gtf transcript_id,TIEScore| awk 'BEGIN{FS=OFS="\t"}{print $1,$2,"5end_supported"}' > fantomCat.lncRNAs.cage.txID.tiescore.tsv
~abreschi/Documents/utils/extract.gtf.tags.sh fantomCat.lncRNAs.nocage.tx.gtf transcript_id,TIEScore| awk 'BEGIN{FS=OFS="\t"}{print $1,$2,"5end_unsupported"}' > fantomCat.lncRNAs.nocage.txID.tiescore.tsv

# join
#cat fantomCat.lncRNAs.fl.txID.tiescore.tsv fantomCat.lncRNAs.nofl.txID.tiescore.tsv > fantomCat.lncRNAs.no+fl.txID.tiescore.tsv
cat fantomCat.lncRNAs.cage.txID.tiescore.tsv fantomCat.lncRNAs.nocage.txID.tiescore.tsv > fantomCat.lncRNAs.no+cage.txID.tiescore.tsv

# plot

echo "
library(ggplot2)
library(scales)
cbPalette <- c(\"#65bbb6\", \"#bef799\")
plot <- read.table('fantomCat.lncRNAs.no+cage.txID.tiescore.tsv', header=F, as.is=T, sep=\"\t\")
colnames(plot)<-c(\"txID\", \"count\", \"type\")
plot\$type=factor(plot\$type, levels=c(\"5end_supported\",\"5end_unsupported\"))
list=c(\"5end_supported\",\"5end_unsupported\")
label=matrix(0,length(list),2) 
for (i in 1:length(list))
{
median=median(plot[plot\$type==list[i], 2],na.rm=T)
len=length(plot[plot\$type==list[i], 2])
t=paste(\"N =\", prettyNum(len,big.mark=\",\",scientific=FALSE), \"\nMedian =\", prettyNum(round(median,2),big.mark=\",\",scientific=FALSE))
label[i,1]=list[i]
label[i,2]=t
}

label=as.data.frame(label)
row.names(label)=NULL
colnames(label)=c('type', 'text')
#png(\"fantomCat.lncRNAs.no+cage.txID.tiescore.png\", bg = \"white\", units=\"in\", width=6, height=6, res=300)
pdf(\"fantomCat.lncRNAs.no+cage.txID.tiescore.pdf\", bg = \"white\", width=6, height=6)
#setEPS()
#postscript(\"fantomCat.lncRNAs.no+cage.txID.tiescore.eps\", width=5, height=5)
ggplot(plot, aes(x=count, color=type)) + geom_density(size=2.5) + scale_color_manual(values=cbPalette)+ xlab(\"TIEScore\") + scale_x_continuous(limits=c(30, 90), breaks=seq(30,90, by=10))+ theme_bw() + ylab (\"Density\")+theme(axis.text.x = element_text(size = 16, colour = \"black\"), axis.text.y = element_text(size = 16,  colour = \"black\"), legend.text = element_text(size=16), plot.title = element_text(size = 22), legend.title =element_text(size=0, color=\"white\"), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 18), legend.text.align=0.5, strip.text.x=element_text(size=18), strip.text.y=element_text(size=18))+ ggtitle('')+ theme(legend.position=\"top\")+ geom_text(data=label[1,], aes(x=80, y=0.07, label=text), size=6, color=\"#65bbb6\")+ geom_text(data=label[2,], aes(x=40, y=0.07, label=text), size=6, color=\"#bef799\")
dev.off()
" | R --slave

cp fantomCat.lncRNAs.no+cage.txID.tiescore.pdf /users/rg/buszczynska/public_html/gencode/tmp/B/review


