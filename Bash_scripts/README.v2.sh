working=/users/rg/buszczynska/Projects/captrapSeq/reviwers/completness

cd $working

cp /users/rg/buszczynska/Projects/captrapSeq/reviwers/ont.brain.samples.tsv .
cp /users/rg/buszczynska/Projects/captrapSeq/reviwers/ends.dist.tsv .
cp /nfs/users/project/gencode_006070/jlagarde/polyAsignals/hg38.polyAsignals.bed .
cp /nfs/users/project/gencode_006070_no_backup/epalumbo/pilot/annotations/cage/hg38_fair+new_CAGE_peaks_phase1and2.bed .

# -----------------------------
# raw reads
# -----------------------------

# get the read gff to see if they are spliced
while read lab sampl
do
    name=`echo $lab\_HpreCap_0+\_$sampl`
    ln -s /nfs/users/project/gencode_006070_no_backup/epalumbo/pilot/output/mappings/readBedToGff/$name.gff.gz .
    zcat $name.gff.gz > $name.gff  
done < ont.brain.samples.tsv

# get polyA reads
while read lab sampl
do
   name=`echo $lab\_HpreCap_0+\_$sampl`
   cp /nfs/users/project/gencode_006070_no_backup/epalumbo/pilot/output/mappings/getPolyAreadsList/$name.polyAreads.list .
done < ont.brain.samples.tsv 

# get raw reads in bed
while read lab sampl
do
    name=`echo $lab\_HpreCap_0+\_$sampl`
    ln -s /nfs/users/project/gencode_006070_no_backup/epalumbo/pilot/output/mappings/readBamToBed/$name.bed.gz .
    zcat $name.bed.gz > $name.bed 
done < ont.brain.samples.tsv

# -------------------
# extract ends from the reads, polyA and CAGE
# -------------------

while read lab sampl
do
    name=`echo $lab\_HpreCap_0+\_$sampl`
    while read end dist
    do
        echo "cat $name.bed | /users/rg/jlagarde/julien_utils_public/extractTranscriptEndsFromBed12.pl $end |sortbed > $name.$end.bed" > $lab.$sampl.$end.endsExtr.sh
        qsub -cwd -N $lab.$sampl.$end.endsExtr -m abe -M barbara.uszczynska@crg.es -q rg-el7,long-centos79,short-centos79  -l h_rt=48:00:00 -l virtual_free=20G -pe smp 1 -e $lab.$sampl.$end.endsExtr.log $lab.$sampl.$end.endsExtr.sh
    done < ends.dist.tsv
done < ont.brain.samples.tsv

# -------------------
# polyA
# -------------------

while read lab sampl
do
    name=`echo $lab\_HpreCap_0+\_$sampl`
    echo $name
    while read end dist
    do
       echo "cat $name.$end.bed | sortbed | ~jlagarde/bin.bkp/bedtools2/bin/bedtools slop -s -l 5 -r 5 -i stdin -g ~jlagarde/chromInfo_hg38.txt | ~jlagarde/bin.bkp/bedtools2/bin/bedtools intersect -u -s -a stdin -b hg38.polyAsignals.bed > $name.$end.bed.vspolyAsignals.bedtsv" > $lab.$sampl.$end.polyA.sh
       qsub -cwd -N $lab.$sampl.$end.polyA -m abe -M barbara.uszczynska@crg.es -q rg-el7,long-centos79,short-centos79 -l h_rt=6:00:00 -l virtual_free=60G -pe smp 1 -e $lab.$sampl.$end.polyA.log $lab.$sampl.$end.polyA.sh
    done < ends.dist.tsv
done < ont.brain.samples.tsv

# -----------------------
# CAGE FANTOM pahse 1+2
# -----------------------

while read lab sampl
do
    name=`echo $lab\_HpreCap_0+\_$sampl`
    echo $name
    while read end dist
    do
        echo "cat $name.$end.bed | sortbed | ~jlagarde/bin.bkp/bedtools2/bin/bedtools slop -s -l 50 -r 50 -i stdin -g ~jlagarde/chromInfo_hg38.txt | ~jlagarde/bin.bkp/bedtools2/bin/bedtools intersect -u -s -a stdin -b hg38_fair+new_CAGE_peaks_phase1and2.bed > $name.$end.bed.vsCage.fantom.bedtsv" > $lab.$sampl.$end.fantomCAGE.sh
        qsub -cwd -N $lab.$sampl.$end.fantomCAGE -m abe -M barbara.uszczynska@crg.es -q rg-el7,long-centos79,short-centos79  -l h_rt=48:00:00 -l virtual_free=20G -pe smp 1 -e $lab.$sampl.$end.fantomCAGE.log $lab.$sampl.$end.fantomCAGE.sh
    done < ends.dist.tsv
done < ont.brain.samples.tsv

# ------------------------------
# get the read ids for each category
# ------------------------------

echo "
while read lab sampl
do
    name=\`echo \$lab\_HpreCap_0+\_\$sampl\`
    cat \$name.bed| cut -f4| sort| uniq > \$name.readID
    cat \$name.5.bed.vsCage.fantom.bedtsv| cut -f4| sed 's/,/\n/g'| sort|uniq > \$name.cage.readID
    cat \$name.3.bed.vspolyAsignals.bedtsv| cut -f4| sed 's/,/\n/g'| sort|uniq > \$name.polyA.readID
    cat \$name.gff| awk '{print \$10}'| sed 's/\"//g'| sed 's/;//g'| sort| uniq -c | awk '\$1==1{print \$2}' > \$name.mono.readID
    cat \$name.gff| awk '{print \$10}'| sed 's/\"//g'| sed 's/;//g'| sort| uniq -c | awk '\$1>1{print \$2}' > \$name.spliced.readID 
    cat \$name.gff| awk '{print \$10}'| sed 's/\"//g'| sed 's/;//g'| sort| uniq > \$name.all.readID 
done < ont.brain.samples.tsv " > readID.conv.sh
qsub -cwd -N readID.conv -m abe -M barbara.uszczynska@crg.es -q rg-el7,long-centos79,short-centos79  -l h_rt=48:00:00 -l virtual_free=40G -pe smp 1 -e readID.conv.log readID.conv.sh


# --------------------------------
# stats polyA reads
# --------------------------------


for typ in all spliced mono
do 
    echo "
    while read lab sampl
    do 
        name=\`echo \$lab\_HpreCap_0+\_\$sampl\`
        total=\`cat \$name.readID| wc -l\`
        All=\`cat \$name.$typ.readID| wc -l\`
        fgrep -f \$name.polyAreads.list \$name.$typ.readID > \$name.$typ.polyAreads.readID
        fgrep -f \$name.cage.readID \$name.$typ.readID > \$name.$typ.cage.polyAreads.readID
        fgrep -f \$name.$typ.polyAreads.readID \$name.$typ.cage.polyAreads.readID > \$name.$typ.FL.polyAreads.readID

        PolyA=\`fgrep -vf \$name.$typ.FL.polyAreads.readID \$name.$typ.polyAreads.readID| wc -l\`
        CAGE=\`fgrep -vf \$name.$typ.FL.polyAreads.readID \$name.$typ.cage.polyAreads.readID | wc -l\`
        FL=\`cat \$name.$typ.FL.polyAreads.readID | wc -l\`
        supported=\$((\$PolyA+\$CAGE+\$FL))
        NA=\`expr \$All - \$supported\` 
        
        #PolyAprop=\`echo \$PolyA / \$All| bc -l\`
        #CAGEprop=\`echo \$CAGE / \$All| bc -l\`
        #FLprop=\`echo \$FL / \$All| bc -l\`
        #NAprop=\`echo \$NA / \$All| bc -l\`
        
        PolyAprop=\`echo \$PolyA / \$total| bc -l\`
        CAGEprop=\`echo \$CAGE / \$total| bc -l\`
        FLprop=\`echo \$FL / \$total| bc -l\`
        NAprop=\`echo \$NA / \$total| bc -l\`
        
        echo -e \"\$lab\tHpreCap\t0+\t\$sampl\tcageAndPolyA\t\$FL\t\$FLprop\t$typ\n\$lab\tHpreCap\t0+\t\$sampl\tcageOnly\t\$CAGE\t\$CAGEprop\t$typ\n\$lab\tHpreCap\t0+\t\$sampl\tpolyAOnly\t\$PolyA\t\$PolyAprop\t$typ\n\$lab\tHpreCap\t0+\t\$sampl\tnoCageNoPolyA\t\$NA\t\$NAprop\t$typ\"
    done < ont.brain.samples.tsv  > ont.brain.end.completenes.raw.reads.stats.polyAreadslist.$typ.tsv" > end.completeness.raw.stats.polyAreadslist.$typ.sh
    qsub -cwd -N end.completeness.raw.stats.polyAreadslist.$typ -m abe -M barbara.uszczynska@crg.es -q rg-el7,long-centos79,short-centos79  -l h_rt=48:00:00 -l virtual_free=40G -pe smp 1 -e end.completeness.raw.stats.polyAreadslist.$typ.log end.completeness.raw.stats.polyAreadslist.$typ.sh
done  


# join
for typ in mono all spliced
do 
    cat ont.brain.end.completenes.raw.reads.stats.polyAreadslist.$typ.tsv
done | sed 's/all/All/g' | sed 's/spliced/Spliced/g'| sed 's/mono/Unspliced/g'| sed 's/01Rep1//g' | sed 's/01Rep3//g'| sed 's/03Rep1//g' | grep -v ont-Crg-CapTrap > ont.brain.end.completenes.raw.reads.stats.polyAreadslist.combined.tsv

echo "
library(ggplot2)
library(scales)
library(wesanderson)
cbPalette=c(\"#a6a6a6\",\"#b3e0ff\", \"#C453C4\",\"#e5b3e5\")
plot <- read.table(\"ont.brain.end.completenes.raw.reads.stats.polyAreadslist.combined.tsv\", header=F, as.is=T)
colnames(plot)<-c('seqTech', 'cap', 'sizeFrac', 'sample', 'categ', 'count', 'percent', 'cls')
plot\$categ=factor(plot\$categ, levels=c(\"noCageNoPolyA\", \"polyAOnly\", \"cageAndPolyA\", \"cageOnly\"))
plot\$seqTech=factor(plot\$seqTech, levels=c(\"ont-Cshl-CapTrap\", \"ont-Crg-CapTrap\", \"ont-Crg-dRNA\", \"ont-Crg-Telop\", \"ont-Cshl-Smarter\"))
pdf(\"ont.brain.end.completenes.raw.reads.stats.polyAreadslist.combined.pdf\", bg = \"white\", width=14, height=6)
ggplot(data=plot, aes(x=cls, y=percent, fill=categ)) + geom_bar(stat=\"identity\") + scale_y_continuous(limits=c(0,1.01), labels=percent) + ylab(\" Percent of raw reads (using polyA reads)\")+ theme_bw() + ggtitle(\"\") + scale_fill_manual(values=cbPalette)+ xlab(\"\")+theme(axis.text.x = element_text(size = 14, colour = \"black\", vjust=0.5), axis.text.y = element_text(size = 14,  colour = \"black\"), legend.text = element_text(size=12), plot.title = element_text(size = 22), legend.title =element_text(size=0, color=\"white\"), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14), legend.text.align=0.5, strip.text.x=element_text(size=18), strip.text.y=element_text(size=18))+ ggtitle('') + theme(legend.position=\"top\") +facet_grid(sample ~ seqTech, scales=\"free_x\")+ geom_text(data=plot, aes(x=cls, y=percent, label =paste(round(percent*100,1),\"(\",comma(count),\")\", sep=\"\")), size=3, position = position_stack(vjust = .5))
dev.off()
"| R --slave


# ---------------------------------------------------------------------------
# stranded HCGMs - completeness
# ---------------------------------------------------------------------------

while read lab sampl
do
    name=`echo $lab\_HpreCap_0+\_$sampl.strandedHCGMs`
    lid=`echo $lab\_HpreCap_0+\_$sampl`
    ln -s /nfs/users/project/gencode_006070_no_backup/epalumbo/pilot/output/mappings/highConfidenceReads/$name.gff.gz .
    zcat $name.gff.gz > $name.gff
    echo "gff2bed_full.pl $name.gff > $name.bed" > $name.bed2gff.sh
    qsub -cwd -N $name.bed2gff -m abe -M barbara.uszczynska@crg.es -q rg-el7,long-centos79,short-centos79  -l h_rt=6:00:00 -l virtual_free=50G -pe smp 1 -e $name.bed2gff.log $name.bed2gff.sh  
done < ont.brain.samples.tsv

# ------------------------------
# hcgms - ends
# ------------------------------

while read lab sampl
do
    name=`echo $lab\_HpreCap_0+\_$sampl.strandedHCGMs`
    while read end dist
    do
        echo "cat $name.bed | /users/rg/jlagarde/julien_utils_public/extractTranscriptEndsFromBed12.pl $end |sortbed  > $name.$end.bed" > $lab.$sampl.$end.hcgms.endsExtr.sh
        qsub -cwd -N $lab.$sampl.$end.hcgms.endsExtr -m abe -M barbara.uszczynska@crg.es -q rg-el7,long-centos79,short-centos79  -l h_rt=48:00:00 -l virtual_free=40G -pe smp 1 -e $lab.$sampl.$end.hcgms.endsExtr.log $lab.$sampl.$end.hcgms.endsExtr.sh
    done < ends.dist.tsv
done < ont.brain.samples.tsv

# -------------------
# polyA
# -------------------

while read lab sampl
do
    name=`echo $lab\_HpreCap_0+\_$sampl.strandedHCGMs`
    echo $name
    while read end dist
    do
       echo "cat $name.$end.bed | sortbed | ~jlagarde/bin.bkp/bedtools2/bin/bedtools slop -s -l 5 -r 5 -i stdin -g ~jlagarde/chromInfo_hg38.txt | ~jlagarde/bin.bkp/bedtools2/bin/bedtools intersect -u -s -a stdin -b hg38.polyAsignals.bed > $name.$end.bed.vspolyAsignals.bedtsv" > $lab.$sampl.$end.hcgms.polyA.sh
       qsub -cwd -N $lab.$sampl.$end.hcgms.polyA -m abe -M barbara.uszczynska@crg.es -q rg-el7,long-centos79,short-centos79 -l h_rt=6:00:00 -l virtual_free=60G -pe smp 1 -e $lab.$sampl.$end.hcgms.polyA.log $lab.$sampl.$end.hcgms.polyA.sh
    done < ends.dist.tsv
done < ont.brain.samples.tsv

# -----------------------
# CAGE FANTOM pahse 1+2
# -----------------------

while read lab sampl
do
    name=`echo $lab\_HpreCap_0+\_$sampl.strandedHCGMs`
    echo $name
    while read end dist
    do
        echo "cat $name.$end.bed | sortbed | ~jlagarde/bin.bkp/bedtools2/bin/bedtools slop -s -l 50 -r 50 -i stdin -g ~jlagarde/chromInfo_hg38.txt | ~jlagarde/bin.bkp/bedtools2/bin/bedtools intersect -u -s -a stdin -b hg38_fair+new_CAGE_peaks_phase1and2.bed > $name.$end.bed.vsCage.fantom.bedtsv" > $lab.$sampl.$end.hcgms.fantomCAGE.sh
        qsub -cwd -N $lab.$sampl.$end.hcgms.fantomCAGE -m abe -M barbara.uszczynska@crg.es -q rg-el7,long-centos79,short-centos79  -l h_rt=48:00:00 -l virtual_free=20G -pe smp 1 -e $lab.$sampl.$end.hcgms.fantomCAGE.log $lab.$sampl.$end.hcgms.fantomCAGE.sh
    done < ends.dist.tsv
done < ont.brain.samples.tsv

# ------------------------------
# get the read ids for each category
# ------------------------------

echo "
while read lab sampl
do
    name=\`echo \$lab\_HpreCap_0+\_\$sampl.strandedHCGMs\`
    cat \$name.bed| cut -f4| sort| uniq > \$name.readID
    cat \$name.5.bed.vsCage.fantom.bedtsv| cut -f4| sed 's/,/\n/g'| sort|uniq > \$name.cage.readID
    cat \$name.3.bed.vspolyAsignals.bedtsv| cut -f4| sed 's/,/\n/g'| sort|uniq > \$name.polyA.readID
    cat \$name.gff| awk '{print \$10}'| sed 's/\"//g'| sed 's/;//g'| sort| uniq -c | awk '\$1==1{print \$2}' > \$name.mono.readID
    cat \$name.gff| awk '{print \$10}'| sed 's/\"//g'| sed 's/;//g'| sort| uniq -c | awk '\$1>1{print \$2}' > \$name.spliced.readID 
    cat \$name.gff| awk '{print \$10}'| sed 's/\"//g'| sed 's/;//g'| sort| uniq > \$name.all.readID 
done < ont.brain.samples.tsv " > readID.conv.hcgms.sh
qsub -cwd -N readID.conv.hcgms -m abe -M barbara.uszczynska@crg.es -q rg-el7,long-centos79,short-centos79  -l h_rt=48:00:00 -l virtual_free=40G -pe smp 1 -e readID.conv.hcgms.log readID.conv.hcgms.sh

# --------------------------------
# stats polyA reads
# --------------------------------


for typ in all spliced mono
do 
    echo "
    while read lab sampl
    do 
        name=\`echo \$lab\_HpreCap_0+\_\$sampl.strandedHCGMs\`
        total=\`cat \$name.readID| wc -l\`
        lid=\`echo \$lab\_HpreCap_0+\_\$sampl\`
        All=\`cat \$name.$typ.readID| wc -l\`
        fgrep -f \$lid.polyAreads.list \$name.$typ.readID > \$name.$typ.polyAreads.readID
        fgrep -f \$lid.cage.readID \$name.$typ.readID > \$name.$typ.cage.polyAreads.readID
        fgrep -f \$name.$typ.polyAreads.readID \$name.$typ.cage.polyAreads.readID > \$name.$typ.FL.polyAreads.readID

        PolyA=\`fgrep -vf \$name.$typ.FL.polyAreads.readID \$name.$typ.polyAreads.readID| wc -l\`
        CAGE=\`fgrep -vf \$name.$typ.FL.polyAreads.readID \$name.$typ.cage.polyAreads.readID | wc -l\`
        FL=\`cat \$name.$typ.FL.polyAreads.readID | wc -l\`
        supported=\$((\$PolyA+\$CAGE+\$FL))
        NA=\`expr \$All - \$supported\` 
        
        #PolyAprop=\`echo \$PolyA / \$All| bc -l\`
        #CAGEprop=\`echo \$CAGE / \$All| bc -l\`
        #FLprop=\`echo \$FL / \$All| bc -l\`
        #NAprop=\`echo \$NA / \$All| bc -l\`
        
        PolyAprop=\`echo \$PolyA / \$total| bc -l\`
        CAGEprop=\`echo \$CAGE / \$total| bc -l\`
        FLprop=\`echo \$FL / \$total| bc -l\`
        NAprop=\`echo \$NA / \$total| bc -l\`
        
        
        echo -e \"\$lab\tHpreCap\t0+\t\$sampl\tcageAndPolyA\t\$FL\t\$FLprop\t$typ\n\$lab\tHpreCap\t0+\t\$sampl\tcageOnly\t\$CAGE\t\$CAGEprop\t$typ\n\$lab\tHpreCap\t0+\t\$sampl\tpolyAOnly\t\$PolyA\t\$PolyAprop\t$typ\n\$lab\tHpreCap\t0+\t\$sampl\tnoCageNoPolyA\t\$NA\t\$NAprop\t$typ\"
    done < ont.brain.samples.tsv  > ont.brain.end.completenes.raw.reads.stats.polyAreadslist.hcgms.$typ.tsv" > end.completeness.raw.stats.polyAreadslist.hcgms.$typ.sh
    qsub -cwd -N end.completeness.raw.stats.polyAreadslist.hcgms.$typ -m abe -M barbara.uszczynska@crg.es -q rg-el7,long-centos79,short-centos79  -l h_rt=48:00:00 -l virtual_free=40G -pe smp 1 -e end.completeness.raw.stats.polyAreadslist.hcgms.$typ.log end.completeness.raw.stats.polyAreadslist.hcgms.$typ.sh
done  


# join
for typ in mono all spliced
do 
    cat ont.brain.end.completenes.raw.reads.stats.polyAreadslist.hcgms.$typ.tsv| awk 'BEGIN{FS=OFS="\t"}{split($1,a,"_"); print a[1], $2, $3, $4, $5, $6, $7, $8}'
done | sed 's/all/All/g' | sed 's/spliced/Spliced/g'| sed 's/mono/Unspliced/g'| sed 's/01Rep1//g' | sed 's/01Rep3//g'| sed 's/03Rep1//g' | grep -v ont-Crg-CapTrap > ont.brain.end.completenes.raw.reads.stats.polyAreadslist.hcgms.combined.tsv

echo "
library(ggplot2)
library(scales)
library(wesanderson)
cbPalette=c(\"#a6a6a6\",\"#b3e0ff\", \"#C453C4\",\"#e5b3e5\")
plot <- read.table(\"ont.brain.end.completenes.raw.reads.stats.polyAreadslist.hcgms.combined.tsv\", header=F, as.is=T)
colnames(plot)<-c('seqTech', 'cap', 'sizeFrac', 'sample', 'categ', 'count', 'percent', 'cls')
plot\$categ=factor(plot\$categ, levels=c(\"noCageNoPolyA\", \"polyAOnly\", \"cageAndPolyA\", \"cageOnly\"))
plot\$seqTech=factor(plot\$seqTech, levels=c(\"ont-Cshl-CapTrap\", \"ont-Crg-CapTrap\", \"ont-Crg-dRNA\", \"ont-Crg-Telop\", \"ont-Cshl-Smarter\"))
pdf(\"ont.brain.end.completenes.raw.reads.stats.polyAreadslist.hcgms.combined.pdf\", bg = \"white\", width=14, height=6)
ggplot(data=plot, aes(x=cls, y=percent, fill=categ)) + geom_bar(stat=\"identity\") + scale_y_continuous(limits=c(0,1.01), labels=percent) + ylab(\" Percent of HCGMs (using polyA reads) \")+ theme_bw() + ggtitle(\"\") + scale_fill_manual(values=cbPalette)+ xlab(\"\")+theme(axis.text.x = element_text(size = 14, colour = \"black\", vjust=0.5), axis.text.y = element_text(size = 14,  colour = \"black\"), legend.text = element_text(size=12), plot.title = element_text(size = 22), legend.title =element_text(size=0, color=\"white\"), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14), legend.text.align=0.5, strip.text.x=element_text(size=18), strip.text.y=element_text(size=18))+ ggtitle('') + theme(legend.position=\"top\") +facet_grid(sample ~ seqTech, scales=\"free_x\")+ geom_text(data=plot, aes(x=cls, y=percent, label =paste(round(percent*100,1),\"(\",comma(count),\")\", sep=\"\")), size=3, position = position_stack(vjust = .5))
dev.off()
"| R --slave

# ---------------------------------------------------------------------------
# stranded HiSS - completeness
# ---------------------------------------------------------------------------
while read lab sampl
do
    name=`echo $lab\_HpreCap_0+\_$sampl.HiSS`
    lid=`echo $lab\_HpreCap_0+\_$sampl`
    #ln -s /nfs/users/project/gencode_006070_no_backup/epalumbo/pilot/output/mappings/highConfidenceReads/HiSS/$name.gff.gz .
    #zcat $name.gff.gz > $name.gff
    echo "gff2bed_full.pl $name.gff > $name.bed" > $name.bed2gff.sh
    qsub -cwd -N $name.bed2gff -m abe -M barbara.uszczynska@crg.es -q rg-el7,long-centos79,short-centos79  -l h_rt=6:00:00 -l virtual_free=50G -pe smp 1 -e $name.bed2gff.log $name.bed2gff.sh  
done < ont.brain.samples.tsv

# ------------------------------
# hiss - ends
# ------------------------------

while read lab sampl
do
    name=`echo $lab\_HpreCap_0+\_$sampl.HiSS`
    while read end dist
    do
        echo "cat $name.bed | /users/rg/jlagarde/julien_utils_public/extractTranscriptEndsFromBed12.pl $end |sortbed  > $name.$end.bed" > $lab.$sampl.$end.hiss.endsExtr.sh
        qsub -cwd -N $lab.$sampl.$end.hiss.endsExtr -m abe -M barbara.uszczynska@crg.es -q rg-el7,long-centos79,short-centos79  -l h_rt=48:00:00 -l virtual_free=40G -pe smp 1 -e $lab.$sampl.$end.hiss.endsExtr.log $lab.$sampl.$end.hiss.endsExtr.sh
    done < ends.dist.tsv
done < ont.brain.samples.tsv

# -------------------
# polyA
# -------------------

while read lab sampl
do
    name=`echo $lab\_HpreCap_0+\_$sampl.HiSS`
    echo $name
    while read end dist
    do
       echo "cat $name.$end.bed | sortbed | ~jlagarde/bin.bkp/bedtools2/bin/bedtools slop -s -l 5 -r 5 -i stdin -g ~jlagarde/chromInfo_hg38.txt | ~jlagarde/bin.bkp/bedtools2/bin/bedtools intersect -u -s -a stdin -b hg38.polyAsignals.bed > $name.$end.bed.vspolyAsignals.bedtsv" > $lab.$sampl.$end.hiss.polyA.sh
       qsub -cwd -N $lab.$sampl.$end.hiss.polyA -m abe -M barbara.uszczynska@crg.es -q rg-el7,long-centos79,short-centos79 -l h_rt=6:00:00 -l virtual_free=60G -pe smp 1 -e $lab.$sampl.$end.hiss.polyA.log $lab.$sampl.$end.hiss.polyA.sh
    done < ends.dist.tsv
done < ont.brain.samples.tsv

# -----------------------
# CAGE FANTOM pahse 1+2
# -----------------------

while read lab sampl
do
    name=`echo $lab\_HpreCap_0+\_$sampl.HiSS`
    echo $name
    while read end dist
    do
        echo "cat $name.$end.bed | sortbed | ~jlagarde/bin.bkp/bedtools2/bin/bedtools slop -s -l 50 -r 50 -i stdin -g ~jlagarde/chromInfo_hg38.txt | ~jlagarde/bin.bkp/bedtools2/bin/bedtools intersect -u -s -a stdin -b hg38_fair+new_CAGE_peaks_phase1and2.bed > $name.$end.bed.vsCage.fantom.bedtsv" > $lab.$sampl.$end.hiss.fantomCAGE.sh
        qsub -cwd -N $lab.$sampl.$end.hiss.fantomCAGE -m abe -M barbara.uszczynska@crg.es -q rg-el7,long-centos79,short-centos79  -l h_rt=48:00:00 -l virtual_free=20G -pe smp 1 -e $lab.$sampl.$end.hiss.fantomCAGE.log $lab.$sampl.$end.hiss.fantomCAGE.sh
    done < ends.dist.tsv
done < ont.brain.samples.tsv

# ------------------------------
# get the read ids for each category
# ------------------------------

echo "
while read lab sampl
do
    name=\`echo \$lab\_HpreCap_0+\_\$sampl.HiSS\`
    cat \$name.bed| cut -f4| sort| uniq > \$name.readID
    cat \$name.5.bed.vsCage.fantom.bedtsv| cut -f4| sed 's/,/\n/g'| sort|uniq > \$name.cage.readID
    cat \$name.3.bed.vspolyAsignals.bedtsv| cut -f4| sed 's/,/\n/g'| sort|uniq > \$name.polyA.readID
    cat \$name.gff| awk '{print \$10}'| sed 's/\"//g'| sed 's/;//g'| sort| uniq -c | awk '\$1==1{print \$2}' > \$name.mono.readID
    cat \$name.gff| awk '{print \$10}'| sed 's/\"//g'| sed 's/;//g'| sort| uniq -c | awk '\$1>1{print \$2}' > \$name.spliced.readID 
    cat \$name.gff| awk '{print \$10}'| sed 's/\"//g'| sed 's/;//g'| sort| uniq > \$name.all.readID 
done < ont.brain.samples.tsv " > readID.conv.hiss.sh
qsub -cwd -N readID.conv.hiss -m abe -M barbara.uszczynska@crg.es -q rg-el7,long-centos79,short-centos79  -l h_rt=48:00:00 -l virtual_free=40G -pe smp 1 -e readID.conv.hiss.log readID.conv.hiss.sh


# --------------------------------
# stats polyA reads
# --------------------------------


for typ in all spliced mono
do 
    echo "
    while read lab sampl
    do 
        name=\`echo \$lab\_HpreCap_0+\_\$sampl.HiSS\`
        total=\`cat \$name.readID| wc -l\`
        lid=\`echo \$lab\_HpreCap_0+\_\$sampl\`
        All=\`cat \$name.$typ.readID| wc -l\`
        fgrep -f \$lid.polyAreads.list \$name.$typ.readID > \$name.$typ.polyAreads.readID
        fgrep -f \$lid.cage.readID \$name.$typ.readID > \$name.$typ.cage.polyAreads.readID
        fgrep -f \$name.$typ.polyAreads.readID \$name.$typ.cage.polyAreads.readID > \$name.$typ.FL.polyAreads.readID

        PolyA=\`fgrep -vf \$name.$typ.FL.polyAreads.readID \$name.$typ.polyAreads.readID| wc -l\`
        CAGE=\`fgrep -vf \$name.$typ.FL.polyAreads.readID \$name.$typ.cage.polyAreads.readID | wc -l\`
        FL=\`cat \$name.$typ.FL.polyAreads.readID | wc -l\`
        supported=\$((\$PolyA+\$CAGE+\$FL))
        NA=\`expr \$All - \$supported\` 
        
        #PolyAprop=\`echo \$PolyA / \$All| bc -l\`
        #CAGEprop=\`echo \$CAGE / \$All| bc -l\`
        #FLprop=\`echo \$FL / \$All| bc -l\`
        #NAprop=\`echo \$NA / \$All| bc -l\`
        
        PolyAprop=\`echo \$PolyA / \$total| bc -l\`
        CAGEprop=\`echo \$CAGE / \$total| bc -l\`
        FLprop=\`echo \$FL / \$total| bc -l\`
        NAprop=\`echo \$NA / \$total| bc -l\`
        
        
        echo -e \"\$lab\tHpreCap\t0+\t\$sampl\tcageAndPolyA\t\$FL\t\$FLprop\t$typ\n\$lab\tHpreCap\t0+\t\$sampl\tcageOnly\t\$CAGE\t\$CAGEprop\t$typ\n\$lab\tHpreCap\t0+\t\$sampl\tpolyAOnly\t\$PolyA\t\$PolyAprop\t$typ\n\$lab\tHpreCap\t0+\t\$sampl\tnoCageNoPolyA\t\$NA\t\$NAprop\t$typ\"
    done < ont.brain.samples.tsv  > ont.brain.end.completenes.raw.reads.stats.polyAreadslist.hiss.$typ.tsv" > end.completeness.raw.stats.polyAreadslist.hiss.$typ.sh
    qsub -cwd -N end.completeness.raw.stats.polyAreadslist.hiss.$typ -m abe -M barbara.uszczynska@crg.es -q rg-el7,long-centos79,short-centos79  -l h_rt=48:00:00 -l virtual_free=40G -pe smp 1 -e end.completeness.raw.stats.polyAreadslist.hiss.$typ.log end.completeness.raw.stats.polyAreadslist.hiss.$typ.sh
done  


# join
for typ in mono all spliced
do 
    cat ont.brain.end.completenes.raw.reads.stats.polyAreadslist.hiss.$typ.tsv| awk 'BEGIN{FS=OFS="\t"}{split($1,a,"_"); print a[1], $2, $3, $4, $5, $6, $7, $8}'
done | sed 's/all/All/g' | sed 's/spliced/Spliced/g'| sed 's/mono/Unspliced/g'| sed 's/01Rep1//g' | sed 's/01Rep3//g'| sed 's/03Rep1//g' | grep -v ont-Crg-CapTrap > ont.brain.end.completenes.raw.reads.stats.polyAreadslist.hiss.combined.tsv

echo "
library(ggplot2)
library(scales)
library(wesanderson)
cbPalette=c(\"#a6a6a6\",\"#b3e0ff\", \"#C453C4\",\"#e5b3e5\")
plot <- read.table(\"ont.brain.end.completenes.raw.reads.stats.polyAreadslist.hiss.combined.tsv\", header=F, as.is=T)
colnames(plot)<-c('seqTech', 'cap', 'sizeFrac', 'sample', 'categ', 'count', 'percent', 'cls')
plot\$categ=factor(plot\$categ, levels=c(\"noCageNoPolyA\", \"polyAOnly\", \"cageAndPolyA\", \"cageOnly\"))
plot\$seqTech=factor(plot\$seqTech, levels=c(\"ont-Cshl-CapTrap\", \"ont-Crg-CapTrap\", \"ont-Crg-dRNA\", \"ont-Crg-Telop\", \"ont-Cshl-Smarter\"))
pdf(\"ont.brain.end.completenes.raw.reads.stats.polyAreadslist.hiss.combined.pdf\", bg = \"white\", width=14, height=6)
ggplot(data=plot, aes(x=cls, y=percent, fill=categ)) + geom_bar(stat=\"identity\") + scale_y_continuous(limits=c(0,1.01), labels=percent) + ylab(\" Percent of HiSS (using polyA reads) \")+ theme_bw() + ggtitle(\"\") + scale_fill_manual(values=cbPalette)+ xlab(\"\")+theme(axis.text.x = element_text(size = 14, colour = \"black\", vjust=0.5), axis.text.y = element_text(size = 14,  colour = \"black\"), legend.text = element_text(size=12), plot.title = element_text(size = 22), legend.title =element_text(size=0, color=\"white\"), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14), legend.text.align=0.5, strip.text.x=element_text(size=18), strip.text.y=element_text(size=18))+ ggtitle('') + theme(legend.position=\"top\") +facet_grid(sample ~ seqTech, scales=\"free_x\")+ geom_text(data=plot, aes(x=cls, y=percent, label =paste(round(percent*100,1),\"(\",comma(count),\")\", sep=\"\")), size=3, position = position_stack(vjust = .5))
dev.off()
"| R --slave

