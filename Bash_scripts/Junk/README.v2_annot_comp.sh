working=/users/rg/buszczynska/Projects/review/human/annot.comp

cd $working

# The goal is to compare different lncRNA annotations: RefSeq, NonCode, MiTranscriptome, FANTOM CAT, Big Transcriptome

mkdir annots

# ------------------------------------------------------------------------------------------------------------------
# Get the gene sets
# ------------------------------------------------------------------------------------------------------------------

# -------------------------------------------------
# GENCODE (v27) (hg38)
# -------------------------------------------------
ln /users/rg/buszczynska/Projects/review/human/gen.v27.lncRNAs.gtf gen.lncRNAs.gtf

# -------------------------------------------------
# NONCODE (v5) (hg38)
# -------------------------------------------------
wget http://www.noncode.org/datadownload/NONCODEv5_human_hg38_lncRNA.gtf.gz .
zcat NONCODEv5_human_hg38_lncRNA.gtf.gz > noncode.lncRNAs.gtf

mv NONCODEv5_human_hg38_lncRNA.gtf.gz annots/

# -------------------------------------------------
# RefSeq (GCF_000001405) (hg38)
# -------------------------------------------------

# copy
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/latest_assembly_versions/GCF_000001405.37_GRCh38.p11/GCF_000001405.37_GRCh38.p11_genomic.gff.gz .

# ungzip
gunzip GCF_000001405.37_GRCh38.p11_genomic.gff.gz

# get the gtf
/users/rg/buszczynska/python_scripts/gffTogtf.py -i GCF_000001405.37_GRCh38.p11_genomic.gff -o refseq.gtf

# get the gtf for lncRNAs
cat refseq.gtf | awk '$3=="gene"'| awk '$14=="\"lncRNA\";"{print $10}'| sed 's/;//g'| sed 's/"//g'| sort| uniq| fgrep -f - refseq.gtf | sed 's/lnc_RNA/transcript/g'> refseq.lncRNAs.gtf

# -------------------------------------------------
# miTranscriptome (hg19)
# -------------------------------------------------

wget http://mitranscriptome.org/download/mitranscriptome.gtf.tar.gz .
# unpack and use mitranscriptome.v2.gtf.gz

# gunzip
gunzip mitranscriptome.v2.gtf.gz

# gtf to gencode like gtf
/users/rg/buszczynska/python_scripts/mitranscriptome.py -i mitranscriptome.v2.gtf -o mitranscriptome.v2.std.gtf

# get the gtf for lncRNAs
cat mitranscriptome.v2.std.gtf | awk '$3=="transcript"'| awk '$14=="\"lncrna\";"{print $10}'| sed 's/;//g'| sed 's/"//g'| sort| uniq| fgrep -f - mitranscriptome.v2.std.gtf > mitrans.lncRNAs.gtf
 
# -------------------------------------------------
# Fantom CAT (hg19)
# -------------------------------------------------

wget http://fantom.gsc.riken.jp/5/suppl/Hon_et_al_2016/data/assembly/lv3_robust/FANTOM_CAT.lv3_robust.gtf.gz .

# gunzip
gunzip FANTOM_CAT.lv3_robust.gtf.gz

# gtf for lncRNAs 
cat FANTOM_CAT.lv3_robust.gtf| awk '$3=="gene"'| awk '{if($14=="\"lncRNA_antisense\";" || $14=="\"lncRNA_divergent\";" || $14=="\"lncRNA_intergenic\";" || $14=="\"lncRNA_sense_intronic\";") print $10}'| sed 's/;//g'| sed 's/"//g'| sort| uniq| fgrep -f - FANTOM_CAT.lv3_robust.gtf > fantomCat.lncRNAs.gtf

# -------------------------------------------------
# BigTranscriptome (hg19)
# -------------------------------------------------

#wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE97nnn/GSE97211/suppl/GSE97211_BIGTranscriptome.gtf.gz .
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE97nnn/GSE97211/suppl/GSE97211_BIGTranscriptome_lncRNA_catalog.gtf.gz .

# gtf for lncRNAs
zcat GSE97211_BIGTranscriptome_lncRNA_catalog.gtf.gz > bigtrans.lncRNAs.gtf

# -------------------------------------------------
# CLS (hg38)
# -------------------------------------------------

# copy
wget https://public_docs.crg.es/rguigo/CLS/data/gencodePlus/hs.GENCODE+.lncRNA.bed .

# change the name
mv hs.GENCODE+.lncRNA.bed cls.hg38.bed12

# get gtf
awk -f /no_backup_isis/rg/0ld_users/sdjebali/Awk/bed12fields2gff.awk cls.hg38.bed12 > cls.lncRNAs.hg38.gtf

# get txs
/users/rg/buszczynska/python_scripts/extract_trans_coord.py -i cls.lncRNAs.hg38.gtf -o cls.lncRNAs.hg38.txs.gtf

# build genes
~jlagarde/bin/bedtools2/bin/bedtools intersect -s -wao -a cls.lncRNAs.hg38.txs.gtf -b cls.lncRNAs.hg38.txs.gtf | ~jlagarde/julien_utils/buildLoci.pl - > cls.lncRNAs.hg38.gene.tmp.gtf
/users/rg/jlagarde/julien_utils/gff2gff.pl cls.lncRNAs.hg38.gene.tmp.gtf > cls.lncRNAs.hg38.gene.gtf

# get the reference 
cat cls.lncRNAs.hg38.gene.gtf| awk 'BEGIN{OFS="\t"}{print $12,$10}'| sort| uniq| sed 's/"//g'| sed 's/;//g' > cls.lncRNAs.hg38.gene.txs.ref

# exchange the geneIDs
/users/rg/buszczynska/python_scripts/exchange.gene_id.py -i cls.lncRNAs.hg38.gtf -o cls.lncRNAs.gtf -r cls.lncRNAs.hg38.gene.txs.ref

# -------------------------------------------------
# Get the full legth CLS (hg38)
# -------------------------------------------------

# get the anchored
wget https://public_docs.crg.es/rguigo/CLS/data/mergedTMs/non-anchored/gtf/hsAll_Cap1_all_noAnchor.compmerge.cage+polyASupported.gtf.gz . 
gunzip hsAll_Cap1_all_noAnchor.compmerge.cage+polyASupported.gtf.gz 

# get the full length ids
~abreschi/Documents/utils/extract.gtf.tags.sh hsAll_Cap1_all_noAnchor.compmerge.cage+polyASupported.gtf transcript_id | sort| uniq > clsFLAll.txIds

# join to get bed12
join.py -a clsFLAll.txIds -b cls.hg38.bed12 -x 1 -y 4 > clsFL.hg38.bed12

# get gtf from bed12
awk -f /no_backup_isis/rg/0ld_users/sdjebali/Awk/bed12fields2gff.awk clsFL.hg38.bed12 > clsFL.lncRNAs.hg38.gtf

# get txs
/users/rg/buszczynska/python_scripts/extract_trans_coord.py -i clsFL.lncRNAs.hg38.gtf -o clsFL.lncRNAs.hg38.txs.gtf

# build genes
~jlagarde/bin/bedtools2/bin/bedtools intersect -s -wao -a clsFL.lncRNAs.hg38.txs.gtf -b clsFL.lncRNAs.hg38.txs.gtf | ~jlagarde/julien_utils/buildLoci.pl - > clsFL.lncRNAs.hg38.gene.tmp.gtf
/users/rg/jlagarde/julien_utils/gff2gff.pl clsFL.lncRNAs.hg38.gene.tmp.gtf > clsFL.lncRNAs.hg38.gene.gtf

# get the reference 
cat clsFL.lncRNAs.hg38.gene.gtf| awk 'BEGIN{OFS="\t"}{print $12,$10}'| sort| uniq| sed 's/"//g'| sed 's/;//g' > clsFL.lncRNAs.hg38.gene.txs.ref

# exchange the geneIDs
/users/rg/buszczynska/python_scripts/exchange.gene_id.py -i clsFL.lncRNAs.hg38.gtf -o clsFL.lncRNAs.gtf -r clsFL.lncRNAs.hg38.gene.txs.ref

#cat /users/rg/jlagarde/projects/encode/scaling/whole_genome/lncRnaCapture/pacBioVsHiSeq/stringtieVsPacBio/transcripts/5pEnds_support/hs.CLS-noAnchor.onTargetLncRna.5pEndsClusters.cage+polyASupported.bed|cut -f4|sed 's/,/\n/g'| sort|uniq > clsFL.txIds

# -------------------------------------------------
# PC genes
# -------------------------------------------------

wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_27/gencode.v27.annotation.gtf.gz .

gunzip gencode.v27.annotation.gtf.gz

# get protein coding trascripts with mRNA start or stop not found
~abreschi/Documents/utils/extract.gtf.tags.sh gencode.v27.annotation.gtf transcript_id,tag| grep -E "mRNA_start_NF|mRNA_end_NF"| cut -f1| sort| uniq > gen27.mrna.start.end.nf.txID

# QC
#cat gencode.v27.annotation.gtf| grep -E "mRNA_start_NF|mRNA_end_NF"| awk '$3=="transcript"{print $12}'| sort| uniq| wc -l

# Check txs biotypes
fgrep -f gen27.mrna.start.end.nf.txID gencode.v27.annotation.gtf| awk '$3=="transcript"{print $18}'|sort| uniq -c

# get the txs
#cat gencode.v27.annotation.gtf| awk '$3=="transcript"'| awk '$18=="\"protein_coding\";"{print $12}'| sed 's/;//g'| sed 's/"//g'| sort| uniq| fgrep -f - gencode.v27.annotation.gtf  > pc.lncRNAs.gtf
cat gencode.v27.annotation.gtf| awk '$3=="transcript"'| awk '$18=="\"protein_coding\";"{print $12}'| sed 's/;//g'| sed 's/"//g'| sort| uniq| fgrep -f - gencode.v27.annotation.gtf| fgrep -vf gen27.mrna.start.end.nf.txID -  > pcConf.lncRNAs.gtf

# ------------------------------------------------------------------------------------------------------------------
#                                               *** ANALYSIS ***
# ------------------------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------
# smaples
# --------------------------------------------------------------------

ls *.lncRNAs.gtf| awk -F "." '{print $1}'| sort| uniq > samples.tsv
cat samples.tsv | grep -v cls > samples.nocls.tsv

# ------------------------------------------------------------------------------------------------------------------
#                                               *** LiftOver ***
# ------------------------------------------------------------------------------------------------------------------

# --------------------------------------------
# BED FILES
# --------------------------------------------

while read lab
do
   gff2bed_full.pl $lab.lncRNAs.gtf | awk 'BEGIN{FS=OFS="\t"}$6!="."' > $lab.bed12
done < samples.tsv 

# --------------------------------------------
# BED FILES (LiftOver to hg38)
# --------------------------------------------

for lab in bigtrans fantomCat mitrans
do
    liftOver $lab.bed12 /users/rg/buszczynska/liftOver/hg19ToHg38.over.chain $lab.hg38.bed12 $lab.unMapped
done

for lab in gen refseq noncode pcConf
do
    cat $lab.bed12 | awk 'BEGIN{FS=OFS="\t"}$6!="."' > $lab.hg38.bed12
done

# --------------------------------------------------------------------
# GTF from hg38.bed12
# --------------------------------------------------------------------

while read lab
do
    awk -f /no_backup_isis/rg/0ld_users/sdjebali/Awk/bed12fields2gff.awk $lab.hg38.bed12 > $lab.lncRNAs.hg38.gtf
done < samples.nocls.tsv 

# --------------------------------------------------------------------
# Get the txs gtf
# --------------------------------------------------------------------

while read lab
do
    /users/rg/buszczynska/python_scripts/extract_trans_coord.py -i $lab.lncRNAs.hg38.gtf -o $lab.lncRNAs.hg38.txs.gtf
done < samples.nocls.tsv


# --------------------------------------------------------------------
# build the gene ids
# --------------------------------------------------------------------

while read lab
do
    echo $lab
    ~jlagarde/bin/bedtools2/bin/bedtools intersect -s -wao -a $lab.lncRNAs.hg38.txs.gtf -b $lab.lncRNAs.hg38.txs.gtf | ~jlagarde/julien_utils/buildLoci.pl - > $lab.lncRNAs.hg38.gene.tmp.gtf
    /users/rg/jlagarde/julien_utils/gff2gff.pl $lab.lncRNAs.hg38.gene.tmp.gtf > $lab.lncRNAs.hg38.gene.gtf
done < samples.nocls.tsv

# --------------------------------------------------------------------
# get the gene ids from the hg38.gtfs
# --------------------------------------------------------------------

while read lab
do
    cat $lab.lncRNAs.hg38.gene.gtf| awk '{print $10}'| sed 's/;//g'| sed 's/"//g' | sort| uniq > $lab.buildLoci.geneIds
done < samples.tsv 

# summary stats (to see if the numbers correspond to what is on gencode web)
while read lab
do
    nb=`cat $lab.buildLoci.geneIds| wc -l`
    echo -e "$lab\t$nb\tbuildLoci"
done < samples.tsv | sed 's/gen/GENCODE/g'| sed 's/refseq/RefSeq/g'|sed 's/fantomCat/FANTOM CAT/g'| sed 's/mitrans/MiTranscriptome/g'| sed 's/bigtrans/BIGTranscriptome/g'| sed 's/noncode/NONCODE/g' | sed 's/cls/GENCODE+/g'| sed 's/pcConf/Protein Coding/g'| sed 's/GENCODE+FL/CLS FL/g' | tee annots.ref.rel.gene.nb.summary.stats.buildLoci.hg38.tsv

# plot
#echo "
#library(ggplot2)
#library(scales)
#plot <- read.table(\"annots.ref.rel.gene.nb.summary.stats.buildLoci.hg38.tsv\", header=F, as.is=T, sep=\"\t\")
#colnames(plot)<-c(\"gene\",\"count\",\"typ\")
#plot\$gene=factor(plot\$gene, levels=c(\"NONCODE\",\"MiTranscriptome\",\"FANTOM CAT\",\"RefSeq\",\"GENCODE\",\"BIGTranscriptome\",\"GENCODE+\"))
#png(\"annots.ref.rel.gene.nb.summary.stats.buildLoci.hg38.human.png\", bg = \"white\", units=\"in\", width=10, height=5, res=300)
#ggplot(data=plot, aes(x=gene, y=count, fill=typ)) + geom_bar(stat=\"identity\") + ylab(\" \")+ theme_bw() + ggtitle(\"\") + scale_y_continuous(limits=c(0, 100000), breaks=seq(0,100000, by=20000))+ scale_fill_manual(values=c(\"#0d6598\"))+ xlab(\"\")+ ylab (\"Number of genes\")+theme(axis.text.x = element_text(size = 12.5, colour = \"black\"), axis.text.y = element_text(size = 14,  colour = \"black\"), legend.text = element_text(size=12), plot.title = element_text(size = 22), legend.title =element_text(size=0, color=\"white\"), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14), legend.text.align=0.5, strip.text.x=element_text(size=18), strip.text.y=element_text(size=18))+ ggtitle('') + theme(legend.position=\"none\") + geom_text(position = \"stack\", aes(x = gene, y = count, label = comma(count), hjust = 0.5, vjust = 0.05), size=4)+facet_wrap(~ typ)
#dev.off()
#" | R --slave

#cp annots.ref.rel.gene.nb.summary.stats.buildLoci.hg38.human.png /users/rg/buszczynska/public_html/gencode/tmp/B/review

# ----------------------------------------------------------------------------------
# get the gene ids from the regular gtfs and check how many genes have undefined str
# ----------------------------------------------------------------------------------

while read lab
do
    ~abreschi/Documents/utils/extract.gtf.tags.sh $lab.lncRNAs.gtf gene_id | sort| uniq > $lab.geneIds
done < samples.tsv 

#while read lab
#do
#    cat $lab.lncRNAs.gtf| awk '$3=="transcript"{print $10}'| sort| uniq | sed 's/;//g'| sed 's/"//g'  > $lab.geneIds
#done < samples.tsv 


# summary stats (to see if the numbers correspond to what is on gencode web) and check how many genes have undefined str
while read lab
do
    nb=`cat $lab.geneIds| wc -l`
    echo -e "$lab\t$nb\tRaw"
    un=`/users/rg/jlagarde/julien_utils/extract_locus_coords.pl $lab.lncRNAs.gtf| awk '$6=="."'| wc -l`
    echo -e "$lab\t$un\tUndefined_str"
done < samples.tsv | sed 's/gen/GENCODE/g'| sed 's/refseq/RefSeq/g'|sed 's/fantomCat/FANTOM CAT/g'| sed 's/mitrans/MiTranscriptome/g'| sed 's/bigtrans/BIGTranscriptome/g'| sed 's/noncode/NONCODE/g' | sed 's/cls/GENCODE+/g'| sed 's/pcConf/Protein Coding/g'| sed 's/GENCODE+FL/CLS FL/g'| tee annots.ref.rel.gene.nb.summary.stats.tsv

# join
cat annots.ref.rel.gene.nb.summary.stats.tsv annots.ref.rel.gene.nb.summary.stats.buildLoci.hg38.tsv >  annots+buildLoci.ref.rel.gene.nb.summary.stats.tsv

# plot
echo "
library(ggplot2)
library(scales)
plot <- read.table(\"annots+buildLoci.ref.rel.gene.nb.summary.stats.tsv\", header=F, as.is=T, sep=\"\t\")
colnames(plot)<-c(\"gene\",\"count\",\"typ\")
plot\$gene=factor(plot\$gene, levels=c(\"NONCODE\",\"MiTranscriptome\",\"FANTOM CAT\",\"RefSeq\",\"GENCODE\",\"BIGTranscriptome\",\"GENCODE+\", \"Protein Coding\", \"CLS FL\"))
plot\$typ=factor(plot\$typ, levels=c(\"Raw\",\"buildLoci\",\"Undefined_str\"))
#png(\"annots+buildLoci.ref.rel.gene.nb.summary.stats.png\", bg = \"white\", units=\"in\", width=10, height=8, res=300)
pdf(\"annots+buildLoci.ref.rel.gene.nb.summary.stats.pdf\", bg = \"white\", width=10, height=9)
#setEPS()
#postscript(\"annots+buildLoci.ref.rel.gene.nb.summary.stats.eps\", family=\"serif\", width=10, height=8)
ggplot(data=plot, aes(x=gene, y=count, fill=typ)) + geom_bar(stat=\"identity\", position=\"dodge\") + ylab(\" \")+ theme_bw() + ggtitle(\"\") + scale_y_continuous(limits=c(0, 100000), breaks=seq(0,100000, by=20000))+ scale_fill_manual(values=c(\"#63979d\", \"#8ec85b\", \"hotpink\"))+ xlab(\"\")+ ylab (\"Number of genes\")+theme(axis.text.x = element_text(size = 16, colour = \"black\", angle=45, hjust=1), axis.text.y = element_text(size = 16,  colour = \"black\"), legend.text = element_text(size=18), plot.title = element_text(size = 22), legend.title =element_text(size=0, color=\"white\"), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 18), legend.text.align=0.5, strip.text.x=element_text(size=18), strip.text.y=element_text(size=18))+ ggtitle('') + theme(legend.position=\"top\") + geom_text(aes(x = gene, y = count, label = comma(count), hjust = 1, vjust = 0.5), size=6, position = position_dodge(width = 0.8), angle=90)
dev.off()
" | R --slave

#cp annots+buildLoci.ref.rel.gene.nb.summary.stats.png /users/rg/buszczynska/public_html/gencode/tmp/B/review
cp annots+buildLoci.ref.rel.gene.nb.summary.stats.pdf /users/rg/buszczynska/public_html/gencode/tmp/B/review
#cp annots+buildLoci.ref.rel.gene.nb.summary.stats.eps /users/rg/buszczynska/public_html/gencode/tmp/B/review

# IMPORTANT NOTE: The number of 17,787 loci for NONCODE is calculated if a transcript that constitutes a start or stop of the locus has undefined strand. If the condition is that the at least one transcript has undefined strand then it we have 18,183 loci. For mitrans is the same 11,330.

#0d6598

# --------------------------------------------------------------
# Get the spliced lengths
# --------------------------------------------------------------

# gtf for exons

while read lab
do
    echo $lab
    cat $lab.lncRNAs.hg38.gtf | awk '$3=="exon"' > $lab.exon.lncRNAs.hg38.gtf  
done < samples.tsv

# get the tx ids
while read lab
do
    cat $lab.lncRNAs.hg38.gene.gtf | awk '$3=="exon"{print $12}' | sed 's/;//g'| sed 's/"//g'| sort| uniq > $lab.buildLoci.txIds
    #nb=`cat $lab.buildLoci.txIds| wc -l`
    #echo -e "$lab\t$nb"
done < samples.tsv

while read lab
do
    cat $lab.lncRNAs.gtf | awk '$3=="exon"{print $12}' | sed 's/;//g'| sed 's/"//g'| sort| uniq > $lab.txIds
    nb=`cat $lab.txIds| wc -l`
    echo -e "$lab\t$nb"
done < samples.tsv

#cat cls.lncRNAs.gtf | awk '$3=="exon"{print $12}' | sed 's/;//g'| sed 's/"//g'| sort| uniq > cls.txIds
#cat clsFL.lncRNAs.gtf | awk '$3=="exon"{print $12}' | sed 's/;//g'| sed 's/"//g'| sort| uniq > clsFL.txIds

## get the spliced lengths in fa (hg19 based annotations)

#for lab in bigtrans fantomCat mitrans
#do
#    extract_spliced_transc_seqs.pl $lab.exon.lncRNAs.gtf /users/rg/buszczynska/GENOMES/hg19/H.sapiens.genome.hg19.main.fa > $lab.spliced.lncRNAs.fa
#done

# get the spliced lengths in fa (hg38 based annotations)

#for lab in gen refseq noncode
#do
#    extract_spliced_transc_seqs.pl $lab.exon.lncRNAs.gtf /users/rg/buszczynska/GENOMES/hg38/human.genome.hg38.fa > $lab.spliced.lncRNAs.fa
#done

while read lab
do
    extract_spliced_transc_seqs.pl $lab.exon.lncRNAs.hg38.gtf /users/rg/buszczynska/GENOMES/hg38/human.genome.hg38.fa > $lab.spliced.lncRNAs.hg38.fa
done < samples.tsv
# sumarize the spliced lengths

echo " 
while read lab
do
    FastaToTbl \$lab.spliced.lncRNAs.hg38.fa| while read tx sequ
    do
        length=\$(echo \$sequ| wc -m)
        echo -e \"\$tx\t\$length\t\$lab\tHuman\" 
    done
done < samples.tsv | sed 's/gen/GENCODE/g'| sed 's/refseq/RefSeq/g'|sed 's/fantomCat/FANTOM CAT/g'| sed 's/mitrans/MiTranscriptome/g'| sed 's/bigtrans/BIGTranscriptome/g'| sed 's/noncode/NONCODE/g' | sed 's/cls/GENCODE+/g'| sed 's/pcConf/Protein Coding/g'| sed 's/GENCODE+FL/CLS FL/g' > annots.ref.rel.spliced.tx.length.summary.stats.tsv "  > spliced.len.sh

qsub -cwd -N spliced.len -q rg-el7,long-sl65,short-sl65,short-sl7,long-sl7  -l h_rt=20:00:00 -l virtual_free=30G -pe smp 1 -e spliced.len.log spliced.len.sh


echo "
library(ggplot2)
library(scales)
cbPalette <- c(\"#E41A1C\", \"#377EB8\", \"#4DAF4A\", \"#984EA3\", \"#edde7e\", \"#FF7F00\")
plot <- read.table('annots.ref.rel.spliced.tx.length.summary.stats.tsv', header=F, as.is=T, sep=\"\t\")
colnames(plot)<-c(\"txID\",\"count\", \"typ\", \"class\")
plot\$typ=factor(plot\$typ, levels=c(\"NONCODE\",\"MiTranscriptome\",\"FANTOM CAT\",\"RefSeq\",\"GENCODE\",\"BIGTranscriptome\",\"GENCODE+\", \"Protein Coding\", \"CLS FL\"))
list=c(\"NONCODE\",\"MiTranscriptome\",\"FANTOM CAT\",\"RefSeq\",\"GENCODE\",\"BIGTranscriptome\")
label=matrix(0,length(list),2) 
for (i in 1:length(list))
{
median=median(plot[plot\$typ==list[i], 2],na.rm=T)
len=length(plot[plot\$typ==list[i], 2])
t=paste(\"N =\", prettyNum(len,big.mark=\",\",scientific=FALSE), \"\nMedian =\", prettyNum(round(median,2),big.mark=\",\",scientific=FALSE))
label[i,1]=list[i]
label[i,2]=t
}

label=as.data.frame(label)
row.names(label)=NULL
colnames(label)=c('typ', 'text')
png(\"annots.ref.rel.spliced.tx.length.summary.stats.human.png\", bg = \"white\", units=\"in\", width=6, height=6, res=300)
ggplot(plot, aes(x=count, color=typ)) + geom_density(size=1.5) + scale_color_manual(values=cbPalette)+ scale_x_continuous(limits=c(0, 5000), breaks=seq(0,5000, by=500))+ xlab(\"Spliced transcript length (nt)\") + theme_bw() + ylab (\"Density\")+theme(axis.text.x = element_text(size = 11, colour = \"black\"), axis.text.y = element_text(size = 12,  colour = \"black\"), legend.text = element_text(size=14), plot.title = element_text(size = 22), legend.title =element_text(size=0, color=\"white\"), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14), legend.text.align=0.5, strip.text.x=element_text(size=18), strip.text.y=element_text(size=18))+ ggtitle('')+ theme(legend.position=\"top\")+ facet_wrap(~class)+ geom_text(data=label[1,], aes(x=2300, y=0.0016, label=text), size=5, color=\"#E41A1C\")+ geom_text(data=label[2,], aes(x=2300, y=0.0012, label=text), size=5, color=\"#377EB8\")+ geom_text(data=label[3,], aes(x=2300, y=0.0008, label=text), size=5, color=\"#4DAF4A\")+ geom_text(data=label[4,], aes(x=4300, y=0.0016, label=text), size=5, color=\"#984EA3\")+ geom_text(data=label[5,], aes(x=4300, y=0.0012, label=text), size=5, color=\"#edde7e\")+ geom_text(data=label[6,], aes(x=4300, y=0.0008, label=text), size=5, color=\"#FF7F00\")

dev.off()
" | R --slave

cp annots.ref.rel.spliced.tx.length.summary.stats.human.png /users/rg/buszczynska/public_html/gencode/tmp/B/review

# --------------------------------------------------------------
# Get the isoform per gene stats
# --------------------------------------------------------------

# get nb of isoforms per gene
while read lab
do
    cat $lab.lncRNAs.hg38.gene.gtf| awk 'BEGIN{OFS="\t"}$3=="exon"{print $10,$12}'| sed 's/;//g'| sed 's/"//g'| sort| uniq| awk '{print $1}'| sort|uniq -c| awk 'BEGIN{OFS="\t"}{print $2,$1}' > $lab.iso.nb.per.locus.lncRNAs.tsv
done < samples.tsv

# summary stats
while read lab
do
    for t in $(seq 1 10)
    do
        total=`cat $lab.buildLoci.geneIds| wc -l`
        count=`cat $lab.iso.nb.per.locus.lncRNAs.tsv| awk -v c=$t '$2==c'| wc -l`
        countP=`echo "$count * 100" / $total| bc -l`
        echo -e "$t\t$countP\t$lab\tHuman"
    done
done < samples.tsv | sed 's/gen/GENCODE/g'| sed 's/refseq/RefSeq/g'|sed 's/fantomCat/FANTOM CAT/g'| sed 's/mitrans/MiTranscriptome/g'| sed 's/bigtrans/BIGTranscriptome/g'| sed 's/noncode/NONCODE/g' | sed 's/cls/GENCODE+/g'| sed 's/pcConf/Protein Coding/g'| sed 's/GENCODE+FL/CLS FL/g' > annots.ref.rel.iso.nb.per.locus.summary.stats.tsv

while read lab
do
    total=`cat $lab.buildLoci.geneIds| wc -l`
    count=`cat $lab.iso.nb.per.locus.lncRNAs.tsv| awk '$2>10'| wc -l`
    countP=`echo "$count * 100" / $total| bc -l`
    echo -e ">10\t$countP\t$lab\tHuman"
done < samples.tsv | sed 's/gen/GENCODE/g'| sed 's/refseq/RefSeq/g'|sed 's/fantomCat/FANTOM CAT/g'| sed 's/mitrans/MiTranscriptome/g'| sed 's/bigtrans/BIGTranscriptome/g'| sed 's/noncode/NONCODE/g' | sed 's/cls/GENCODE+/g'| sed 's/pcConf/Protein Coding/g'| sed 's/GENCODE+FL/CLS FL/g' >> annots.ref.rel.iso.nb.per.locus.summary.stats.tsv

echo "
library(ggplot2)
library(RColorBrewer)
cbPalette <- c(\"#a6761d\", \"#e6ab02\", \"#FF7F00\", \"#984EA3\", \"#4DAF4A\", \"#377EB8\", \"#E41A1C\",\"#999999\", \"#e7298a\")
plot <- read.table('annots.ref.rel.iso.nb.per.locus.summary.stats.tsv', header=F, as.is=T, sep=\"\t\")
colnames(plot)<-c(\"labExpId\",\"count\",\"category\", \"class\")
plot\$labExpId=factor(plot\$labExpId, levels=c(\"1\",\"2\",\"3\",\"4\",\"5\",\"6\",\"7\",\"8\",\"9\",\"10\",\">10\"))
plot\$category=factor(plot\$category, levels=c(\"NONCODE\",\"MiTranscriptome\",\"FANTOM CAT\",\"RefSeq\",\"GENCODE\",\"BIGTranscriptome\",\"GENCODE+\", \"Protein Coding\", \"CLS FL\"))
png(\"annots.ref.rel.iso.nb.per.locus.summary.stats.human.png\", bg = \"white\", units=\"in\", width=10, height=6, res=300)
ggplot(plot, aes(x=labExpId, y=count, fill=category)) + geom_bar(stat=\"identity\", position=position_dodge()) + scale_fill_manual(values=cbPalette) +ylim(0,100)+ xlab(\"Number of isoform(s)\")+ theme_bw()+ ylab (\"Percent of genes [%]\")+theme(axis.text.x = element_text(size = 11, colour = \"black\"), axis.text.y = element_text(size = 12,  colour = \"black\"), legend.text = element_text(size=14), plot.title = element_text(size = 22), legend.title =element_text(size=0, color=\"white\"), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14), legend.text.align=0.5, strip.text.x=element_text(size=18), strip.text.y=element_text(size=18))+ ggtitle('')+ theme(legend.position=\"top\")+ facet_wrap(~category)
dev.off()
" | R --slave

cp annots.ref.rel.iso.nb.per.locus.summary.stats.human.png /users/rg/buszczynska/public_html/gencode/tmp/B/review

# -------------------------------------------------------------- 
# Exons per transcript
# --------------------------------------------------------------

# summary stats
while read lab
do
    for t in $(seq 1 10)
    do
       total=`cat $lab.buildLoci.txIds| wc -l`
       count=`cat $lab.exon.lncRNAs.hg38.gtf| awk '{print $12}'| sort| uniq -c| awk -v c=$t '$1==c'| wc -l`
       countP=`echo "$count * 100" / $total| bc -l`
       echo -e "$lab\t$t\t$countP\tHuman"
    done 
done < samples.tsv | sed 's/gen/GENCODE/g'| sed 's/refseq/RefSeq/g'|sed 's/fantomCat/FANTOM CAT/g'| sed 's/mitrans/MiTranscriptome/g'| sed 's/bigtrans/BIGTranscriptome/g'| sed 's/noncode/NONCODE/g' | sed 's/cls/GENCODE+/g'| sed 's/pcConf/Protein Coding/g'| sed 's/GENCODE+FL/CLS FL/g' > annots.ref.rel.iso.nb.exons.summary.stats.tsv

# for median
while read lab
do
   cat $lab.exon.lncRNAs.hg38.gtf| awk '{print $12}'| sort| uniq -c| awk -v c=$lab 'BEGIN{OFS="\t"}{print $2,$1,c,"Human"}'
done < samples.tsv | sed 's/;//g'| sed 's/"//g'| sed 's/gen/GENCODE/g'| sed 's/refseq/RefSeq/g'|sed 's/fantomCat/FANTOM CAT/g'| sed 's/mitrans/MiTranscriptome/g'| sed 's/bigtrans/BIGTranscriptome/g'| sed 's/noncode/NONCODE/g' | sed 's/cls/GENCODE+/g'| sed 's/pcConf/Protein Coding/g'| sed 's/GENCODE+FL/CLS FL/g' > annots.ref.rel.iso.nb.exons.median.tsv

# plot

echo "
library(ggplot2)
library(RColorBrewer)
cbPalette <- c(\"#E41A1C\", \"#377EB8\", \"#4DAF4A\", \"#984EA3\", \"#e6ab02\", \"#FF7F00\",\"#a6761d\",\"#999999\", \"#e7298a\")
plot <- read.table('annots.ref.rel.iso.nb.exons.summary.stats.tsv', header=F, as.is=T, sep=\"\t\")
colnames(plot)<-c(\"groupId\",\"labExpId\",\"count\",\"category\")
plot\$labExpId=factor(plot\$labExpId, levels=c(\"1\",\"2\",\"3\",\"4\",\"5\",\"6\",\"7\",\"8\",\"9\",\"10\"))
plot\$groupId=factor(plot\$groupId, levels=c(\"NONCODE\",\"MiTranscriptome\",\"FANTOM CAT\",\"RefSeq\",\"GENCODE\",\"BIGTranscriptome\",\"GENCODE+\", \"Protein Coding\", \"CLS FL\"))
list=c(\"NONCODE\",\"MiTranscriptome\",\"FANTOM CAT\",\"RefSeq\",\"GENCODE\",\"BIGTranscriptome\",\"GENCODE+\", \"Protein Coding\", \"CLS FL\")
support=read.table('annots.ref.rel.iso.nb.exons.median.tsv', header=F, as.is=T, sep=\"\t\")
colnames(support)<-c(\"txId\",\"exon\",\"gen\",\"org\")
label=matrix(0,length(list),2) 
for (i in 1:length(list))
{
median=median(support[support\$gen==list[i], 2],na.rm=T)
len=length(support[support\$gen==list[i], 2])
t=paste(\"N =\", prettyNum(len,big.mark=\",\",scientific=FALSE), \"\nMedian =\", prettyNum(round(median,2),big.mark=\",\",scientific=FALSE))
label[i,1]=list[i]
label[i,2]=t
}

label=as.data.frame(label)
row.names(label)=NULL
colnames(label)=c('groupId', 'text')

png(\"annots.ref.rel.iso.nb.exons.summary.stats.human.png\", bg = \"white\", units=\"in\", width=10, height=8, res=300)
ggplot(plot, aes(x=labExpId, y=count, fill=groupId)) + geom_bar(stat=\"identity\", position=position_dodge())+ylim(0,100) + scale_fill_manual(values=cbPalette) + xlab(\"Number of exons\")+ theme_bw()+ ylab (\"Percent of transcripts (%)\")+theme(axis.text.x = element_text(size = 11, colour = \"black\"), axis.text.y = element_text(size = 12,  colour = \"black\"), legend.text = element_text(size=14), plot.title = element_text(size = 22), legend.title =element_text(size=0, color=\"white\"), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14), legend.text.align=0.5, strip.text.x=element_text(size=18), strip.text.y=element_text(size=18))+ ggtitle('')+ theme(legend.position=\"top\")+ facet_wrap(~ category)+ geom_text(data=label[1,], aes(x=8, y=70, label=text), size=5, color=\"#E41A1C\")+ geom_text(data=label[2,], aes(x=8, y=70, label=text), size=5, color=\"#377EB8\")+ geom_text(data=label[3,], aes(x=8, y=70, label=text), size=5, color=\"#4DAF4A\")+ geom_text(data=label[4,], aes(x=8, y=70, label=text), size=5, color=\"#984EA3\")+ geom_text(data=label[5,], aes(x=8, y=70, label=text), size=5, color=\"#e6ab02\")+ geom_text(data=label[6,], aes(x=8, y=70, label=text), size=5, color=\"#FF7F00\")+ geom_text(data=label[7,], aes(x=8, y=70, label=text), size=5, color=\"#a6761d\")+ geom_text(data=label[8,], aes(x=8, y=70, label=text), size=5, color=\"#999999\")+ geom_text(data=label[9,], aes(x=8, y=70, label=text), size=5, color=\"#e7298a\")+facet_wrap(~ groupId)
dev.off()
" | R --slave

cp annots.ref.rel.iso.nb.exons.summary.stats.human.png /users/rg/buszczynska/public_html/gencode/tmp/B/review

# --------------------------------------------
# GFFCOMPARE
# --------------------------------------------

#while read lab
#do 
#     cat $lab.lncRNAs.gtf| awk '{if($3=="transcript" || $3=="exon") print $0}'| sort| uniq > $lab.lncRNAs.tx.gtf
#done < samples.tsv 


while read lab
do 
    echo $lab
    #/users/rg/buszczynska/bin/gffcompare/gffcompare -N -M -r gen.lncRNAs.tx.gtf $lab.lncRNAs.tx.gtf -o $lab.genOverlap
    /users/rg/buszczynska/bin/gffcompare/gffcompare -N -M -r gen.lncRNAs.hg38.gtf $lab.lncRNAs.hg38.gtf -o $lab.genOverlap
done < samples.tsv 

#reassemble all data

while read lab
do
    echo -e "annot\ttyp\tCategory\tvalue" > annots.genOverlap.allLevels.stats.tsv
    for level in `echo Baselevel Exonlevel Intronchainlevel Intronlevel Locuslevel Transcriptlevel`
    do
        file="/users/rg/buszczynska/Projects/review/human/annot.comp/$lab.genOverlap"
        #echo $level
        Sn=`cat $file |grep "level:" |sed 's/ //g'| sed 's/:/\t/'|sed 's/|$//'|sed 's/|/\t/g' | awk -v l=$level '$1==l' |cut -f2`
        Sp=`cat $file |grep "level:" |sed 's/ //g'| sed 's/:/\t/'|sed 's/|$//'|sed 's/|/\t/g' | awk -v l=$level '$1==l' |cut -f3`
        echo -e "$lab\t$level\tSn\t$Sn\n$lab\t$level\tSp\t$Sp"  
    done
done < samples.tsv | sort| uniq | sed 's/gen/GENCODE/g'| sed 's/refseq/RefSeq/g'|sed 's/fantomCat/FANTOM CAT/g'| sed 's/mitrans/MiTranscriptome/g'| sed 's/bigtrans/BIGTranscriptome/g'| sed 's/noncode/NONCODE/g' | sed 's/cls/GENCODE+/g'| sed 's/pcConf/Protein Coding/g'| sed 's/GENCODE+FL/CLS FL/g'| sed 's/Sn/Sensitivity/g'| sed 's/Sp/Precision/g'| sed 's/level//g' >> annots.genOverlap.allLevels.stats.tsv 

# plot
echo "
library(ggplot2)
library(scales)
cbPalette <- c(\"#E41A1C\", \"#377EB8\", \"#4DAF4A\", \"#984EA3\", \"#e6ab02\", \"#FF7F00\",\"#a6761d\",\"#999999\", \"#e7298a\")
dat <- read.table(\"annots.genOverlap.allLevels.stats.tsv\", header=T, as.is=T, sep=\"\t\")
dat\$annot=factor(dat\$annot, levels=c(\"NONCODE\",\"MiTranscriptome\",\"FANTOM CAT\",\"RefSeq\",\"GENCODE\",\"BIGTranscriptome\",\"GENCODE+\", \"Protein Coding\", \"CLS FL\"))
png(\"annots.genOverlap.allLevels.stats.human.png\", bg = \"white\", units=\"in\", width=8, height=6, res=300)
ggplot(dat, aes(x=typ, y=value)) + geom_point(aes(color=annot, shape=Category), size=5, alpha=0.7)+xlab(\"\")+ scale_colour_manual (values=cbPalette, name=\"Accuracy measure\", breaks=c(\"Sn\", \"Sp\"), labels=c(\"Sensitivity\",\"Precision\"))+ theme_bw()+ ylab(\"Sensitivity | Precision (%)\") + scale_y_continuous() + expand_limits(y=c(0,100))+ theme(legend.position=\"bottom\", legend.background = element_rect(fill=\"gray90\", size=.5, linetype=\"dotted\"))+  guides(color = guide_legend(nrow=1))+facet_wrap(~ annot)+theme(axis.text.x = element_text(size = 11, colour = \"black\", angle=45, hjust = 1), axis.text.y = element_text(size = 12,  colour = \"black\"), legend.text = element_text(size=14), plot.title = element_text(size = 22), legend.title =element_text(size=0, color=\"white\"), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14), legend.text.align=0.5, strip.text.x=element_text(size=18), strip.text.y=element_text(size=18))+ ggtitle('')+ theme(legend.position=\"bottom\")
dev.off()
" | R --slave

cp annots.genOverlap.allLevels.stats.human.png /users/rg/buszczynska/public_html/gencode/tmp/B/review

# --------------------------------------------
# Custom tracks
# --------------------------------------------

while read lab
do
   { echo -e "track name=$lab\tdescription=\"$lab\"\tcolor=38,173,251"; cat $lab.hg38.bed12 ; } > /users/rg/buszczynska/public_html/gencode/custom.tracks/$lab.hg38.bed12
done < samples.tsv 

# ------------------------------------------------------------------------------
# FANTOM CAT hg38
# ------------------------------------------------------------------------------

# FANTOM CAT transcripts
{ echo -e "track name=fantomCat.lncRNAs\tdescription=\"fantomCat.lncRNAs\"\tcolor=2,182,121"; /users/rg/jlagarde/julien_utils/extract_locus_coords.pl fantomCat.lncRNAs.hg38.gtf ; } > fantomCat.lncRNAs.hg38.bed

# BUILDLOCI
{ echo -e "track name=fantomCat.lncRNAs.buildLoci\tdescription=\"fantomCat.lncRNAs.buildLoci\"\tcolor=238,115,164"; /users/rg/jlagarde/julien_utils/extract_locus_coords.pl fantomCat.lncRNAs.hg38.gene.gtf ; } > fantomCat.lncRNAs.loci.hg38.buildLoci.bed

# loci orginal fantom
/users/rg/jlagarde/julien_utils/extract_locus_coords.pl fantomCat.lncRNAs.gtf > fantomCat.lncRNAs.hg19.loci.bed

# lift up
liftOver fantomCat.lncRNAs.hg19.loci.bed /users/rg/buszczynska/liftOver/hg19ToHg38.over.chain fantomCat.lncRNAs.hg38.loci.bed fantomCat.lncRNAs.hg38.loci.unMapped

# FANTOM CAT loci
{ echo -e "track name=fantomCat.loci\tdescription=\"fantomCat.loci\"\tcolor=115,226,238"; cat fantomCat.lncRNAs.hg38.loci.bed ; } > fantomCat.lncRNAs.loci.hg38.bed


# ----------------------------------------
# Join the old and new gene ids
# ----------------------------------------

# gene, txids for old
cat fantomCat.lncRNAs.gtf| awk 'BEGIN{OFS="\t"}$3=="transcript"{print $10,$12}'| sed s'/;//'g| sed 's/"//g' > fantomCat.hg19.gene.txIDs

# gene, txids for old
cat fantomCat.lncRNAs.hg38.gene.gtf| awk 'BEGIN{OFS="\t"}{print $10,$12}'| sed s'/;//'g| sed 's/"//g'| sort| uniq > fantomCat.hg38.buildLoci.gene.txIDs

# join 

join.py -a fantomCat.hg19.gene.txIDs -b fantomCat.hg38.buildLoci.gene.txIDs -x 2 -y 2 | cut -f1,3| sort| uniq| cut -f2| sort| uniq -c | awk '$1>1'| wc -l
join.py -a fantomCat.hg19.gene.txIDs -b fantomCat.hg38.buildLoci.gene.txIDs -x 2 -y 2 | cut -f1,3| sort| uniq| cut -f1| sort| uniq -c | awk '$1>1'| wc -l

# ------------------------------------------------------------------------------
# CLS
# ------------------------------------------------------------------------------

# FANTOM CAT transcripts
{ echo -e "track name=cls.lncRNAs\tdescription=\"cls.lncRNAs\"\tcolor=2,182,121"; /users/rg/jlagarde/julien_utils/extract_locus_coords.pl cls.lncRNAs.hg38.gtf ; } > /users/rg/buszczynska/public_html/gencode/custom.tracks/cls.lncRNAs.hg38.bed

# BUILDLOCI
{ echo -e "track name=cls.lncRNAs.buildLoci\tdescription=\"cls.lncRNAs.buildLoci\"\tcolor=238,115,164"; /users/rg/jlagarde/julien_utils/extract_locus_coords.pl cls.lncRNAs.hg38.gene.gtf ; } > /users/rg/buszczynska/public_html/gencode/custom.tracks/cls.lncRNAs.loci.hg38.buildLoci.bed

# ------------------------------------------------------------------------------
# NONCODE hg38
# ------------------------------------------------------------------------------

# NONCODE transcripts
{ echo -e "track name=noncode.lncRNAs\tdescription=\"noncode.lncRNAs\"\tcolor=32,52,164"; /users/rg/jlagarde/julien_utils/extract_locus_coords.pl noncode.lncRNAs.hg38.gtf ; } > /users/rg/buszczynska/public_html/gencode/custom.tracks/noncode.lncRNAs.hg38.bed

# BUILDLOCI
{ echo -e "track name=noncode.lncRNAs.buildLoci\tdescription=\"noncode.lncRNAs.buildLoci\"\tcolor=243,226,30"; /users/rg/jlagarde/julien_utils/extract_locus_coords.pl noncode.lncRNAs.hg38.gene.gtf ; } > /users/rg/buszczynska/public_html/gencode/custom.tracks/noncode.lncRNAs.loci.hg38.buildLoci.bed

# ORginal loci
{ echo -e "track name=noncode.loci\tdescription=\"noncode.loci\"\tcolor=115,226,238"; /users/rg/jlagarde/julien_utils/extract_locus_coords.pl noncode.lncRNAs.gtf ; } > /users/rg/buszczynska/public_html/gencode/custom.tracks/noncode.lncRNAs.loci.hg38.bed

# ----------------------------------------
# Join the old and new gene ids
# ----------------------------------------

# gene, txids for old
cat noncode.lncRNAs.gtf| awk 'BEGIN{OFS="\t"}$3=="transcript"{print $10,$12}'| sed s'/;//'g| sed 's/"//g' > noncode.hg19.gene.txIDs

# gene, txids for old
cat noncode.lncRNAs.hg38.gene.gtf| awk 'BEGIN{OFS="\t"}{print $10,$12}'| sed s'/;//'g| sed 's/"//g'| sort| uniq > noncode.hg38.buildLoci.gene.txIDs

# join 
join.py -a noncode.hg38.gene.txIDs -b noncode.hg38.buildLoci.gene.txIDs -x 2 -y 2 | cut -f1,3| sort| uniq| cut -f1| sort| uniq -c | awk '$1>1'| head

# ------------------------------------------------------------------------------
# MITRANS hg38
# ------------------------------------------------------------------------------

# miTranscriptome transcripts
{ echo -e "track name=mitrans.lncRNAs\tdescription=\"mitrans.lncRNAs\"\tcolor=26,161,245"; /users/rg/jlagarde/julien_utils/extract_locus_coords.pl mitrans.lncRNAs.hg38.gtf ; } > /users/rg/buszczynska/public_html/gencode/custom.tracks/mitrans.lncRNAs.hg38.bed

# BUILDLOCI
{ echo -e "track name=mitrans.lncRNAs.buildLoci\tdescription=\"mitrans.lncRNAs.buildLoci\"\tcolor=181,151,202"; /users/rg/jlagarde/julien_utils/extract_locus_coords.pl mitrans.lncRNAs.hg38.gene.gtf ; } > /users/rg/buszczynska/public_html/gencode/custom.tracks/mitrans.lncRNAs.loci.hg38.buildLoci.bed

# loci orginal fantom
/users/rg/jlagarde/julien_utils/extract_locus_coords.pl mitrans.lncRNAs.gtf > mitrans.lncRNAs.hg19.loci.bed

# lift up
liftOver mitrans.lncRNAs.hg19.loci.bed /users/rg/buszczynska/liftOver/hg19ToHg38.over.chain mitrans.lncRNAs.hg38.loci.bed mitrans.lncRNAs.hg38.loci.unMapped

# miTranscriptome loci
{ echo -e "track name=mitrans.loci\tdescription=\"mitrans.loci\"\tcolor=154,243,30"; cat mitrans.lncRNAs.hg38.loci.bed ; } > /users/rg/buszczynska/public_html/gencode/custom.tracks/mitrans.lncRNAs.loci.hg38.bed


# ----------------------------------------
# Join the old and new gene ids
# ----------------------------------------

# gene, txids for old
cat mitrans.lncRNAs.gtf| awk 'BEGIN{OFS="\t"}$3=="transcript"{print $10,$12}'| sed s'/;//'g| sed 's/"//g' > mitrans.hg19.gene.txIDs

# gene, txids for old
cat mitrans.lncRNAs.hg38.gene.gtf| awk 'BEGIN{OFS="\t"}{print $10,$12}'| sed s'/;//'g| sed 's/"//g'| sort| uniq > mitrans.hg38.buildLoci.gene.txIDs

# join 
join.py -a mitrans.hg19.gene.txIDs -b mitrans.hg38.buildLoci.gene.txIDs -x 2 -y 2 | cut -f1,3| sort| uniq| cut -f1| sort| uniq -c | awk '$1>1'| head

# ------------------------------------------------------------------------------
# FANTOM CAT hg19
# ------------------------------------------------------------------------------

#FANTOM CAT txs 
#awk -f /no_backup_isis/rg/0ld_users/sdjebali/Awk/bed12fields2gff.awk fantomCat.bed12 > fantomCat.lncRNAs.hg19.gff

#{ echo -e "track name=fantomCat.lncRNAs\tdescription=\"fantomCat.lncRNAs\"\tcolor=2,182,121"; /users/rg/jlagarde/julien_utils/extract_locus_coords.pl fantomCat.lncRNAs.hg19.gff ; } > fantomCat.lncRNAs.hg19.bed


#FANTOM CAT loci 
#{ echo -e "track name=fantomCat.loci\tdescription=\"fantomCat.loci\"\tcolor=115,226,238"; /users/rg/jlagarde/julien_utils/extract_locus_coords.pl fantomCat.lncRNAs.gtf ; } > fantomCat.lncRNAs.loci.hg19.bed

# FANTOM CAT FOR BULID LOCI
# get gff witout gene records (don't contain transcript_ids)
#/users/rg/jlagarde/julien_utils/gff2gff.pl fantomCat.lncRNAs.gtf| awk '$3=="exon"' > fantomCat.lncRNAs.hg19.gtf

# buildLoci
#~jlagarde/bin/bedtools2/bin/bedtools intersect -s -wao -a fantomCat.lncRNAs.hg19.gtf -b fantomCat.lncRNAs.hg19.gtf | ~jlagarde/julien_utils/buildLoci.pl - > fantomCat.lncRNAs.hg19.gene.tmp.gtf
#/users/rg/jlagarde/julien_utils/gff2gff.pl fantomCat.lncRNAs.hg19.gene.tmp.gtf > fantomCat.lncRNAs.hg19.gene.gtf

# get custom track for genes
#{ echo -e "track name=fantomCat.lncRNAs.buildLoci\tdescription=\"fantomCat.lncRNAs.buildLoci\"\tcolor=238,115,164"; /users/rg/jlagarde/julien_utils/extract_locus_coords.pl fantomCat.lncRNAs.hg19.gene.gtf ; } > fantomCat.lncRNAs.loci.hg19.buildLoci.bed

 