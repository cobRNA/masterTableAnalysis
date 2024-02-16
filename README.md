# masterTableAnalysis
Short scripts and notebooks written mostly in Python and R to analyze Gencode - CLS masterTable






# Clustering
bedtools merge -d 5 -s -c 4,6 -o distinct -i bigtrans.5.bed.vsCage.fantom.bedtsv > merged.bigtrans_5


# Anchoring
./anchorTranscriptsEnds.pl bigtrans.hg38.gff bigtrans.tx.cage bigtrans.tx.polyA merged.bigtrans_5 merged.bigtrans_3 > test.gff
