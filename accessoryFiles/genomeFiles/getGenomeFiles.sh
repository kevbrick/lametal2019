#!/bin/bash
wget http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/chromFa.tar.gz
tar -xzvf chromFa.tar.gz

cat chr??.fa chr?.fa >mm10_genome.fa
cat chr??.fa chr?.fa KMetStatPanel.fa >mm10_KmetStat_genome.fa

samtools faidx mm10_genome.fa
samtools faidx mm10_KmetStat_genome.fa

rm chr*fa
rm chromFa.tar.gz

