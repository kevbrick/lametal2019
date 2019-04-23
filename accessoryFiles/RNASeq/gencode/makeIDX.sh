perl ../../scripts/catGTFgeneNames.pl <gencode.vM20.annotation.gtf >catNames.gtf
../../scripts/gtf_to_fasta catNames.gtf ../../genomeFiles/mm10_genome.fa gencode.vM20.transcripts.fa
perl -pi -e 's/^\>\d+\s+(\S+)\s+chr.+/>$1/' gencode.vM20.transcripts.fa

kallisto index -i gencode.vM20.kallistoIDX --make-unique gencode.vM20.transcripts.fa
