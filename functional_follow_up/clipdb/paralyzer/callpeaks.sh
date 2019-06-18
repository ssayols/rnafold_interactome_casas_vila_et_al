#!/bin/bash
# as in clipdb (postar2), use 2 methods:

#
# PARalyzer, par-clip specific
#
# needs first to align with bowtie to the 2bit genome (ucsc)
# wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa
# wget ftp://hgdownload.soe.ucsc.edu/goldenPath/currentGenomes/Saccharomyces_cerevisiae/bigZips/sacCer3.2bit
# ./twoBitToFa sacCer3.2bit sacCer3.fa
# bowtie-build --threads 16 sacCer3.fa sacCer3
ml bowtie samtools 
#samtools view -h ../mapped/PAB1_GSM1442553.bam | awk '$3 != "*"' > PAB1_GSM1442553.sam
zcat ../rawdata/PAB1_GSM1442553.trimmed.fastq.gz | bowtie -q -v 2 -m 10 --best --strata ../ref/sacCer3 - > PAB1_GSM1442553.bowtie.aln
java -Xmx16G -cp /fsimb/groups/imb-buttergr/=EVOREG/rnalib/bin/PARalyzer_v1_5/src/ PARalyze PAB1_PARalyzer_Parameters.ini

