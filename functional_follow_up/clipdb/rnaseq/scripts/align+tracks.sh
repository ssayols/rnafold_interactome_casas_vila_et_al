#!/bin/bash
set -euxo pipefail

export PROJECT=/fsimb/groups/imb-buttergr/=EVOREG/rnalib/functional_follow_up/clipdb/rnaseq
#export REF=/fsimb/common/genomes/saccharomyces_cerevisiae/ucsc/saccer3/canonical/index/bowtie/genome
export REF=/fsimb/groups/imb-buttergr/=EVOREG/rnalib/functional_follow_up/clipdb/ref/sacCer3
export CHRSIZES=${REF}.chrom.sizes
export BOWTIE=/fsimb/groups/imb-bioinfocf/common-tools/dependencies/bowtie/1.1.2/bowtie
export SAMTOOLS=/fsimb/groups/imb-bioinfocf/common-tools/dependencies/samtools/1.3.1/samtools

# tracks
function track {
  BAM=$1

  ml picard bedtools/2.25.0 kentUtils
  TOTAL_MAPPED=$(picard BamIndexStats I=${BAM} | awk '{ SUM += $5 } END { print SUM }')
  SCALE=$(echo "1000000/$TOTAL_MAPPED" | bc -l)
  genomeCoverageBed -bg -split -scale ${SCALE} -ibam $BAM -g $CHRSIZES | sort -k1,1 -k2,2n > ${BAM%.bam}_scaled.bedgraph
  bedGraphToBigWig ${BAM%.bam}_scaled.bedgraph $CHRSIZES ${BAM%.bam}_scaled.bw
  rm ${BAM%.bam}_scaled.bedgraph
  mv ${BAM%.bam}_scaled.bw ${PROJECT}/tracks
}
export -f track

# align + tracks
function align {
  SRR=$1
  F=$2

  ml sratoolkit
  fastq-dump --gzip -Z -A ${SRR} | \
    tee ${PROJECT}/rawdata/${F}.fastq.gz | \
    gzip -cd | \
    $BOWTIE -q -p 4 -S -v 2 -m 1 --best --strata $REF - | \
    $SAMTOOLS view -bShu | \
    $SAMTOOLS sort -@4 -T /tmp/${F}.bam > ${PROJECT}/mapped/${F}.bam
  samtools index ${PROJECT}/mapped/${F}.bam
  track ${PROJECT}/mapped/${F}.bam
}
export -f align

# run
parallel -j4 --xapply align {1} {2} ::: SRR1066657 SRR1066658 SRR1066659 SRR1066660 ::: WT_NR_A_GSM1299413 WT_NR_B_GSM1299414 WT_CR_A_GSM1299415 WT_CR_B_GSM1299416

