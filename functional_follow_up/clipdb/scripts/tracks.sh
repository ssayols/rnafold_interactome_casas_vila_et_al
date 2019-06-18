#!/bin/bash
PROJECT=/fsimb/groups/imb-buttergr/=EVOREG/rnalib/functional_follow_up/clipdb
CHRSIZES=${PROJECT}/ref/sacCer3.chrom.sizes

function process {
    BAM=$1
    ml picard bedtools/2.25.0 kentUtils
    TOTAL_MAPPED=$(picard BamIndexStats I=${BAM} | awk '{ SUM += $5 } END { print SUM }')
    SCALE=$(echo "1000000/$TOTAL_MAPPED" | bc -l)
    genomeCoverageBed -bg -split -scale ${SCALE} -ibam $BAM -g $CHRSIZES | sort -k1,1 -k2,2n > ${BAM%.bam}_scaled.bedgraph
    bedGraphToBigWig ${BAM%.bam}_scaled.bedgraph $CHRSIZES ${BAM%.bam}_scaled.bw
    rm ${BAM%.bam}_scaled.bedgraph
    mv ${BAM%.bam}_scaled.bw ${PROJECT}/tracks
}

export -f process
for f in ${PROJECT}/mapped/*.bam; do
#  bsub -n1 `process "$f"`
  echo $f
  process "$f"
done

