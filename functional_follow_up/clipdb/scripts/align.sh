#!/bin/bash
set -euxo pipefail
PROJECT=/fsimb/groups/imb-buttergr/=EVOREG/rnalib/functional_follow_up/clipdb
#REF=/fsimb/common/genomes/saccharomyces_cerevisiae/ucsc/saccer3/canonical/index/bowtie/genome
REF=/fsimb/groups/imb-buttergr/=EVOREG/rnalib/functional_follow_up/clipdb/ref/sacCer3
BOWTIE=/fsimb/groups/imb-bioinfocf/common-tools/dependencies/bowtie/1.1.2/bowtie
SAMTOOLS=/fsimb/groups/imb-bioinfocf/common-tools/dependencies/samtools/1.3.1/samtools
MACHINES=(imbc1 imbc4 imbc5)

for f in ${PROJECT}/rawdata/*.trimmed.fastq.gz; do
  F=$(basename ${f%.trimmed.fastq.gz}.bam)
#  echo "zcat $f | $BOWTIE -q -p 4 -S -v 1 -m 1 -y --best --strata $REF - | $SAMTOOLS view -bShu |  $SAMTOOLS sort -@4 -T /tmp/${F}.tmp > ${PROJECT}/mapped/${F}" | bsub -J $F -W5:00 -n4 -R'span[ptile=4]' -app Reserve2G -o $F.out -e $F.err
  # new set of bowtie options, as recommended by paralyzer (see webpage and ../paralyzer/callpeaks.sh)
#  echo "zcat $f | $BOWTIE -q -p 8 -S -v 2 -m 10 --best --strata $REF - | $SAMTOOLS view -bShu |  $SAMTOOLS sort -@4 -T /tmp/${F}.tmp > ${PROJECT}/mapped/${F}" | bsub -J $F -W5:00 -n4 -R'span[ptile=4]' -app Reserve2G -o $F.out -e $F.err
  # bsub is not available anymore!
  zcat $f | $BOWTIE -q -p 8 -S -v 2 -m 10 --best --strata $REF - | $SAMTOOLS view -bShu |  $SAMTOOLS sort -@4 -T /tmp/${F}.tmp > ${PROJECT}/mapped/${F} && samtools index ${PROJECT}/mapped/${F}
done
