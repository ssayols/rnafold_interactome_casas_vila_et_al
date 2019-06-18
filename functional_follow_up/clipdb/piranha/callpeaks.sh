#!/bin/bash
# as in clipdb (postar2), use 2 methods:

#
# piranha
#
../../../bin/piranha/Piranha -s -b 50 -d Poisson ../mapped/PAB1_GSM1442553.bam | sort -k1,1 -k2,2n > PAB1_GSM1442553.peaks.bed
../../../bin/piranha/Piranha -s -b 50 -d Poisson ../mapped/YRA1_GSM1442559.bam | sort -k1,1 -k2,2n > YRA1_GSM1442559.peaks.bed

