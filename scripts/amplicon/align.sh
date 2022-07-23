#!/bin/bash -e
#
# Script for aligning amplicon reads to SARS-CoV-2 reference
ALIGNER="$1"
READS="$2"
OUTPUT_SAM="$3"
CORES=`grep -c ^processor /proc/cpuinfo`

# Use minimap2 as aligner (default)
if [[ "${ALIGNER}" == "minimap2" ]]; then
  minimap2 -ax sr "wuhCor1.fa" "${READS}" > $OUTPUT_SAM
fi

#TODO: Add ViralMSA
