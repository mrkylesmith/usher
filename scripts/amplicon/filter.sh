#!/bin/bash
#
# Combine a directory of many FASTQ files into one, remove duplicates.

INPUT_DIR="$1"  # Directory containing many FASTQ files to combine 
OUTPUT_DIR="$2"
LOGGING="$OUTPUT_DIR/log"
TIMING="$OUTPUT_DIR/log/runtime_logging"
CORES=`grep -c ^processor /proc/cpuinfo`
START=$SECONDS
mkdir -p $OUTPUT_DIR/log

#######################################
# Determine if the input file 
# is compressed or not
# Arguments:
#  Input file
# Outputs:
#  Return true if file is compressed (.gz),
#  false if not compressed
#######################################
is_compressed() {
	if [[ $1 =~ \.gz$ ]]; then
    return true
  else
    return false
  fi
}

#######################################
# Get number of reads in a FASTQ file
# Arguments:
#  Input FASTQ file 
# Outputs:
#   Number of reads in that FASTQ file
#######################################
num_reads() {
	if [[ $1 =~ \.gz$ ]]; then
    lines=`zcat $1 | wc -l`
  else
    lines=`wc -l < $1`
  fi
  reads=$(($lines / 4))
	return $reads
}

# Combine all FASTQ files in the given directory into one large FASTQ file
if [[ -d $INPUT_DIR ]]; then
  for file in $INPUT_DIR/*; do
      echo $file
      cat $file >> combined_amplicons.fastq
  done
elif [[ -f $INPUT_DIR ]];then
	      if [[ $INPUT_DIR =~ \.gz$ ]]; then
            echo "Single compressed FASTQ file given as input: $INPUT_DIR"
				    cp $INPUT_DIR combined_amplicons.fastq.gz
        else
            echo "Single uncompressed FASTQ file given as input: $INPUT_DIR"
            cp $INPUT_DIR combined_amplicons.fastq
				fi
else
				echo "ERROR: Check FASTQ input file type"
fi

# Deduplicate FASTQ reads by sequence content
# Information on deduplicated reads (if any) will be logged in $OUTPUT_DIR/log directory
seqkit rmdup --by-seq --ignore-case -d $LOGGING/removed_duplicates.fastq -D $LOGGING/duplicated.detail.txt combined_amplicons.fastq > amplicons.dedup.fastq

# Log runtime
END=$(( SECONDS - START ))
echo "Time to remove FASTQ duplicates: $END seconds" >> $TIMING

# Move aggregated and deduplicated amplicon reads FASTQ file into data directory
mv combined_amplicons.fastq $OUTPUT_DIR/data
mv amplicons.dedup.fastq $OUTPUT_DIR/data
