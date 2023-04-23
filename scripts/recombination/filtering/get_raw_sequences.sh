#!/bin/bash -x
# 
# Script to retrieve fasta sequences for all nodes in descendants.txt

# Command line args
TREE_DATE=""       # -d
INPUT=""           # -i
OUTPUT=""          # -o
BUCKET_ID=""       # -b
DESCENDENTS=""
help() {
    echo
    echo "-b, --bucket_id                          GCP bucket_id name to copy the files from"
		echo "-d, --date                               tree date to get sequences for (format: year-month-day)"
    echo "-i, --input                              input directory containing the following files:"
		echo
		echo "    allDescendants.txt                   file containing all descendants to get raw sequences for"
		echo "    gisaid_fullNames_date.fa.xz          fasta file containing all gisaid sequences for the given tree"
		echo "    cog_all.fasta.xz                     fasta file containing all cogUK sequences for the given tree"
		echo "    genbank.fa.xz                        fasta file containing all ncbi sequences for the given tree"
		echo
    echo "-o, --output                             output directory"
    echo
}

fastaSeqCount () {
    xcat $1 \
    | grep ^\> | wc -l
}
export -f fastaSeqCount


xcat () {
    local inFile=$1
    if [ "${inFile##*.}" == "xz" ]; then
        xzcat $inFile
    elif [ "${inFile##*.}" == "gz" ]; then
        zcat $inFile
    else
        cat $inFile
    fi
}
export -f xcat

fastaNames () {
    xcat $1 \
    | grep ^\> | sed -re 's/^>//;'
}
export -f fastaNames


# Command line flag processing
while getopts b:d:i:o: flag
do
        case "${flag}" in
                b) BUCKET_ID=${OPTARG} ;;
                d) TREE_DATE=${OPTARG} ;;
                i) INPUT=${OPTARG} ;;
                o) OUTPUT=${OPTARG} ;;
                *) echo "Invalid option: -$flag" ;;
        esac
done
DESCENDENTS="$INPUT/allDescendants.txt"

# Check correct number of command line arguments give
if [ "$#" -ne 8 ]; then
		echo "ERROR: Incorrect number of parameters given ($# arguments given)"
		# Display command line options
	  help
		exit 1
fi

if [ "$BUCKET_ID" == "None" ]
then
    echo "Running QC and filtration checks on-premise";
    # NCBI sequences
		cp genbank.fa.xz $INPUT
    # CogUK sequences
		cp cog_all.fasta.xz $INPUT
    # GISAID sequences, protobuf and metadata, if running private tree
    cp gisaidAndPublic.$TREE_DATE.masked.pb $INPUT
    cp gisaidAndPublic.$TREE_DATE.metadata.tsv.gz $INPUT
    cp metadata_batch_$TREE_DATE.tsv.gz $INPUT
		cp gisaid_fullNames_$TREE_DATE.fa.xz $INPUT
else
		echo "BUCKET_ID given, running on GCP instance";
    # Copy all sequence files needed into local directory
    # NCBI sequences
		gsutil cp gs://$BUCKET_ID/genbank.fa.xz $INPUT
    # CogUK sequences
		gsutil cp gs://$BUCKET_ID/cog_all.fasta.xz $INPUT
    # GISAID sequences, protobuf and metadata
    gsutil cp gs://$BUCKET_ID/gisaidAndPublic.$TREE_DATE.masked.pb $INPUT
    gsutil cp gs://$BUCKET_ID/gisaidAndPublic.$TREE_DATE.metadata.tsv.gz $INPUT
    gsutil cp gs://$BUCKET_ID/metadata_batch_$TREE_DATE.tsv.gz $INPUT
		gsutil cp gs://$BUCKET_ID/gisaid_fullNames_$TREE_DATE.fa.xz $INPUT
fi

# Create output log files
LOG_FILE="log_$TREE_DATE.log"

mkdir -p $OUTPUT
sort -u $DESCENDENTS > $OUTPUT/allDescendants.uniq.txt

#echo "Number of duplicates: " > $OUTPUT/$LOG_FILE
#echo `(wc -l $DESCENDENTS) - (wc -l $OUTPUT/allDescendants.uniq.txt)` >> $OUTPUT/$LOG_FILE

grep EPI_ISL $OUTPUT/allDescendants.uniq.txt | awk -F\| '{print $2 "\t" $0;}' > $OUTPUT/allDescendants.epiToName
#wc -l $OUTPUT/allDescendants.epiToName


grep -v EPI_ISL $OUTPUT/allDescendants.uniq.txt \
| egrep '[A-Z]{2}[0-9]{6}\.[0-9]+' \
| awk -F\| '{ if ($3 == "") { print $1 "\t" $0; } else { print $2 "\t" $0; } }' \
> $OUTPUT/allDescendants.gbAccToName
#wc -l $OUTPUT/allDescendants.gbAccToName


grep -v EPI_ISL $OUTPUT/allDescendants.uniq.txt \
| egrep -v '[A-Z]{2}[0-9]{6}\.[0-9]+' \
| egrep '^(England|Northern|Scotland|Wales)' \
| awk -F\| '{print $1 "\t" $0;}' \
> $OUTPUT/allDescendants.cogToName
#wc -l $OUTPUT/allDescendants.cogToName

chmod +x filtering/tawk
zcat $INPUT/metadata_batch_$TREE_DATE.tsv.gz \
| grep -Fwf <(cut -f 1 $OUTPUT/allDescendants.epiToName) \
| filtering/tawk '{print $3, $1 "|" $3 "|" $5;}' \
> $OUTPUT/allDescendants.epiToFastaName
#wc -l $OUTPUT/allDescendants.epiToFastaName

faSomeRecords <( xzcat $INPUT/gisaid_fullNames_$TREE_DATE.fa.xz ) \
<(cut -f 2 $OUTPUT/allDescendants.epiToFastaName) $OUTPUT/allDescendants.gisaid.fa

faSomeRecords <( xzcat $INPUT/genbank.fa.xz ) \
<(cut -f 1 $OUTPUT/allDescendants.gbAccToName) $OUTPUT/allDescendants.genbank.fa

faSomeRecords <( xzcat $INPUT/cog_all.fasta.xz ) \
<(cut -f 1 $OUTPUT/allDescendants.cogToName) $OUTPUT/allDescendants.cog.fa

join -t$'\t' <(sort $OUTPUT/allDescendants.epiToFastaName) <(sort $OUTPUT/allDescendants.epiToName) \
| cut -f 2,3 > $OUTPUT/allDescendants.rename

cat $OUTPUT/allDescendants.gbAccToName $OUTPUT/allDescendants.cogToName >> $OUTPUT/allDescendants.rename

wait
# Concatenate and rename fasta
cat $OUTPUT/allDescendants.gisaid.fa $OUTPUT/allDescendants.genbank.fa $OUTPUT/allDescendants.cog.fa \
| faRenameRecords stdin $OUTPUT/allDescendants.rename $OUTPUT/allDescendants.fa

fastaNames $OUTPUT/allDescendants.fa | sort > $OUTPUT/namesFound
comm -23 $OUTPUT/allDescendants.uniq.txt $OUTPUT/namesFound
