# Convert SAM file to MSA format.
#
#
import pysam
import sys
from util import *

def sam_to_msa(samfile, msafile, fastq_name_mapping):
  """
  Function to convert a SAM file containing alignment information to an MSA file.
  Imports from util.py, containing various helper functions to generate MSA file.

  Args:
    samfile: A path to the SAM file.
    msafile: A path to write MSA file to.

  Returns:
    Void function. Writes to given msafile file path.
  """
  # Generate config string for reference
  reference_file = "wuhCor1.fa"
  ref = generate_reference_contig(reference_file)
  
  # Load SAM and MSA files 
  samfile = pysam.AlignmentFile(samfile, "r")
  msafile = open(msafile,"w")
  ref_name = "NC_045512v2"

  # Write SARS-CoV2 reference first to MSA file
  msafile.write(">" + ref_name + "\n" + ref + "\n")

  # Amplicon reads will be named uniquely "amplicon_$i", 
  query_name = "amplicon_"
  # File containing mapping from original amplicon FASTQ naming 
  # to new naming given for placement
  fastq_name_mapping_file = open(fastq_name_mapping, "w")

  i = 1
  for read in samfile:
    # If no CIGAR string present in SAM file, don't consider read
    if read.cigartuples == None:
        print("CIGAR string not present. Discarding read with id: {}".format(read.query_name))
        continue

    # If "H" hardclipping present in CIGAR string, don't consider read
    if "H" in read.cigarstring:
        print("CIGAR string contains hard clipping ('H'). Discarding read with id: {}".format(read.query_name))
        continue

    # Discard query sequences with any "Ns" in them
    if "N" in read.cigarstring:
        print("CIGAR string contains skipped region ('N'). Discarding read with id: {}".format(read.query_name))
        continue

    # Generate new alignment template for read
    alignment_template = align_seq_template(read.reference_start)

    # "H" and "N" CIGAR strings will be filtered to this point already, handle just "S"
    if read.cigartuples[0][0] == 4:
      alignment_template = align_seq_template(read.reference_start - read.cigartuples[0][1])
   
    # Get the query sequence  and give read unique name
    query = read.query_sequence
    read_name = query_name + str(i)
    i+=1

    query_index = 0
    # If query happens to align off the end of reference sequence, skip it
    if read.reference_end > len(ref):
        print("Error. Query aligns off end of reference. Skipping this read.")
        continue

    # Check if first code in CIGAR string is not a match ("M")
    #NOTE: If beginning of the query is not a match,
    # reset reference coordinate to part of the query sequence that is aligning to reference
    first = True
    for code in read.cigartuples:
      operator = code[0]
      length = code[1]

      # Get updated alignment template for each code in CIGAR string
      alignment_template, query_index = cigar_table(operator, length, query, alignment_template, query_index)
      first = False
   
      # Read contains CIGAR string that is not handled by cigar_table function.
      if alignment_template == None:
          print("READ contains CIGAR string code that has not been handled.")
          exit(1)
    
    # Write updated alignment template to MSA file for the read
    alignment_template += "-" * (len(ref) - len(alignment_template))
    if len(alignment_template) != len(ref):
        continue
    msafile.write(">" + read_name + '\t' + str(query_index) + '\t' + str(read.reference_start + 1) + '\t' + str(read.reference_end) + "\n")
    msafile.write(alignment_template + "\n")

    # Map original FASTQ amplicon name to new assigned amplicon name
    fastq_name_mapping_file.write(read.query_name + "\t" + read_name + "\n")

  msafile.close()
  samfile.close()
  fastq_name_mapping_file.close()
