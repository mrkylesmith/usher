# Helper functions for sam2msa program.
#
#
import pysam

def generate_reference_contig(reference_file):
  """
  Function to add SARS-CoV2 reference as first sequence in MSA file.

  Args:
    reference_file: SARS-CoV2 reference

  Returns:
    SARS-CoV2 reference in MSA format to add to MSA file.
  """
  reference = open(reference_file,"r")
  lines = reference.readlines()
  reference.close()
  ref = ""
  count=0
  for line in lines:
      read = line.split(" ")
      if(count==0):
          refName = read[0].split(">")[1].split("\n")[0]
          count+=1
          continue
      else:
          ref += read[0].split("\n")[0]
  # Return reference as a single long string
  return ref

def get_cigar(read):
  """
  Get the operator and length for CIGAR string of read.

  Args:
    read: Given query sequence

  Returns:
    Operator and length given from CIGAR string
  """
  return read.cigartuples[0][0], read.cigartuples[0][1]

def get_ref_name(samfile):
  """
  Get the name of the reference sequence.

  Args:
    samfile: Path to samfile to extract reference name from.

  Returns:
    Return name of the reference sequence used for alignment.
  """
  for read in samfile:
      return read.reference_name

def align_seq_template(reference_start):
  """
  Generate alignment template with "-" before query read aligns.

  Args:
    reference_start: Position in reference where query read starts aligning.

  Returns:
    Alignment template up to where query starts aligning to reference.
  """
  aligned_seq = ""
  aligned_seq += "-" * reference_start
  return aligned_seq


def cigar_table(operator, length, query, alignment_template, query_index):
  """
  Modify alignment template for query read given CIGAR string.

  Args:
    operator: Type of operation in CIGAR string: "M", "I", "D", "H", "S", "N"
    length: Length to apply given operation for CIGAR string
    alignment_template: Alignment of query to reference to update
    query_index: Position in query that has been aligned.

  Returns:
    Updated alignment template after given CIGAR operation/length is applied
  """
    # "M" Match; nucleotide is present in the reference.
    if operator == 0:
       # Add matching query sequence directly to alignment template
       alignment_template += query[query_index:length + query_index]
       # Update query index
       query_index += length
       # Return the new alignment_template, specific to this read
       return alignment_template, query_index

    # "I" Insertion; the nucleotide is present in the read but not in the reference.
    elif operator == 1:
       query_index += length
       return alignment_template, query_index

    # "D" Deletion; the nucleotide is present in the reference but not in the read.
    elif operator == 2:
       # Add gaps to query sequence where nucleotide/s are missing in read
       alignment_template += "-" * length
       return alignment_template, query_index

    # "N" Skipped region; a region of nucleotides is not present in the read.
    elif operator == 3:
       alignment_template += "-" * length
       # Update query index
       query_index += length
       return alignment_template, query_index

    # "S" Soft Clipping;  the clipped nucleotides are present in the read.
    elif operator == 4:
       alignment_template += query[query_index:length + query_index]
       query_index += length
       return alignment_template, query_index

    # "H" Hard Clipping; the clipped nucleotides are not present in the read.
    elif operator == 5:
       # Move query index b/c nucleotides not present in read, alignment template stays the same
       query_index += length
       return alignment_template, query_index

    #TODO: Add more CIGAR string options if necessary.
    return None
