#!/bin/usr/python3
#
# Script to generate VCF files for given of amplicon FASTQ read files.
'''
get_variants: Script to generate VCF file for given set of amplicon FASTQ read files. 
'''

import argparse
import errno
import os
from os import makedirs
from os.path import abspath, expanduser, isdir, isfile
import gzip
import sys
import subprocess
from util import *
from sam_to_msa import sam_to_msa

# Downloaded reference files
mask_sites = "problematic_sites_sarsCov2.vcf"
covid_ref = "wuhCor1.fa"
mask_sites_link = "https://raw.githubusercontent.com/W-L/ProblematicSites_SARS-CoV2/master/problematic_sites_sarsCov2.vcf"
covid_ref_link = "http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/public-latest.all.masked.pb.gz"

# Check if reference files exist, otherwise download 
if (os.path.exists(mask_sites) == False):
    print()
    print("Downloading {}".format(mask_sites))
    subprocess.run(["wget", mask_sites_link])
    print()
else:
    print("problematic_sites_sarsCov2.vcf found in current directory.")

# Check if reference files exist, otherwise download 
if (os.path.exists(covid_ref) == False):
    print()
    print("Downloading {}".format(covid_ref))
    subprocess.run(["wget", covid_ref_link])
    print()
else:
    print("SARS-Cov2 reference file (wuhCor1.fa) found in current directory.")

# Check correct masking sites file downloaded
if not isfile(mask_sites):
    print()
    print("ERROR: 'problematic_sites_sarsCov2.vcf' file not found. Make sure to download this file before running this pipeline.")
    print("This file can be downloaded here:")
    print()
    print("wget https://raw.githubusercontent.com/W-L/ProblematicSites_SARS-CoV2/master/problematic_sites_sarsCov2.vcf")
    print()
    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), mask_sites)

# Check correct covid reference genome downloaded
if not isfile(covid_ref):
    print()
    print("ERROR: 'wuhCor1.fa' SARS-CoV2 reference file not found. Make sure to download this file before running this pipeline.")
    print("This file can be downloaded and uncompressed here:")
    print()
    print("wget https://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/bigZips/wuhCor1.fa.gz")
    print("gunzip wuhCor1.fa.gz")
    print()
    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), covid_ref)

# Option parsing
parser = argparse.ArgumentParser(description=__doc__)
req_grp = parser.add_argument_group(title='Required')
parser.add_argument("-a", "--align",  default="minimap2", type=str, help="Specify aligner to use. default = minimap2")
parser.add_argument("-i", "--reads", required=True, type=str, help="Give input directory containing FASTQ files with amplicon reads to place.")
parser.add_argument("-o", "--output", required=True, type=str, help="Give output directory to put intermediate data and output results.")
args = parser.parse_args()

# Command line args
input_reads = args.reads                                           # Directory containing input FASTQ amplicon reads
aligner = args.align
output = abspath(expanduser(args.output))
# Create /data directory inside given output directory
subdir = abspath(output + "/data")

# Intermediate output files
amplicon_reads = "{}/amplicons.dedup.fastq".format(subdir)         # Aggreated/deduplicated FASTQ file
samfile = "{}/amplicons.sam".format(subdir)                        # Output alignment file
amplicons_msa = "{}/amplicons.fa".format(subdir)                   # Output MSA of all amplicon reads
output_vcf = "{}/amplicons.vcf".format(subdir)                     # VCF generated from faToVcf
fastq_name_mapping = "{}/fastq_name_mapping.txt".format(subdir)    # Mapping FASTQ -> placement naming

# Check if directory containing amplicon reads exists.
if not isdir(input_reads):
    print()
    print("ERROR: Amplicon reads directory not found.  Make sure all the amplicon read files are located in {} directory.".format(input_reads))
    print()
    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), input_reads)

#TODO: Add viralMSA as alternate aligner.
# Check if valid aligner given
if (aligner != "minimap2"):
    # Print options function
    print()
    print("ERROR: Invalid aligner given. See options: `python3 get_variants.py --help`  ")
    print()
    raise argparse.ArgumentTypeError("Invalid aligner specified, please see options for specifying aligner with -a flag above.")

# Check that given output directory doesn't exist already
if isdir(output) or isfile(output):
    print("Output directory given: {}  already exists.".format(output))
    exit(1)

print("Reading FASTQ files from {} directory: ".format(input_reads))
print("Using {} aligner".format(aligner))

# Create output directory for results and data directory to dump intermediate file outputs
makedirs(output)
makedirs(subdir)

# Make sure filtering script has permissions 
subprocess.run(["chmod", "+x", "filter.sh"])
    
# Call filter script: ./filter.sh <input_reads> <output>
run_filter = [ "./filter.sh", input_reads, output ]
p = subprocess.run(run_filter)
if p.returncode != 0:
    print("ERROR MESSAGE: ", p)
    print("EXIT CODE: ", p.returncode)
    raise Exception("Exit code: {} from filter.sh script".format(p.returncode))
print("FASTQ reads aggregation and deduplication finished")


# Make sure alignment script has permissions 
subprocess.run(["chmod", "+x", "align.sh"])

# Call alignment script: ./aligner <aligner> <amplicon_reads> <samfile>
run_aligner = [ "./align.sh", aligner, amplicon_reads, samfile ]
p = subprocess.run(run_aligner)
if p.returncode != 0:
    print("ERROR MESSAGE: ", p)
    print("EXIT CODE: ", p.returncode)
    raise Exception("Exit code: {} from align.sh script".format(p.returncode))
print("Alignment with {} finished and alignment file written to: {}".format(aligner, samfile))

# Check that SAM file was created
if (os.path.exists(samfile) == False):
    print("ERROR: Alignment file, {} not found.".format(samfile))
    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), samfile)

# Generate MSA file from SAM file
sam_to_msa(samfile, amplicons_msa, fastq_name_mapping)

# Check MSA file was generated
if (os.path.exists(amplicons_msa) == False):
    print("ERROR: MSA file, {} not found.".format(amplicons_msa))
    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), amplicons_msa)

print("MSA file generated.  Starting variant calling.")

# Use amplicon MSA file to call variants 
generate_msa = ["faToVcf", "-maskSites=problematic_sites_sarsCov2.vcf", amplicons_msa, output_vcf]
subprocess.run(generate_msa)
p = subprocess.run(generate_msa)
if p.returncode != 0:
    print("ERROR MESSAGE: ", p)
    print("EXIT CODE: ", p.returncode)
    raise Exception("Exit code: {} from faToVcf".format(p.returncode))
print("Variant calling with faToVcf finished and VCF file written to: {}".format(output_vcf))

# Check VCF file exists
if (os.path.exists(output_vcf) == False):
    print("ERROR: Variant call file, {} not found.".format(output_vcf))
    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), output_vcf)

print("SUCCESS: Variant call file generated for input amplicon reads: {}".format(output_vcf))
print("Output SAM file written here: {}".format(samfile))
print("Output MSA file written here: {}".format(amplicons_msa))
print("Output VCF file written here: {}".format(output_vcf))
