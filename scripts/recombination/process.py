#!/bin/python3
#
# Script to run ripples and filtration pipeline on remote GCP machine.

import time
import random
import sys
import subprocess
import datetime
import os
import os.path
import errno
import timeit
from datetime import timedelta

version = sys.argv[1]
mat = sys.argv[2]
start_range = sys.argv[3]
end_range = sys.argv[4]
# Subdir within "results" directory to place results
out = sys.argv[5]
bucket_id = sys.argv[6]
# Remote GCP Storage Bucket location to put end results
results = sys.argv[7]
reference = sys.argv[8]
date = sys.argv[9]
# Log to output runtime information for each partition
logging = sys.argv[10]  
num_descendants = sys.argv[11]
pipeline_dir = os.getcwd()

def convert(n):
    return str(datetime.timedelta(seconds = n))

def parse_ripples_command(version, mat, start, end, num_descendants):
    # Expecting ripples output (recombination.txt and descendents.txt)
    # in recombination/filtering to start this pipeline
    #command = [version, "-i", mat, "-n", "2", "-S", start, "-E", end, "-d", "filtering/data"]
    command = [version, "-i", mat, "-n", num_descendants, "-S", start, "-E", end, "-d", "filtering/data"]
    return command

# Start total runtime for instance
start = timeit.default_timer()

# Check starting directory is correct
if (os.path.exists("process.py") == False):
    print()
    print("ERROR: Process.py not found in current directory. Check starting directory in container.") 
    print()
    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), "process.py")

# Copy over protobuf, raw sequence file and reference from GCP Storage into irectory inside container
print("Copying input MAT: {} from GCP Storage bucket into local directory on remote machine.".format(mat))
if not os.path.exists(mat):
  subprocess.run(["gsutil", "cp", "gs://{}/{}".format(bucket_id, mat), pipeline_dir])
if not os.path.exists(reference):
  subprocess.run(["gsutil", "cp", "gs://{}/{}".format(bucket_id, reference), pipeline_dir])

# start runtime for RIPPLES
start_ripples = timeit.default_timer()

# Run ripples on current GCP instance
cmd = parse_ripples_command(version, mat, start_range, end_range, num_descendants)
with open("ripples_stdout", "w") as  ripples_stdout:
  with open("ripples_stderr", "w") as  ripples_stderr:
    subprocess.run(cmd,stderr=ripples_stderr,stdout=ripples_stdout)
# Stop timer for RIPPLES
stop_ripples = timeit.default_timer()
# Start runtime for filtration pipeline
start_filtration = timeit.default_timer()
filtration = ["./run_ripples_filtration.sh", mat, str(date), reference, results, out, bucket_id]
# Run putative recombinants through post-processing filtration pipeline
subprocess.run(filtration)

# Job finished, log runtime, copy to GCP logging directory
stop_filtration = timeit.default_timer()
stop = timeit.default_timer()

runtime_log = open("runtime_{}.log".format(out), "w")
runtime_log.write("Timing for recombination detection nodes:{}{}{}".format('\t', out, '\n'))
runtime_log.write("Time for ripples job:{}{}  (Hours:Minutes:Seconds){}".format('\t', str(timedelta(seconds=stop_ripples - start_ripples)), '\n'))
runtime_log.write("Time for filtration pipeline:{}{}  (Hours:Minutes:Seconds){}".format('\t', str(timedelta(seconds=stop_filtration - start_filtration)), '\n'))
runtime_log.write("Total runtime for nodes {}:{}{}  (Hours:Minutes:Seconds){}".format(out, '\t', str(timedelta(seconds=stop - start)), '\n'))
runtime_log.close()

subprocess.run(["gsutil", "cp", "runtime_{}.log".format(out), "gs://{}/{}".format(bucket_id, logging)])
subprocess.run(["gsutil", "cp", "check_mutant_log", "gs://{}/{}".format(bucket_id, logging)])
