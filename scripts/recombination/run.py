#!/bin/python3
#
# Launch script to run parallel ripples jobs on GCP
import subprocess
import sys
import time
import os
import os.path
import pathlib
import re
import json
import time
import timeit
import datetime
from datetime import timedelta
import yaml


def get_config():
    with open('ripples.yaml') as f:
      data = yaml.load(f, Loader=yaml.FullLoader)
      return data 

def auth():
    cmd1 = ["gcloud", "config", "set", "project", project_id]
    cmd2 = ["gcloud", "auth", "activate-service-account", "--key-file", key_file]
    print("Setting up GCP access through gcloud.")
    subprocess.run(cmd1)
    subprocess.run(cmd2)

def get_partitions(long_branches, instances):
    partitions = []
    per_instance = long_branches // instances
    k = 0
    for i in range(1, instances+1):
        # Last partition gets extra 
        if i == instances:
            partitions.append((k, long_branches))
            break
        partitions.append((k, k + per_instance))
        k += per_instance + 1
    return partitions

def convert(n):
    return str(datetime.timedelta(seconds = n))

def parse_command(mat, start, end, out, logging):
    command = "python3 process.py {} {} {} {} {} {} {} {} {} {} {}".format(version,
            mat, start, end, out, bucket_id, results, reference, date, logging, num_descendants)
    return command

# Takes in .gz newick and metadata 
def parse_chron_command(newick, metadata, steps):
    command = "chronumental --tree {} --dates {} --steps {}".format(newick, metadata, steps)
    return command

def gcloud_run(command, log):
    cmd = ["gcloud", "beta", "lifesciences", "pipelines", "run",
        "--location", "us-central1",
        "--regions", "us-central1",
        "--machine-type", machine_type, 
        "--boot-disk-size", boot_disk_size,
        "--logging", "gs://{}/{}".format(bucket_id, log),
        "--docker-image", docker_image,
        "--command-line", command,
        #"--outputs", results,
        "--format=json"]

    print(" ".join(cmd))
    out = subprocess.check_output(cmd)
    result = json.loads(out)
    name = result['name']
    id = re.search("operations\/(\d+)", name).groups()[0]
    return {'operation_id': id, 'result': result}

def gcloud_describe(operation_id):
     cmd = ["gcloud", "beta", "lifesciences",
            "operations", "describe", "--format=json", operation_id]
     out = subprocess.check_output(cmd)
     result = json.loads(out)
     done = 'done' in result and result['done'] == True
     return {'done': done, 'result': result}

# Configs and credentials from yaml
config = get_config()

# Authenticate GCP account
bucket_id = config["bucket_id"]
project_id = config["project_id"]  
key_file = config["key_file"]

# Activate credentials for GCP Console and gcloud util
auth()

# Set ripples job config options
#docker_image = "mrkylesmith/ripples_pipeline:latest"
#TODO: Using dev image at the moment
docker_image = "mrkylesmith/ripples_pipeline_dev:latest"
boot_disk_size = str(config["boot_disk_size"])
instances = config["instances"] # Number of remote machines to parallelize ripples across 
machine_type = config["machine_type"]
logging = config["logging"]
version = config["version"]
mat = config["mat"]
newick = config["newick"]
metadata = config["metadata"]
date = config["date"]
reference = config["reference"]
num_descendants = config["num_descendants"]

# Remote location in GCP Storge Bucket to copy filtered recombinants
#TODO: Make sure this folder is created in Storage bucket ahead of time.
# Use gcloud client libary
results = "gs://{}/{}".format(bucket_id, config["results"])

# Copy over protobuf from GCP storage bucket to local container
current = str(os.getcwd())
if not os.path.isfile("{}/{}".format(current,mat)):
  print("Copying input MAT: {} from GCP Storage bucket into local directory in container.".format(mat))
  subprocess.run(["gsutil", "cp", "gs://{}/{}".format(bucket_id, mat), current])
else:
    print("Input MAT found in local directory.")

# Check logging file created on GCP bucket  #TODO: Check the logging file exists on GCP bucket, otherwise create it ...
if not logging.endswith("/"):
    logging += "/"

if num_descendants == None:
  init = "ripplesInit -i {}".format(mat)
elif isinstance(num_descendants, int):
  init = "ripplesInit -i {} -n {}".format(mat, num_descendants)
else:
    print("Check ripples.yaml file configuration for num_descendants. Provide int number of descendents to consider or leave field blank to use default value 2.")
    exit(1)

# Run ripples init scripts 
try:
  # Get number of long branches to search
  long_branches = int(subprocess.getoutput(init))
except:
  print("Empty tree input or error with input metadata. Check 'ripples.yaml' config file")
  exit(1)

print("Found {} long branches in {} to split across {} instances (number of instances set in 'ripples.yaml').".format(long_branches, mat, instances))

# Num of long branches to search per instance
branches_per_instance = long_branches//instances        

partitions = get_partitions(long_branches, instances)
print("long branches: {}, instances: {}, branches_per_instance: {}".format(long_branches, instances, branches_per_instance))
print("partitions: {}".format(partitions))


processes = []
completed = []
for partition in partitions:

    start_range = str(partition[0])
    end_range = str(partition[1])
    out = "{}_{}".format(start_range, end_range)
    log = logging + out + ".log"

    # The following command gets executed on remote machine: 
    # python3 process.py <version> <mat> <start> <end> <out> <bucket_id> <results> <reference> <date> <logging> <num_descendants>
    command = parse_command(mat, start_range, end_range, out, logging)

    info = gcloud_run(command, log)
    processes.append({'partition': partition, 'operation_id': info['operation_id']})
    completed.append(False)

# Start total runtime for all instances running in parallel 
start_time = timeit.default_timer()

while not all(completed):
   i = 0
   for process in processes:
     partition = process['partition']
     operation_id = process['operation_id']
     info = gcloud_describe(operation_id)
     done = info['done']
     if done:
       completed[i] = True
     print("partition: {}, operation_id: {}, done: {}".format(partition, operation_id, done))
     i+=1
   time.sleep(1)

print("All instance jobs have finished.  Aggregating results from remote machines.")

# All job have completed, log total runtime, copy to GCP logging directory
stop_time = timeit.default_timer()
runtime_log = open("aggregate_runtime.log", "w")
runtime_log.write("Timing for recombination detection tree date:{}{}{}".format('\t', date, '\n'))
runtime_log.write("Total runtime searching {} tree with {} long branches:{}{}  (Hours:Minutes:Seconds){}".format(date, long_branches, '\t', str(timedelta(seconds=stop_time - start_time)), '\n'))
runtime_log.close()
subprocess.run(["gsutil", "cp", "aggregate_runtime.log", "gs://{}/{}".format(bucket_id, logging)])

current = str(os.getcwd())

local_results = current + "/{}".format(config["results"])
# Create local output directory 
subprocess.run(["mkdir", "-p", local_results])

# Create local temporary directory to copy remote results and aggregate 
temp = current + "/merge_results/"
subprocess.run(["mkdir", "-p", temp])

# Make sure temp directory was created correctly
if (os.path.isdir(temp) == False):
    print("Local results directory not created. Check naming error")
    raise FileNotFoundError

# Copy over all subdirectories from GCP Bucket to local temp directory
remote_results = "gs://{}/{}/*".format(bucket_id, config["results"])
subprocess.run(["gsutil", "cp", "-r", remote_results, temp])

# File to aggregate all detected recombination events
recombinants = open(local_results + "/recombinants_{}.txt".format(date), "w")
unfiltered_recombinants = open(local_results + "/unfiltered_recombinants_{}.txt".format(date), "w")

# Aggregate the results from all remote machines and combine into one final file in 'results/' dir
# TODO: gcloud client function
for directory in os.listdir(temp):
    subdir = temp + directory 
    print("SUBDIR: ", subdir)
    # Guarantee to be only 3 files in each subdirectory
    files = os.listdir(subdir)
    for file in files:
      # Skip over recombination.tsv and descendents.tsv files for aggregating recombination events
      if "filtered_recombinants.txt" in file:
        f1 = open(subdir + "/" +  file, "r")
        for line in f1:
          # One detected recombinant per line, aggregate all lines in each file
          recombinants.write(line)
        f1.close()
      if "recombination.tsv" in file:
        f2 = open(subdir + "/" +  file, "r")
        for line in f2:
          # One detected recombinant per line, aggregate all lines in each file
          unfiltered_recombinants.write(line)
        f2.close()
      #TODO: Aggregate all descendents
recombinants.close()
unfiltered_recombinants.close()

# Remove temp directory 
#subprocess.run(["rm", "-r", temp])

print("Final recombination event results written to {}/recombinants_{}.txt".format(local_results,date))
