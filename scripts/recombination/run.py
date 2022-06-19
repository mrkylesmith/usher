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
import datetime
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
    per_instance = long_branches//instances
    partitions = [ (i ,min(i + per_instance - 1, long_branches - 1)) for i in range(0, long_branches, per_instance) ]
    return partitions

def convert(n):
    return str(datetime.timedelta(seconds = n))

def parse_command(mat, start, end, out):
    command = "python3 process.py {} {} {} {} {} {}".format(version,
            mat, start, end, out, bucket_id, results, reference, raw_sequences)
    return command

# Takes in .gz newick and metadata 
def parse_chron_command(newick, metadata, steps):
    comand = "chronumental --tree {} --dates {} --steps {}".format(newick, metadata, steps)
    return command

def gcloud_run(command):
    cmd = ["gcloud", "beta", "lifesciences", "pipelines", "run",
        "--location", "us-central1",
        "--regions", "us-central1",
        "--machine-type", machine_type, 
        "--boot-disk-size", boot_disk_size,
        "--logging", logging,
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
docker_image = "mrkylesmith/ripples-pipeline:latest"
boot_disk_size = config["boot_disk_size"]
instances = config["instances"] # Number of remote machines to parallelize ripples across 
machine_type = config["machine_type"]
logging = "gs://{}/logging/{}".format(bucket_id, config["logging"])
version = config["version"]
mat = config["mat"]
newick = config["newick"]
metadata = config["metadata"]
date = config["date"]
reference = config["reference"]
raw_sequences = config["raw_sequences"]

# Remote location in GCP Storge Bucket to copy filtered recombinants
#NOTE: Make sure this folder is created in Storage bucket ahead of time.
results = "gs://{}/{}".format(bucket_id, config["results"])

# Copy over protobuf from GCP storage bucket to local container
current = str(os.getcwd())
if not os.path.isfile("{}/{}".format(current,mat)):
  print("Copying input MAT: {} from GCP Storage bucket into local directory in container.".format(mat))
  subprocess.run(["gsutil", "cp", "gs://{}/{}".format(bucket_id, mat), current])
else:
    print("Input MAT found in local directory.")

# Run ripples init scripts 
init = "ripplesInit -i {}".format(mat)
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
print("long branches: {}, instances: {}, w: {}".format(long_branches, instances, branches_per_instance))
print("partitions: {}".format(partitions))

processes = []
for partition in partitions:

    start = str(partition[0])
    end = str(partition[1])
    out = "{}_{}".format(start, end)
   
    # python3 process.py <version> <tree.pb> <start> <end> <bucket_id> <output_dir> 
    command = parse_command(mat, start, end, out) 

    info = gcloud_run(command)
    processes.append({'partition': partition, 'operation_id': info['operation_id']})

completed = 0
while completed < instances:
   for process in processes:
     partition = process['partition']
     operation_id = process['operation_id']
     info = gcloud_describe(operation_id)
     done = info['done']
     if done:
       completed += 1
     print("partition: {}, operation_id: {}, done: {}".format(partition, operation_id, done))
   print("{} of {} completed".format(completed, instances))
   time.sleep(1)

print("All instance jobs have finished.  Aggregating results from remote machines.")

current = str(os.getcwd())

local_results = current + "/{}".format(config["results"])
# Create local output directory (this might have already been done???)
subprocess.run(["mkdir", "-p", local_results])

# Create local temporary directory to copy remote results and aggregate 
temp = current + "/merge_results/"
subprocess.run(["mkdir", "-p", temp])

# File to aggregate all detected recombination events
recombinants = open(local_results + "/recombinants_{}.txt".format(date), "w")

# Make sure temp directory was created correctly
if (os.path.isdir(temp) == False):
    print("Local results directory not created. Check naming error")
    raise FileNotFoundError

# Copy over all subdirectories from GCP Bucket to local temp directory
remote_results = "gs://{}/{}/*".format(bucket_id, config["results"])
subprocess.run(["gsutil", "cp", "-r", remote_results, temp])

# Aggregate the results from all remote machines and combine into one final file in 'results/' dir
for directory in os.listdir(temp):
    subdir = temp + directory 
    print("SUBDIR: ", subdir)
    # Guarantee to be only 2 files in each subdirectory
    files = os.listdir(subdir)
    for file in files:
      # Skip over descendents.tsv files for aggregating recombination events
      if "recombination" not in file:
          continue
      f = open(subdir + "/" +  file, "r")
      for line in f:
          # One detected recombinant per line, aggregate all lines in each file
          recombinants.write(line)
      f.close()

recombinants.close()

# Remove temp directory 
subprocess.run(["rm", "-r", temp])

print("Final recombination event results written to {}/recombinants_{}.txt".format(local_results,date))
