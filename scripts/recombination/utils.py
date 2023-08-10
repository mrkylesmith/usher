#!/bin/python3
# Utility functions for rivet-backend.py
import wget
import os
import subprocess
from termcolor import colored
import sys
import time
import os.path
import pathlib
import re
import json
import time
import timeit
import datetime
import yaml
import shutil
import gzip
from datetime import timedelta
from multiprocessing import Process
from google.cloud import storage

def get_config():
    if not os.path.isfile("ripples.yaml"):
        print(colored("[ERROR] FILE NOT FOUND. Make sure you copy the ripples.yaml config from the 'template/' directory into current directory."), "red")
        exit(1)
    with open('ripples.yaml') as f:
      data = yaml.load(f, Loader=yaml.FullLoader)
      return data 

def auth(project_id, key_file):
    cmd1 = ["gcloud", "config", "set", "project", project_id]
    cmd2 = ["gcloud", "auth", "activate-service-account", "--key-file", key_file]
    print("Setting up GCP access through gcloud.")
    subprocess.run(cmd1)
    subprocess.run(cmd2)

def uncompress_gz_file(file):
    """
    """
    if file.endswith(".gz"):
        with gzip.open(file, 'rb') as in_file:
            with open(os.path.splitext(file)[0], 'wb') as out_file:
                shutil.copyfileobj(in_file, out_file)

def non_public_tree_warning(msg):
    """
    """
    warning_msg_input = input(msg)
    if warning_msg_input.lower() == 'yes' or warning_msg_input.lower() == 'y':
        print(colored("Proceeding with RIVET backend pipeline using non-public tree.", "green"))
    elif warning_msg_input.lower() == 'no' or warning_msg_input.lower() == 'n':
        print("Exiting pipeline. Edit pipeline settings in config file.")
        exit(1)
    else:
        err_msg = "Type yes(y) or no(n):  "
        non_public_tree_warning(err_msg)

def check_config_format(config, LOCAL_PIPELINE):
    """
    """
    if not isinstance(config["num_descendants"], int):
      return (True, "num_descendants", int)
    # If running GCP job
    if not LOCAL_PIPELINE:
        # Check to make sure all necessary GCP config settings are provided
        parameters = [("bucket_id", str), ("project_id", str), ("key_file", str), 
                      ("instances", int), ("boot_disk_size", int),
                      ("machine_type", str)]
        for p in parameters:
            # If the field is empty or wrong type, exit and print error message
            if config[p[0]] is None or not isinstance(config[p[0]], p[1]):
                return (True, p[0], p[1])

def aggregate_results(bucket_id, results_dir, date):
    """
    Aggregate results locally from remote GCP node instances
    """
    current = str(os.getcwd())
    local_results = current + "/{}".format(results_dir)
    # Create local results output directory
    if os.path.isdir(local_results) is False:
        os.makedirs(local_results)
        os.makedirs("{}/results".format(local_results))

    # Create local temporary directory to copy remote results and aggregate 
    temp = local_results + "/merge_results/"
    os.makedirs(temp)

    # Make sure temp directory was created correctly
    if (os.path.isdir(temp) == False):
        print(colored("[ERROR]: Local results directory not created. Check naming error", 'red'))
        raise FileNotFoundError

    # Copy over all subdirectories from GCP Bucket to local temp directory
    remote_results = "gs://{}/{}/*".format(bucket_id, results_dir)
    subprocess.run(["gsutil", "cp", "-r", remote_results, temp])

    # File to aggregate all detected recombination events
    recombinants = open(local_results + "/results/filtered_recombinants_{}.txt".format(date), "w")
    unfiltered_recombinants = open(local_results + "/results/unfiltered_recombinants_{}.txt".format(date), "w")
    descendants = open(local_results + "/results/descendants_{}.txt".format(date), "w")

    # Aggregate the results from all remote machines and combine into one final file in 'results/' dir
    for directory in os.listdir(temp):
        subdir = temp + directory 
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
          if "descendants.tsv" in file:
            f3 = open(subdir + "/" +  file, "r")
            for line in f3:
              # One detected recombinant per line, aggregate all lines in each file
              descendants.write(line)
            f3.close()
    recombinants.close()
    descendants.close()
    unfiltered_recombinants.close()

def run_ripples_locally(version, mat, num_descendants):
    """
    """
    # Expecting ripples output (recombination.txt and descendents.txt)
    # in recombination/filtering to start this pipeline
    #command = [version, "-i", mat, "-n", "2", "-S", start, "-E", end, "-d", "filtering/data"]
    command = [version, "-i", mat, "-n", str(num_descendants), "-d", "filtering/data"]
    return command

def generate_translation(translation_outfile, mat, reference, NCBI_GENES):
    """
    """
    print(colored("Generating amino acid translations.", 'green'))
    subprocess.run(["matUtils", "summary", "--translate", translation_outfile, "-i", mat, "-g", NCBI_GENES, "-f", reference])

def build_taxonium_tree(mat, metadata, date, RESULTS, taxonium_config):
    """
    Build Taxonium tree jsonl file 
    """
    taxonium_file = str(date) + ".taxonium.jsonl.gz"

    # Check `usher_to_taxonium` executable reachable in path
    if(shutil.which('usher_to_taxonium') == None):
         print(colored("[ERROR]: (rivet) environment not activated.  Activate (rivet) Conda env", 'red'))
         exit(1)
    # Check if GENEBANK_FILE exists
    if (os.path.exists("hu1.gb") == False):
        print("Downloading GenBank file containing reference genome")
        wget.download("https://raw.githubusercontent.com/theosanderson/taxonium/master/taxoniumtools/test_data/hu1.gb")
        print()
        print(colored("GenBenk file download successful", 'green'))
    if os.path.exists(taxonium_config) is None:
        print(colored("{} file not found in current directory", "red"))
        exit(1)
    print(colored("Writing Taxonium tree file to: {}".format(taxonium_file), 'green'))
    build_taxonium_tree = [
                "usher_to_taxonium", "--input", mat, "--output", taxonium_file, 
                "--metadata", metadata, "--clade_type", 'nextstrain,pango',
                "--genbank", "hu1.gb", "--config_json", taxonium_config, 
                "--name_internal_nodes", "--columns", "genbank_accession,country,date,pango_lineage_usher"
                ]
    try:
        subprocess.run(build_taxonium_tree)
        shutil.move(taxonium_file, "{}/results/{}".format(RESULTS, taxonium_file))
        print(colored("[SUCCESS]: Taxonium file written to: {}".format(taxonium_file), 'green'))
    except subprocess.CalledProcessError as e:
        print(e.output)
        exit(1)


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

def create_bucket_folder(project_id, bucket_id, path_from_bucket_root, key_file):
    # Set environment variable with path to service account key file
    os.environ['GOOGLE_APPLICATION_CREDENTIALS'] = key_file
    # Check to make sure folder ends with backslash, add if not
    if not path_from_bucket_root.endswith("/"):
        path_from_bucket_root += "/"
    client = storage.Client(project=project_id)
    bucket = client.get_bucket(bucket_id)
    blob = bucket.blob(path_from_bucket_root)
    result = blob.upload_from_string('')

def convert(n):
    return str(datetime.timedelta(seconds = n))

def parse_command(version, mat, start, end, out, logging, bucket_id, results, reference, date, num_descendants):
    command = "python3 process.py {} {} {} {} {} {} {} {} {} {} {}".format(version,
            mat, start, end, out, bucket_id, results, reference, date, logging, num_descendants)
    return command

def run_chronumental(newick, metadata):
    """
    """
    reference_node = "DP0803|LC571037.1|2020-02-17"
    steps = "100"
    #TODO: Print message and write log file to results directory
    log_file = open("chronumental_stdout.log", "w")
    command = ["chronumental", "--tree", newick, "--dates", metadata, "--reference_node", reference_node, "--steps", steps]
    p = subprocess.Popen(command, stdout=log_file)
    return p

def newick_from_mat(mat, out_newick):
    """
    """
    command = ["matUtils", "extract", "-i", mat, "-t", out_newick]
    p = subprocess.run(command)
    return p

# Generate run command for final recombination list and ranking
def post_process(mat, filtration_results_file, chron_dates_file, date, recomb_output_file, metadata, output_dir):
    cmd = ["post_filtration", "-i", mat, "-f", filtration_results_file, "-c", chron_dates_file, "-d", date, "-r", recomb_output_file, "-w", "-m", metadata, "-o", output_dir]
    return cmd

def gcloud_run(command, log, bucket_id, machine_type, docker_image, boot_disk_size):
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

def gcloud_describe(operation_id, key_file):
     # Refresh service account credentials
     auth = ["gcloud", "auth", "activate-service-account", "--quiet", "--key-file", key_file]
     print("Refreshing Service Account credentials.")
     subprocess.run(auth)
     print("Checking job status")
     cmd = ["gcloud", "beta", "lifesciences",
            "operations", "describe", "--format=json", operation_id]
     out = subprocess.check_output(cmd)
     result = json.loads(out)
     done = 'done' in result and result['done'] == True
     return {'done': done, 'result': result}

