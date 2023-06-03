#!/bin/python3
#
# Launch script to run parallel ripples jobs on GCP
import subprocess
from multiprocessing import Process
from google.cloud import storage
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
from termcolor import colored
import yaml
import shutil
import wget
import utils


def main():
  # Configs and credentials from yaml
  config = utils.get_config()
  bucket_id = config["bucket_id"]
  project_id = config["project_id"]  
  # Authentication keys for GCP account
  key_file = config["key_file"]
  # Required parameters
  version = config["version"]
  mat = config["mat"]
  newick = config["newick"]
  metadata = config["metadata"]
  date = config["date"]
  reference = config["reference"]
  num_descendants = config["num_descendants"]
  LOCAL_PIPELINE = False
  PUBLIC_TREE = config["public_tree"]
  PATH = str(os.getcwd())
  RESULTS = PATH + "/" + config["results"]
  LOGGING = PATH + "/{}/logging".format(config["results"])

  # Check that executables found
  if shutil.which("ripplesInit") is None or shutil.which("ripplesUtils") is None:
      print(colored("[ERROR] ENV not setup correctly. Activate rivet Conda environment in 'install/rivet_env.yml'", "red"))
      exit(1)

  # Check `usher_to_taxonium` executable
  if shutil.which("usher_to_taxonium") is None:
       print(colored("[ERROR]: rivet environment not activated.  Activate (rivet) Conda env in 'install/rivet_env.yml'", 'red'))
       exit(1)

  if config["generate_taxonium"] is True:
      taxonium_file = "taxonium_config.json"
      if not os.path.isfile("{}/{}".format(PATH, taxonium_file)):
          print(colored("[ERROR] Taxonium config file '{}' not found in current directory".format(taxonium_file), 'red'))
          exit(1)

  if config["public_tree"] is not True:
      print(colored("[WARNING] Pipeline use recommended for running public UShER trees only. Errors will occur without authorized access to sequencing data.", "yellow"))
      msg = "Do you wish to continue pipeline with non-public UShER tree? (yes/no):  "
      utils.non_public_tree_warning(msg)

  if not bucket_id and not project_id:
      print(colored("Running RIVET backend pipeline locally.", 'green'))
      LOCAL_PIPELINE = True

      config_err = utils.check_config_format(config, LOCAL_PIPELINE)
      if config_err is not None and config_err[0]:
          print("[ERROR] Check ripples.yaml file configuration field {}. Must be type: ".format(config_err[1]), config_err[2])
          exit(1)

      if config["results"] is not None:
          print("Creating output {} folder to save all backend pipeline outputs".format(config["results"]))
          # Create results directory on local machine
          os.makedirs(RESULTS, exist_ok=True)
          # Create logging directory on local machine
          os.makedirs(LOGGING, exist_ok=True)
      else:
          print(colored("[ERROR] Fill in 'results' and 'logging' fields in YAML file where output files will be written", 'red'))
          exit(1)

      # Check that necessary sequence files are present in current directory
      if not os.path.isfile("{}/{}".format(PATH, mat)):
          print(colored("[ERROR] {} not found in current directory: {}".format(mat, PATH), "red"))
          exit(1)
      if not os.path.isfile("{}/{}".format(PATH, metadata)):
          print(colored("[ERROR] {} not found in current directory: {}".format(metadata, PATH), "red"))
          exit(1)

      if not os.path.exists(reference):
          print("Downloading SARS-CoV-2 reference genome.")
          wget.download("https://storage.googleapis.com/public_trees/reference.fa")
          os.rename(reference, "reference.fa")
          print()

      # Check if public sequence files are present in current directory, retrieve them if not
      if not os.path.isfile("{}/{}".format(PATH, "genbank.fa.xz")):
          try:
              print("Downloading genbank.fa.xz sequence file")
              wget.download("https://hgwdev.gi.ucsc.edu/~angie/sarscov2phylo/ncbi.{}/genbank.fa.xz".format(date))
          except:
              print(colored("[ERROR] genbank.fa.xz public seqeunce file not downloaded. Check link", "red"))
              exit(1)
      if not os.path.isfile("{}/{}".format(PATH, "cog_all.fasta.xz")):
          try:
              print("Downloading cog_all.fasta.xz sequence file")
              wget.download("https://hgwdev.gi.ucsc.edu/~angie/sarscov2phylo/cogUk.{}/cog_all.fasta.xz".format(date))
          except:
              print(colored("[ERROR] cog_all.fasta.xz public seqeunce file not downloaded. Check link", "red"))
              exit(1)

      if not PUBLIC_TREE:
          if not os.path.isfile("{}/{}".format(PATH, "metadata_batch_{}.tsv.gz".format(date))):
              print(colored("[ERROR] {} not found in current directory: {}".format("metadata_batch_{}.tsv.gz".format(date), PATH), "red"))
              exit(1)
          if not os.path.isfile("{}/{}".format(PATH, "gisaid_fullNames_{}.fa.xz".format(date))):
              print(colored("[ERROR] {} not found in current directory: {}".format("gisaid_fullNames_{}.fa.xz".format(date), PATH), "red"))
              exit(1)

      # Install 3SEQ locally into filtering directory
      if not os.path.isfile("{}/{}".format(PATH, "/filtering/3seq/3seq")) or not os.path.isfile("{}/{}".format(PATH, "/filtering/3seq/my3seqTable700")):
          os.chdir(PATH + "/filtering")
          try:
              subprocess.run(["./3seq_install.sh"])
          except:
              print(colored("[ERROR] 3SEQ not installed correctly.", "red"))
              os.chdir(PATH)
              exit(1)
          # Change back to parent recombination directory
          os.chdir(PATH)

  else:
      print(colored("Running RIVET backend pipeline on GCP, configuring setup.", 'green'))

      config_err = utils.check_config_format(config, LOCAL_PIPELINE)
      if config_err is not None and config_err[0]:
          print(colored("[ERROR] Check ripples.yaml file configuration field: '{}'. Must be type: {}".format(config_err[1], config_err[2]), "red"))
          exit(1)

      # Check `gsutil` utility reachable in path
      if shutil.which("gsutil") is None:
          print(colored("[ERROR]: rivet environment not activated.  Activate (rivet) Conda env in 'install/rivet_env.yml'", 'red'))
          exit(1)

      # Activate credentials for GCP Console and gcloud util
      utils.auth(project_id, key_file)

      # Set ripples job config options
      docker_image = str(config["docker_image"])
      boot_disk_size = str(config["boot_disk_size"])
      # Number of remote machines to parallelize ripples across 
      instances = config["instances"]
      machine_type = config["machine_type"]
      
      # Remote location in GCP Storge Bucket to copy filtered recombinants
      results = "gs://{}/{}".format(bucket_id, config["results"])

      # Check logging file created on GCP bucket
      #if not logging.endswith("/"):
      #    logging += "/"

      # Create remote logging folder
      utils.create_bucket_folder(project_id, bucket_id, "logging", key_file)
      print("Created empty GCP storage bucket folder for logging: logging/")

      # Create remote results folder
      utils.create_bucket_folder(project_id, bucket_id, config["results"], key_file)
      print("Created empty GCP storage bucket folder for results: {}".format(config["results"]))

      # Copy over protobuf and metadata file from GCP storage bucket to local container
      if not os.path.isfile("{}/{}".format(PATH, mat)):
          print("Copying input MAT: {} from GCP Storage bucket into local directory in container.".format(mat))
          subprocess.run(["gsutil", "cp", "gs://{}/{}".format(bucket_id, mat), PATH])
      else:
          print("Input MAT found in local directory.")

      if not os.path.isfile("{}/{}".format(PATH, metadata)):
          print("Copying input MAT metadata: {} from GCP Storage bucket into local directory in container.".format(metadata))
          subprocess.run(["gsutil", "cp", "gs://{}/{}".format(bucket_id, metadata), PATH])
      else:
          print("Input MAT metadata found in local directory.")

  # Generate newick tree input for Chronumental, if doesn't exist already
  print("Generating newick file from MAT using matUtils extract")
  utils.newick_from_mat(mat, newick)
  print(colored("Newick tree generated: {}".format(newick), "green"))

  # Launch Chronumental job locally
  p1 = Process(target=utils.run_chronumental, args=(newick, metadata))
  p1.start()
  print(colored("Chronumental job finished: chronumental_dates_{}.tsv".format(metadata), "green"))

  print("Finding the number of long branches.")
  # Using default num_descendants
  if num_descendants == None:
      init = "ripplesInit -i {}".format(mat)
  elif isinstance(num_descendants, int):
      init = "ripplesInit -i {} -n {}".format(mat, num_descendants)
  else:
      print("[ERROR] Check ripples.yaml file configuration for num_descendants. Provide int number of descendents to consider or leave field blank to use default value 2.")
      exit(1)

  # Run ripples init scripts 
  try:
      # Get number of long branches to search
      long_branches = int(subprocess.getoutput(init))
  except:
      print(colored("[ERROR] Empty tree input or error with input metadata. Check 'ripples.yaml' config file", "red"))
      # Exit early, kill Chronumental process
      p1.terminate()
      exit(1)
  print(colored("Found {} long branches in {}.".format(long_branches, mat), "green"))

  # Running pipeline on GCP
  if not LOCAL_PIPELINE:
      # Num of long branches to search per instance
      branches_per_instance = long_branches//instances        
      partitions = utils.get_partitions(long_branches, instances)
      print("long branches: {}, instances: {}, branches_per_instance: {}".format(long_branches, instances, branches_per_instance))
      print("partitions: {}".format(partitions))

      processes = []
      completed = []
      for partition in partitions:
          start_range = str(partition[0])
          end_range = str(partition[1])
          out = "{}_{}".format(start_range, end_range)
          log = "logging/" + out + ".log"

          # The following command gets executed on remote machine: 
          # python3 process.py <version> <mat> <start> <end> <out> <bucket_id> <results> <reference> <date> <logging> <num_descendants>
          command = utils.parse_command(version, mat, start_range, end_range, out, "logging", bucket_id, results, reference, date, num_descendants)

          info = utils.gcloud_run(command, log, bucket_id, machine_type, docker_image, boot_disk_size)
          processes.append({'partition': partition, 'operation_id': info['operation_id']})
          completed.append(False)
      # Start total runtime for all instances running in parallel 
      start_time = timeit.default_timer()
      while not all(completed):
          i = 0
          for process in processes:
              partition = process['partition']
              operation_id = process['operation_id']
              info = utils.gcloud_describe(operation_id, key_file)
              done = info['done']
              if done:
                  completed[i] = True
              print("partition: {}, operation_id: {}, done: {}".format(partition, operation_id, done))
              i+=1
          # Query active GCP jobs every 2 mins
          time.sleep(120)
      print(colored("All instance jobs have finished.  Aggregating results from remote machines.", "green"))

      # All job have completed, log total runtime, copy to GCP logging directory
      stop_time = timeit.default_timer()
      runtime_log = open("aggregate_runtime.log", "w")
      runtime_log.write("Timing for recombination detection tree date:{}{}{}".format('\t', date, '\n'))
      runtime_log.write("Total runtime searching {} tree with {} long branches:{}{}  (Hours:Minutes:Seconds){}".format(date, long_branches, '\t', str(timedelta(seconds=stop_time - start_time)), '\n'))
      runtime_log.close()
      subprocess.run(["gsutil", "cp", "aggregate_runtime.log", "gs://{}/logging".format(bucket_id)])
      # Aggregate results from each GCP instance locally
      utils.aggregate_results(bucket_id, config["results"], date)
  
  if LOCAL_PIPELINE:
      print(colored("Beginning recombination inference on premise", "green", attrs=['reverse']))
      # Start total runtime for local job
      start = timeit.default_timer()

      # No partitioning of long branches for running on premise
      start_range = 0 
      end_range = long_branches
      out = "{}_{}".format(start_range, end_range)
      log = LOGGING + "/" + out + ".log"

      # start runtime for RIPPLES
      start_ripples = timeit.default_timer()

      # Run ripples job locally
      cmd = utils.run_ripples_locally(version, mat, num_descendants)
      #TODO: add verbose option to config file for stdout 
      #with open("{}/ripples_stdout".format(LOGGING), "w") as  ripples_stdout:
      #    with open("{}/ripples_stderr".format(LOGGING), "w") as  ripples_stderr:
      #        subprocess.run(cmd, stderr=ripples_stderr, stdout=ripples_stdout)
      with open("{}/ripples_stderr".format(LOGGING), "w") as  ripples_stderr:
          subprocess.run(cmd, stderr=ripples_stderr)

      # Stop timer for RIPPLES
      stop_ripples = timeit.default_timer()
      print(colored("Recombination inference finished in {}{}  (Hours:Minutes:Seconds){}".format('\t', str(timedelta(seconds=stop_ripples - start_ripples)), '\n') , "green", attrs=['reverse']))

      print(colored("Starting QC and filtration pipeline.", "green"))
      # Start runtime for filtration pipeline
      start_filtration = timeit.default_timer()

      # Run putative recombinants through post-processing filtration pipeline
      filtration = ["./run_ripples_filtration.sh", mat, str(date), reference, RESULTS, out, str(bucket_id), LOGGING]
      subprocess.run(filtration)

      # Job finished, log runtime
      stop_filtration = timeit.default_timer()
      stop = timeit.default_timer()
      runtime_log = open("{}/aggregate_runtime.log".format(LOGGING), "w")
      runtime_log.write("Timing for recombination detection tree date:{}{}{}".format('\t', date, '\n'))
      runtime_log.write("Total runtime searching {} tree with {} long branches:{}{}  (Hours:Minutes:Seconds){}".format(date, long_branches, '\t', str(timedelta(seconds=stop - start)), '\n'))
      runtime_log.close()


      # Final filtered results file
      if not os.path.isdir("{}/results".format(RESULTS)):
          os.makedirs("{}/results".format(RESULTS), exist_ok=True)
      if os.path.isfile("{}/filtering/data/filtered_recombinants.txt".format(PATH)):
          os.rename("{}/filtering/data/filtered_recombinants.txt".format(PATH),"{}/results/filtered_recombinants_{}.txt".format(RESULTS, date))
      print(colored("QC filtration pipeline complete", "green"))

  # Move logging directory into results directory
  #shutil.move(LOGGING, RESULTS)

  # Check to make sure Chronumental job finished successfully
  while(p1.is_alive()):
      print("Chronumental job not finished running yet.")
      time.sleep(20)

  # Check if metadata file is unzipped, if not unzip it for post_filtration pipeline
  utils.uncompress_gz_file(metadata)

  # Assumes Chronumental inferred dates file output to current directory
  filtration_results_file = "{}/results/filtered_recombinants_{}.txt".format(RESULTS, date)
  chron_dates_file = "chronumental_dates_{}.tsv".format(metadata)
  recomb_output_file = "{}/results/final_recombinants_{}.txt".format(RESULTS, date)

  print(colored("Ranking recombinant results and aggregating additional information.", "green"))
  # Rank recombinants and generate final recombinant node information output file
  cmd = utils.post_process(mat, filtration_results_file, chron_dates_file, str(date), recomb_output_file, metadata, "{}/results".format(RESULTS))
  try:
      subprocess.run(cmd)
  except:
      print("[Error] Ranking unsuccessful, check formatting, rerun ranking with 'post_filtration --help'", "red")
      exit(1)
  node_to_extract_file = "{}/results/all_trio_nodes.txt".format(RESULTS)
  #cmd="sort -u <(cut -f 1 {0} |tail -n +2 ) <(cut -f 2 {0} |tail -n +2 ) <(cut -f 3 {0} |tail -n +2 ) > {1}".format(recomb_output_file,node_to_extract_file)
  #subprocess.run(["bash","-c",cmd])
  # Create VCF file of all trio sequences
  subprocess.run(["matUtils","extract","-i",mat,"-s",node_to_extract_file,"-v","{}/results/trios.vcf".format(RESULTS)])

  # If running public tree, build Taxonium tree jsonl file 
  if PUBLIC_TREE and config["generate_taxonium"] is True:
      utils.build_taxonium_tree(mat, metadata, date, RESULTS)

  # Copy over final results file to GCP storage bucket
  if not LOCAL_PIPELINE:
      subprocess.run(["gsutil", "cp", recomb_output_file, results])
      subprocess.run(["gsutil", "cp", "{}/trios.vcf".format(RESULTS), results])
      subprocess.run(["gsutil", "cp", "{}/samples_descendants.txt.xz".format(RESULTS), results])

  print(colored("[Success] RIVET Pipeline complete. Final recombination event results written to {}".format(recomb_output_file), "green"))

if __name__ == '__main__':
    main()
