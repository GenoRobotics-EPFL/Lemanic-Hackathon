from lib.general_helpers.run_command import run_command
from lib.general_helpers.configure_loggers import *
from lib.demultiplexing.demultiplexing_helpers import *
import os
import os.path as ospath
from Bio import SeqIO
import pandas as pd

def isONclust_pipeline(input_file_path: str, output_folder: str, logger):
    """
    pipeline for clustering reads from a multiplexed fastq file using isONclust tool. Turns one fastq files into multiple ones, named arbitrarily using integers starting at 0.

    Parameters
    ----------
    input_name: name of test to include in logger
    input_file_path(str): path to fastq file to demultiplex
    logger(Logging.logger): logger 
    """
    
    if not ospath.exists(output_folder):
        os.makedirs(output_folder)
    tsv_path = ospath.join(output_folder, "final_clusters.tsv")

    logger.info(f"Running demultiplexing on {input_file_path}...")
    isONclust_command = f"isONclust --ont --fastq {input_file_path} --outfolder {output_folder}"
    logger.info(f"Running command: {isONclust_command}...")
    run_command(isONclust_command,logger)

    logger.info(f"Separating fastq into clusters for {input_file_path}...")
    clustering_command = f"isONclust write_fastq --clusters {tsv_path} --fastq {input_file_path} --outfolder {output_folder}  --N 1"
    logger.info(f"Running command: {clustering_command}...")
    run_command(clustering_command,logger)
    info_df = get_info_on_clustering(output_folder)
    # keep_best_clusters(info_df,output_folder)


    
    

    
