from lib.general_helpers.configure_loggers import configure_demultiplexing_logger
from lib.demultiplexing.demultiplexing_pipelines.demultiplexing_isONclust_pipeline import isONclust_pipeline
from lib.general_helpers.configure_loggers import *
import os
import os.path as ospath

def run_demultiplexing(input_name:str, input_file_path: str, output_folder: str, demultiplexing_method: str ="isONclust", logger = None):
    
    if logger == None:
        logger = configure_demultiplexing_logger(output_folder, input_name)
        print(f"Logging set up at {output_folder}/{input_name}_demultiplexing_pipeline_log.log")
    
    if demultiplexing_method == "isONclust":
        logger.info(f"Running demultiplexing pipeline on {input_file_path} using the {demultiplexing_method} method")
        isONclust_pipeline(input_name,input_file_path,output_folder, logger=logger)
    
    