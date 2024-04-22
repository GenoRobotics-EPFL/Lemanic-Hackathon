import os
import os.path as ospath
import pandas as pd
from standard_pipeline import pipeline_single
import numpy as np
from lib.identification.identification_pipelines.identification_processing import get_best_species_from_xml
import argparse
from Bio import SeqIO
import logging
from lib.general_helpers.configure_loggers import get_logger
import argparse

def run_expedition(test_name: str, expedition_folder: str, logger: logging.Logger = None, threshold_fastq: int = 10):
    """
    run the consensus-identification pipeline on an expedition folder. 

    Parameters
    -----------
    test_name (str): name given to the results folder
    expedition_folder (str): name of the input expedition folder
    logger (logging.Logger, optional): logger to track execution
    threshold_fastq (int, optional): number of minimum reads in a fastq file to run the pipeline.
    """
    
    expedition_result_folder_path = ospath.join("output",test_name)
    if not ospath.exists(expedition_result_folder_path):
        os.makedirs(expedition_result_folder_path)
    
    if logger == None:
        logger = get_logger(test_name, expedition_result_folder_path, expedition_folder)
    
    expedition_folder_path = ospath.join("data", expedition_folder)
    init_db = pd.read_csv(ospath.join(expedition_folder_path, "sample_info.csv"))
    df = init_db[["Barcode", "Species", "gene"]]
    df.insert(3, "matK", ["Not sequenced"]*df.shape[0])
    df.insert(4, "rbcL", ["Not sequenced"]*df.shape[0])
    df.insert(5, "psbA-trnH", ["Not sequenced"]*df.shape[0])
    df.insert(6, "ITS", ["Not sequenced"]*df.shape[0])

    for barcode in df["Barcode"]:
        

        for root, dirs, files in os.walk(expedition_folder_path):
            if root.endswith(str(f"barcode{barcode}")):
                input_barcode_path = root

        for _,_, files in os.walk(input_barcode_path):
            for file in files:
                if file.endswith(".fastq"):
                    input_fastq_path = ospath.join(input_barcode_path, file)
                    input_fastq_filename=file
        
        too_small = True
        n=0
        for read in SeqIO.parse(input_fastq_path, "fastq"):
            n+=1
            if n >= threshold_fastq:
                too_small = False
                break
        
        if too_small:
            logging.error(f"Barcode{barcode} folder has only {n} reads, lower than the minimum threshold of {threshold_fastq} to run pipeline. Review the fastq folder or lower the threshold.")
            df.loc[df["Barcode"] == barcode, ["matK", "rbcL", "psbA-trnH", "ITS"]] = "fastq too small"
        else:
            logger.info(f"Running pipeline for Barcode{barcode}...\n")
            genes = df.loc[df.Barcode == barcode, "gene"]
            genes = genes.iloc[0]
            genes = genes.split(",")
            genes = {gene.strip() for gene in genes}
            pipeline_single(f"barcode{barcode}", input_fastq_path= input_fastq_path, db= genes, windows=False, expedition_name=test_name)

            identification_results_path = ospath.join(expedition_result_folder_path, f"barcode{barcode}" ,"identification")
            for gene in genes:
                result_tuple = get_best_species_from_xml(ospath.join(identification_results_path, f"{gene}.txt"))
                df.loc[df["Barcode"]==barcode, gene] = f"species: {result_tuple[0]}, with score {result_tuple[1][0]} and evalue {result_tuple[1][1]})"

    
    logger.info("Finished!")
    df.to_excel(ospath.join(expedition_result_folder_path, test_name+".xlsx"))
    

def main():
    run_expedition("test_jardin_botanique_on_all_genes","expedition_jardin_botanique")

if __name__ == "__main__":
    main()
