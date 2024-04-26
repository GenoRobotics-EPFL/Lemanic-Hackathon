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
    
    
from logging import getLogger
from Bio import SeqIO
from Bio.Align import sam
from lib.general_helpers.run_command import run_command
import pandas as pd
import shutil


def demultiplexing_allingnment(input_fastq, reference_paths, output_folder):
    """
    Demultiplexing the fastq file containing the reads.
    It saves one fastq file for each gene, demultiplexing the reads.

    Args:
        input_fastq: input fastq file containing the multiplexed reads
        reference_paths: paths with reference sequences for each gene
        output_folder: output folder path for saving each final fastq file
    """
    logger = getLogger('demultiplexing')

    # Get genes' names from the reference paths
    genes = []
    for path in reference_paths:
        genes.append(path.split('.')[0])

    # Deliting output directiory
    shutil.rmtree("Alignment", ignore_errors=False)

    # Create SAM files for each gene
    alignment_paths = []
    for gene, reference_path in zip(genes, reference_paths):
        alignment_path = f"Alignment/{gene}.sam"
        minimap2_command = f"minimap2 -ax map-ont {reference_path} {input_fastq} > {alignment_path}"
        _, minimap2_time = run_command(minimap2_command, logger, windows=False)
        alignment_paths.append(alignment_path)

    # Create empy dictionary to store information from SAM files
    data = {}
    for gene in genes:
        data[f'{gene}_seq_id'] = []
        data[f'{gene}_flag'] = []
        data[f'{gene}_mapq'] = []

    # Get the sequence id, flag and map quality from SAM files
    for gene, alignment_path in zip(genes, alignment_paths):
        print(alignment_path)
        # Get the iterator on the alignments
        alignment_iterator = sam.AlignmentIterator(alignment_path)
        for alignment in alignment_iterator:
            if alignment.sequences[1].id in data[f'{gene}_seq_id']:  # Check it is not duplicated
                if bin(alignment.flag)[2:].zfill(12)[0] == 1 and bin(alignment.flag)[2:].zfill(12)[3] == 1:
                    pass
            else:
                data[f'{gene}_seq_id'].append(alignment.sequences[1].id)
                data[f'{gene}_flag'].append(alignment.flag)
                data[f'{gene}_mapq'].append(alignment.mapq)

    # Create the dataframe from the dictionary
    sam_out = pd.DataFrame.from_dict(data)
    # Get the gene with the best map quality score
    sam_out['mapq_max'] = sam_out.iloc[:, 2::3].idxmax(axis=1)
    sam_out['unmapped'] = sam_out.iloc[:, 2::3].max(axis=1)

    # Opening output files for each gene
    output_files = []
    for gene in genes:
        output_files.append(open(f'{output_folder}/{gene}', 'w'))

    # Saving each read in the correct fastq file
    index = 0
    for read in SeqIO.parse(input_fastq, "fastq"):
        if sam_out.iloc[index]['unmapped'] == 0:  # Flag = 0 -> unmapped read, so we discard it
            pass
        else:
            gene_name = sam_out.iloc[index]['mapq_max'][0:-5]
            position = genes.index(gene_name)
            SeqIO.write(read, output_files[position], 'fastq')
        index += 1

    # Closing the output files
    for output_file in output_files:
        output_file.close()

    # Deliting output directiory
    shutil.rmtree("Alignment", ignore_errors=False)

