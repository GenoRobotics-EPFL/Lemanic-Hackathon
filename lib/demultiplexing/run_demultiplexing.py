from lib.general_helpers.configure_loggers import configure_demultiplexing_logger
from lib.demultiplexing.demultiplexing_pipelines.demultiplexing_isONclust_pipeline import isONclust_pipeline
from lib.general_helpers.configure_loggers import *
from lib.identification.identification import run_identification
from lib.identification.identification_pipelines.identification_processing import *
from lib.general_helpers.download_GenBank import *
from lib.demultiplexing.demultiplexing_helpers import delete_files_except
import os
import os.path as ospath
from logging import getLogger
from Bio import SeqIO
from Bio.Align import sam
from lib.general_helpers.run_command import run_command
from lib.consensus.consensus import run_consensus_pipeline_80_20_best_sequence
import pandas as pd
import shutil

def run_demultiplexing(input_name:str, input_file_path: str, output_folder: str, demultiplexing_method: str ="complete", logger = None, genes=["matK","rbcL","psbA-trnH","ITS"]):
    
    if logger == None:
        logger = configure_demultiplexing_logger(output_folder, input_name)
        print(f"Logging set up at {output_folder}/{input_name}_demultiplexing_pipeline_log.log")
    
    if demultiplexing_method == "isONclust":
        logger.info(f"Running demultiplexing pipeline on {input_file_path} using the {demultiplexing_method} method")
        isONclust_pipeline(input_file_path,output_folder, logger=logger)
    
    if demultiplexing_method == "complete":
        logger.info(f"Running demultiplexing pipeline on {input_file_path} using the {demultiplexing_method} method")
        logger.info(f"Launching first step: isONclust")
        isONclust_pipeline(input_file_path,output_folder, logger=logger)
        cluster_files = os.listdir(output_folder)
        for file in cluster_files:
            if file.endswith(".fastq"):
                logger.info(f"Running consensus pipeline on {ospath.join(output_folder,file)}...")
                cluster_consensus_output_folder=ospath.join(output_folder,f"{file[:-6]}","consensus")
                os.makedirs(cluster_consensus_output_folder,exist_ok=True)
                run_consensus_pipeline_80_20_best_sequence(f"{file[:-6]}",ospath.join(output_folder,file),cluster_consensus_output_folder,logger)
        
        for root,dirs,_ in os.walk(output_folder):
            for dir in dirs:
                consensus_folder_path = ospath.join(root,dir,"consensus")
                files= os.listdir(consensus_folder_path)
                for file in files:
                    if "final_consensus" in file and not ospath.getsize(ospath.join(consensus_folder_path,file))==0:
                        input_name = ospath.basename(dir)
                        final_consensus_path = ospath.join(consensus_folder_path)
                        identification_path = ospath.join(root,dir,"identification")
                        os.makedirs(identification_path,exist_ok=True)
                        run_identification(input_name,input_path=final_consensus_path, db=genes,logger=logger,output_dir=identification_path)
            break

        
        reference_seq_folder_path= ospath.join(output_folder,"reference_seq")
        os.makedirs(reference_seq_folder_path)
        logger.info("selecting best sequences to use as reference...")
        logger.info(f"genes of interest are {genes}")
        for gene in genes:
            best_sequence_id= None
            best_evalue=100
            best_alignment_score=0
            for root,dirs,files in os.walk(output_folder):
                for file in files:
                    if gene in file:
                        current_best_id,current_best_evalue,current_best_alignment_score =get_best_seqid_from_blastn_xml(ospath.join(root,file))
                        if current_best_evalue < best_evalue or (current_best_evalue == best_evalue and current_best_alignment_score > best_alignment_score):
                            best_evalue=current_best_evalue
                            best_sequence_id=current_best_id
                            best_alignment_score=current_best_alignment_score
            logger.info(f"Best sequence for {gene} has id {best_sequence_id}, with e-value {best_evalue} and alignment score {best_alignment_score}.")
            logger.info("Downloading it from NCBI through Entrez...")
            download_sequence_from_id([best_sequence_id],ospath.join(reference_seq_folder_path,f"{gene}.fasta"))
        
        
        reference_paths=[]
        ref_files= os.listdir(reference_seq_folder_path)
        for ref_file in ref_files:
            if ref_file.endswith(".fasta"):
                reference_paths.append(ospath.join(reference_seq_folder_path,ref_file))
        logger.info("Launching second step: Alignment-based classification")
        demultiplexing_alignment(input_file_path,reference_paths,output_folder,logger)
    
    


def demultiplexing_alignment(input_fastq, reference_paths, output_folder,logger):
    """
    Demultiplexing the fastq file containing the reads.
    It saves one fastq file for each gene, demultiplexing the reads. CAREFUL SAM =0

    Args:
        input_fastq: input fastq file containing the multiplexed reads
        reference_paths: paths with reference sequences for each gene
        output_folder: output folder path for saving each final fastq file
    """
    # Get genes' names from the reference paths
    genes = []
    for path in reference_paths:
        genes.append(ospath.basename(path).split('.')[0])
    # Create SAM files for each gene
    alignment_paths = []
    os.makedirs(ospath.join(output_folder,"Alignment"))
    for gene, reference_path in zip(genes, reference_paths):
        logger.info(f"Mapping input fastq file {input_fastq} to reference sequence {reference_path}")
        alignment_path = ospath.join(output_folder,f"Alignment/{gene}.sam")
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
        output_files.append(open(f'{output_folder}/{gene}.fastq', 'w'))

    # Saving each read in the correct fastq file
    logger.info("Seperating reads based on alignment...")
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
    files_to_keep=[]
    for gene in genes:
        files_to_keep.append(f'{gene}.fastq')
    delete_files_except(output_folder,files_to_keep)
    
    #shutil.rmtree(ospath.join(output_folder,"Alignment"), ignore_errors=False)

