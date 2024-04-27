from lib.general_helpers.download_GenBank import download_sequence_from_id
from lib.general_helpers.generate_fastq import generate_fastq
from lib.general_helpers.generate_fastq import multiplex
import os
import os.path as ospath
import argparse
from Bio import SeqIO
import shutil

def main():

    destination_folder = "Solanum_lycopersicum_fake_data_total"
    ids= {"matK":"OQ303562.1","rbcL":"HF572813.1", "psbA-trnH":"FJ493314.1", "ITS":"KC213749.1"}
    error_rates= [0.01,0.03,0.05,0.1,0.15,0.2,0.3,0.4]
    multiplexed = True


    destination_folder_path= ospath.join("data", destination_folder)
    if ospath.exists(destination_folder_path):
        shutil.rmtree(destination_folder_path)
    if not ospath.exists(destination_folder_path):
        os.makedirs(destination_folder_path)

    record_parser = download_sequence_from_id(ids.values())
    sequences=dict()
    for gene,record in zip(ids.keys(),record_parser):
        sequences[gene] =record.seq


    for error_rate in error_rates:
        specific_error_folder= ospath.join(destination_folder_path,f"Solanum_lycopersicum_{error_rate*100:.0f}")
        os.makedirs(specific_error_folder)
        for gene in ids.keys():
            print("in")
            new_demultiplexed_fastq_path = ospath.join(specific_error_folder, str(ids[gene])+f"_{gene}" +f"_Solanum_lycopersicum_{error_rate*100:.0f}.fastq")
            generate_fastq(sequences[gene], 1000, new_demultiplexed_fastq_path,mutation_prob=error_rate)

        if multiplexed:
            multiplex_path=ospath.join(specific_error_folder,"multiplexed.fastq")
            multiplex([specific_error_folder],multiplex_path)


if __name__ == "__main__":
    main()