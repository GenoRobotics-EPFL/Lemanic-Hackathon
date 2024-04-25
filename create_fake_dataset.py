from lib.general_helpers.download_GenBank import download_sequence_from_id
from lib.general_helpers.generate_fastq import generate_fastq
from lib.general_helpers.generate_fastq import multiplex
import os
import os.path as ospath
import argparse
from Bio import SeqIO

def main():

    destination_folder = "test_fake4"
    ids= ["ON409953.1", "AF288129.1", "AF288126.1", "AF288125.1" ]
    multiplexed = True



    destination_folder_path= ospath.join("data", destination_folder)
    if not ospath.exists(destination_folder_path):
        os.makedirs(destination_folder_path)

    record_parser = download_sequence_from_id(ids)

    for id, record in zip(ids,record_parser):
        new_demultiplexed_fastq_path = ospath.join(destination_folder_path, str(id)+".fastq")
        generate_fastq(record.seq, 100, new_demultiplexed_fastq_path)

    if multiplexed:
        multiplex_path=ospath.join(destination_folder_path,"multiplexed.fastq")
        multiplex([destination_folder_path],multiplex_path)


if __name__ == "__main__":
    main()