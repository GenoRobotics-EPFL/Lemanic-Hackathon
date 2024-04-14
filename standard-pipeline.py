import os
import sys
from lib.consensus.consensus import run_consensus
from lib.identification.identification import run_identification
import argparse

# def the main function

def pipeline_single(input_name: str, input_fastq_path: str, db: str, windows: bool, expedition_name:str = None):

    # preprocessing

    # run_preprocessing(input_fastq_filename, input_fastq_path, output_dir)

    # consensus

    consensus_method = "80_20_best_sequence"
    run_consensus(input_name, input_fastq_path, consensus_method, expedition_name= expedition_name, windows=windows)

    # identification

    return run_identification(input_name=input_name, expedition_name= expedition_name, db=db)

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_name", required =True)
parser.add_argument("-path", "--input_fastq_path", required =True)
parser.add_argument("-exp", "--expedition_name", required =False)
parser.add_argument("-db","--db", required = False)
parser.add_argument("-w", "--windows", required = False)

def main():
    args = parser.parse_args()
    windows = True if args.windows.lower() in ["true", "t", "yes", "y"] else False
    return pipeline_single(input_name =args.input_name, input_fastq_path=os.path.normpath(args.input_fastq_path), db=args.db, windows=windows, expedition_name=args.expedition_name)

if __name__ == "__main__":
    main()