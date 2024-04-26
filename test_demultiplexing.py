from lib.demultiplexing.run_demultiplexing import run_demultiplexing
from lib.general_helpers.run_command import run_command
from lib.general_helpers.configure_loggers import *
from lib.demultiplexing.demultiplexing_helpers import *
import os
import os.path as ospath
from Bio import SeqIO
import pandas as pd

#logger = get_logger("my_logger",os.getcwd(),"any")
#run_demultiplexing("test_first","/home/ghassan_unix/genorobotics/Lemanic-Hackathon/data/fake_multiplexed/Allium_Ursinum_mutiplexed.fastq",ospath.join("output","first"),logger=logger)





def run_test_demultiplexing():
    input_multiplex_path = ospath.join("data","fake_multiplexed")
    output_multiplex_path = ospath.join("output","test_on_fake_multiplex")
    if not ospath.exists(output_multiplex_path):
        os.makedirs(output_multiplex_path)
    logger = get_logger("demultiplexing",output_multiplex_path,"test_on_fake_multiplex")

    for root,dirs,files in os.walk(input_multiplex_path):
        for file in files:
            if file.endswith(".fastq"):
                run_demultiplexing(file[:-6],ospath.join(root,file),ospath.join(output_multiplex_path,file[:-6]))
    """
    for root,dirs,files in os.walk(output_multiplex_path)      :
        for dir in dirs:
            print(dir)
            df = get_info_on_clustering(ospath.join(root,dir))
            print(df)
            print(get_nb_clusters_per_gene(df))
        break
    """

run_test_demultiplexing()

    
            









        
