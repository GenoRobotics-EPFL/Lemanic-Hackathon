from lib.general_helpers.run_command import run_command
from lib.general_helpers.configure_loggers import *
import os
import os.path as ospath
from Bio import SeqIO
import pandas as pd
import shutil


def get_info_on_cluster(fastq_file_path):
    """"
    get information of the fastq file of one cluster generated with isONclust, number of total reads and reads associated with each gene. 
    The input fastq file should be labeled in the description of reads to include the gene it comes from.

    Parameters
    -----------
    fastq_file_path(str): input LABELED fastq path

    Returns
    ----------
    tuple: (dictionary with genes as keys and number of keys as values, number of total reads)
    """
    parser = SeqIO.parse(fastq_file_path, "fastq")
    nb_reads=0
    nb_reads_per_gene = {"matK":0,"ITS":0,"rbcL":0,"psbA-trnH":0}
    for record in parser:
        nb_reads+=1
        genes= ["matK","ITS","rbcL","psbA-trnH"]
        for gene in genes:
            if gene in record.description:
                nb_reads_per_gene[gene]+=1
    return nb_reads_per_gene,nb_reads

def get_info_on_clustering(demultiplexing_folder_path):
    """
    Get information about all the clusters in a folder generated by isONclust. folder should contain no fastq files other than those generated by isONclust

    Parameters
    -----------
    demultiplexing_folder_path: folder path

    Returns
    ------------
    final_df: sorted pandas dataframe with information about number of reads for each cluster
    """
    rows_list =[]
    for root,dirs,files in os.walk(demultiplexing_folder_path):
        for file in files:
            if file.endswith("fastq"):
                info_dict ={}
                info_dict["cluster"]= "total" if file=="sorted.fastq" else f"{file[:-6]}"
                nb_reads_per_gene,nb_reads=get_info_on_cluster(ospath.join(root,file))
                info_dict["nb_reads"]=nb_reads
                for gene in ["matK","ITS","rbcL","psbA-trnH"]:
                    info_dict[gene]=nb_reads_per_gene[gene]
                rows_list.append(info_dict)
    result_df = pd.DataFrame(rows_list)
    final_df = result_df.sort_values("nb_reads", ascending=False)
    final_df= final_df.reset_index(drop=True)
    return final_df

def get_dominating_genes_per_cluster(result_df,threshold_percentage=0.01):
    """
    return the dominating gene for each cluster, aka the gene with the largest number of reads.

    Parameters
    ----------
    result_df(pandas.Dataframe): pandas df with a summary of results of clustering, as generated by get_info_on_clustering
    threshold_percentage(optional, int = 0.01): percentage of total reads a cluster should contain to be considered in analysis

    Returns
    ----------
    Pandas.Series.Series: series with same row indexes and dominating gene as values
    """
    total_nb_reads = result_df.loc[result_df.cluster=="total","nb_reads"]
    total_nb_reads=total_nb_reads.iloc[0]
    threshold= total_nb_reads * threshold_percentage
    working_df= result_df[result_df.nb_reads>=threshold]

    return  working_df.iloc[:,2:].idxmax(axis=1)


def get_nb_clusters_per_gene(result_df,threshold_percentage=0.01):
    """
    get the number of clusters where each gene dominates

    Parameters
    ----------
    result_df(pandas.Dataframe): pandas df with a summary of results of clustering, as generated by get_info_on_clustering
    threshold_percentage(optional, int = 0.01): percentage of total reads a cluster should contain to be considered in analysis

    Returns
    ----------
    pandas.Series.Series: Series with the genes as rows and the count as values
    """
    dominating_genes= get_dominating_genes_per_cluster(result_df,threshold_percentage)
    return dominating_genes.value_counts()

def keep_best_clusters(df, demultiplexing_folder_path,threshold_percentage=0.01):
    """
    delete all of the cluster fastq files which have an insufficient number of reads. 

    Parameters
    ----------    
    result_df(pandas.Dataframe): pandas df with a summary of results of clustering, as generated by get_info_on_clustering
    demultiplexing_folder_path: path to the folder containing all cluster fastq files
    threshold_percentage(optional, int = 0.01): percentage of total reads a cluster should contain to be kept.

    """

    total_nb_reads = df.loc[df.cluster=="total","nb_reads"].loc[0]
    keep_files=[]
    for i in range(len(df)):
        if df.loc[i,"nb_reads"]>= total_nb_reads*threshold_percentage:
            cluster_number = df.loc[i,"cluster"]
            keep_files.append(f"{cluster_number}.fastq")
    for root,dirs,files in os.walk(demultiplexing_folder_path):
        for file in files:
            if file.endswith(".log"):
                keep_files.append(file)
    delete_files_except(demultiplexing_folder_path,keep_files)
        
def delete_files_except(folder_path, keep_files):
    """
    delete all files in folder except for files ending in .log and those passed as argument

    Parameters
    ----------
    folder_path(str): path to the folder to clean
    keep_files(list): list of file names, not files paths, to keep
    
    """
    # Get a list of all files and directories in the folder
    for root, dirs, files in os.walk(folder_path, topdown=False):
        for name in files:
            file_path = os.path.join(root, name)
            # Check if the file is not in the list of files to keep
            if name not in keep_files and ".log" not in name:
                # Delete the file
                os.remove(file_path)

        for name in dirs:
            dir_path = os.path.join(root, name)
            # Delete the directory and its contents
            shutil.rmtree(dir_path)


def evaluation_of_isoncluster(input_folder="output/test_on_fake_multiplex_just_isonclust"):
    # Get the list of directories in input folder
    list_of_directories = os.listdir(input_folder)

    evaluation_dictionary = {'total_nb_reads': [], 'nb_of_non_clustered_reads': [], 'number_of_clusters': [],
                             'avg_cluster_accuracy': []}

    names_of_species = []

    # Iterating on the directories
    for directory in list_of_directories:
        if directory.endswith("_mutiplexed"):  # For every directory that ends with multplexed
            final_df = get_info_on_clustering(input_folder + "/" + directory)
            final_df['dominating_gene'] = get_dominating_genes_per_cluster(final_df)
            final_df.drop(final_df[final_df['dominating_gene'].isna()].index, inplace=True)
            for index, row in final_df.iterrows():
                final_df.loc[index, 'accuracy'] = row[row['dominating_gene']] / row['nb_reads']

            evaluation_dictionary['total_nb_reads'].append(final_df['nb_reads'].iloc[0])
            evaluation_dictionary['nb_of_non_clustered_reads'].append(final_df['nb_reads'].iloc[0] -
                                                                      final_df['nb_reads'].iloc[1:].sum())
            evaluation_dictionary['number_of_clusters'].append(len(final_df) - 1)
            evaluation_dictionary['avg_cluster_accuracy'].append(final_df['accuracy'].iloc[1:].mean())
            names_of_species.append("_".join(directory.split("_")[:-1]))

    final = pd.DataFrame.from_dict(evaluation_dictionary)
    final.index = names_of_species
    return final
