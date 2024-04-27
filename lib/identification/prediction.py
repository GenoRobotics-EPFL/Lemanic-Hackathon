import math

from Bio.Blast import NCBIXML
from Bio import Entrez
import pandas as pd
from Bio import SeqIO
from cmath import log10, sqrt


def read_blastn(xmlName):
    """read the output of a blastn file and get X best result
            for now gets all the best results
            the output format is:
            [id code, id name, score, e value]"""

    result_handle = open(xmlName, 'r')
    blast_records = NCBIXML.parse(result_handle)
    sequences = []
    codes = []

    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                sequence_id = alignment.title  # Extract sequence ID
                score = hsp.score  # score
                e_value = hsp.expect  # E value

                parts = sequence_id.split()
                sequenceCode = parts[0]  # accesion to id
                sequenceDescription = ""
                codes.append(sequenceCode)
                sequences.append([sequenceCode, sequenceDescription, score, e_value])

    codesToSpecies = get_species_from_ids(codes)
    for i, pairs in enumerate(codesToSpecies):
        sequences[i][1] = pairs[1]

    return sequences


def get_species_name_from_gene_id(gene_id):  ## adapt
    Entrez.email = "ilgazarslanplus@gmail.com"

    handle = Entrez.efetch(db="nuccore", id=gene_id, rettype="gb", retmode="xml")
    record = Entrez.read(handle)
    handle.close()
    return record[0]['GBSeq_organism']


def get_species_from_ids(ids: list):
    """
    Return the species from GenBank entries with certain ids.

    Parameters
    ----------
    ids: list of ids to search

    Returns
    ----------
    dictionary that matches ids to organism names

    """
    Entrez.email = "ghassan.abboud@epfl.ch"
    species_list = list()
    handle = Entrez.efetch(db="nuccore", id=ids, retmode="txt", rettype="gb")
    records = SeqIO.parse(handle, "genbank")
    for id in ids:
        one_record = next(records)
        species = one_record.annotations["organism"]
        species_list.append([id, species])
    return species_list


def get_best_species(*gene_data_parameter):
    """read BLASTn outputs for the 4 genes and get a list of the best species with e-values and percent identities.
    return a pandas dataframe with all the best_species with their scores"""

    # Create an empty DataFrame
    df = pd.DataFrame()
    pd.set_option('display.float_format', '{:.2e}'.format)
    weightDict = {}

    # Iterate over each gene and populate the DataFrame
    for gene_data_item in gene_data_parameter:
        (gene, gene_data, weight) = (gene_data_item[0], gene_data_item[1], len(gene_data_item[1]))
        gene_df = pd.DataFrame(gene_data, columns=["id", "species", "score", "e_value"])
        # gene_df = gene_df[["species", "score", "e_value"]] ,this line removes the id part !!!!
        gene_df.set_index("species", inplace=True)
        gene_df.columns = [(gene, col) for col in gene_df.columns]
        weightDict[gene, "e_value"] = weight  # the amount of genes for each gene

        df = pd.merge(df, gene_df, left_index=True, right_index=True, how="outer")

    total_weight = sum(weightDict.values())
    for key, value in weightDict.items():
        if value == 0:
            weightDict[key] = 1
        else:
            weightDict[key] = value / total_weight
    #weightDict = {key: value / total_weight for key, value in weightDict.items()}

    print(weightDict)
    scoreColumns = df.columns[1::3].tolist()
    eValueColumns = df.columns[2::3].tolist()

    for column in scoreColumns:
        minScore = df[column].min()
        df[column].fillna(minScore, inplace=True)

    for column in eValueColumns:
        maxE_value = df[column].max()
        df[column].fillna(maxE_value, inplace=True)

    return df, weightDict


def calculate_metric(pandasFrame, weightDict):
    """take as input the output of get_best_species and calculate for each species a metric that determines how close
    it is."""

    score_columns = pandasFrame.iloc[:, 1::3]  # scores (1 mod3 columns)
    # e_value_columns = pandasFrame.iloc[:, 2::3]  # evalues (2 mod 3 columns)

    pandasFrame['score_average'] = score_columns.mean(axis=1)

    for col, weight in weightDict.items():
        pandasFrame[col] *= weight
    pandasFrame['e_value_weighted'] = pandasFrame[list(weightDict.keys())].sum(axis=1)

    #pandasFrame['e_value_p'] = pandasFrame['e_value_weighted'].apply(lambda x: (-1 * log10(x)).real)

    pandasFrame.sort_values(by=['species', 'e_value_weighted'], inplace=True) # sorts ascending order
    pandasFrame = pandasFrame[~pandasFrame.index.duplicated(keep='first')]

    return pandasFrame, pandasFrame["score_average"].tolist(), pandasFrame["e_value_weighted"].tolist()


def give_prediction(*gene_data_files):
    """call all functions and give the best species"""

    # put them in the needed format and take the worst case for each
    gene_data_list = []

    for gene_data_file in gene_data_files:
        gene_name = gene_data_file.split("_")[0].upper()
        gene_data_list.append((gene_name, read_blastn(gene_data_file)))

    dataframe, weightDict = get_best_species(*gene_data_list)
    final_frame, scores, e_value = calculate_metric(dataframe, weightDict)

    names = final_frame.index.tolist()

    bestSpecie = ""

    if (len(set(e_value))) != 1:  # e values are different

        while len(set(e_value)) != 1:
            min_evalue_gene = min(e_value)
            indices_to_remove = [i for i, e_val in enumerate(e_value) if e_val != min_evalue_gene]
            e_value = [e_val for i, e_val in enumerate(e_value) if i not in indices_to_remove]
            scores = [score for i, score in enumerate(scores) if i not in indices_to_remove]
            names = [name for i, name in enumerate(names) if i not in indices_to_remove]

        minpos = e_value.index(min(e_value))
        bestSpecie = names[minpos]

    if (len(set(e_value))) == 1:  # e values are same, look at the score
        maxpos = scores.index(max(scores))
        bestSpecie = names[maxpos]

    return f'The best match is with the species : {bestSpecie}'


print(give_prediction("psbA-trnH.txt","rbcL.txt"))
