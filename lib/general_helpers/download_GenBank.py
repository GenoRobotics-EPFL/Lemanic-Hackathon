from Bio import Entrez
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
Entrez.email = "ghassan.abboud@epfl.ch"

def download_sequence(species, gene_name, dst = None, start_length=None, stop_length= None, permissive_search = True):
    """
    download sequence from GenBank through the Entrez database. 

    Parameters:
    ----------
    species(str): name of species
    gene_name(str): name of gene
    dst(str,Path-like): destination file path, if None does not write to a file, only returns 
    start_length(int): minimum length of sequence
    stop_length(int): maximum length of sequence
    id(list): list of NCBi ids of sequences to download. If provided, overrides gene_name and species.
    permissive_search(bool, default = True): when True, if Advanced NCBI query returns nothing, replace it with a less precise general query.
    
    Returns
    -----------
    Bio.SeqRecord.SeqRecord: sequence found from GenBank

    """
    

    search_term = f"{gene_name}[Gene Name] AND {species}[Organism]"
    if start_length!= None or stop_length!= None:
            search_term += f" {start_length}:{stop_length}[Sequence Length]"
    handle = Entrez.esearch(db = "nucleotide", term = search_term, retmax = 1, idtype = "acc")
    search_result = Entrez.read(handle)
    handle.close()
    id = search_result["IdList"]
    if len(id)==0 and permissive_search:
        search_term = f"{gene_name} {species} {start_length}:{stop_length}[Sequence Length]"
        handle = Entrez.esearch(db = "nucleotide", term = search_term, retmax = 1, idtype = "acc")
        search_result = Entrez.read(handle)
        handle.close()
        id = search_result["IdList"]
    if len(id)==1:
        return next(download_sequence_from_id(id,dst))
    else: 
         raise ValueError("Your search did not return any results")

def download_sequence_from_id(id: list, dst = None):
    """
    download sequence from GenBank through the Entrez database using unique identifiers.

    Parameters
    ----------
    id (list): list of IDs of sequences to download from GenBank
    dst(str): destination file path

    returns
    ----------
    Bio.SeqIO.FastaIO.FastaIterator: iterator over the records corresponding to the input ids
    
    """
    handle = Entrez.efetch(db="nucleotide", id=id, retmode = "fasta", rettype = "fasta")
    record_parser = SeqIO.parse(handle, "fasta")
    if dst != None:
        with open(dst, mode="a") as writer:
            SeqIO.write(record_parser,writer,"fasta")
        return SeqIO.parse(dst,"fasta")
    return record_parser

def get_species_from_ids(ids:list):
    """
    Return the species from GenBank entries with certain ids.

    Parameters
    ----------
    ids: list of ids to search

    Returns
    ----------
    dictionary that matches ids to organism names
    
    """
    species_dict= dict()
    handle = Entrez.efetch(db="nucleotide", id=ids, retmode = "gb", rettype = "gb")
    records = SeqIO.parse(handle,"genbank")
    for id in ids:
        one_record = next(records)
        species_dict[id]=one_record.annotations["organism"]
    return species_dict
         

     