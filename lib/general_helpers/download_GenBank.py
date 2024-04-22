from Bio import Entrez
Entrez.email = "ghassan.abboud@epfl.ch"

def download_sequence(species, gene_name, dst, start_length=None, stop_length= None, permissive_search = True):
    """
    download sequence from GenBank through the Entrez database. 

    Parameters:
    ----------
    species(str): name of species
    gene_name(str): name of gene
    dst(str,Path-like): destination file path
    start_length(int): minimum length of sequence
    stop_length(int): maximum length of sequence
    id(list): list of NCBi ids of sequences to download. If provided, overrides gene_name and species.
    permissive_search(bool, default = True): when True, if Advanced NCBI query returns nothing, replace it with a less precise general query.
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
        download_sequence_from_id(id,dst)

def download_sequence_from_id(id: list, dst):
    """
    download sequence from GenBank through the Entrez database using unique identifiers.

    Parameters
    ----------
    id (list): list of IDs of sequences to download from GenBank
    dst: destination file path
    
    """
    handle = Entrez.efetch(db="nucleotide", id=id, retmode = "fasta", rettype = "fasta")
    sequence = handle.read()
    handle.close()
    with open(dst, mode="a") as writer:
        writer.write(sequence)
     