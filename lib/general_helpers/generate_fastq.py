from Bio import SeqIO
import random
from Bio.Seq import Seq
import numpy as np
import os
import os.path as ospath



def generate_weighted_sequence_variant(sequence, weights=[0.25, 0.25, 0.25, 0.25]):
    """
    Flips an initial sequence as it could happen when using the MinIon tool

    Arguments:
    sequence (Seq): sequence of interest
    weights(list of 4 float): probability of each sequence version (initial sequence, reverse complement, complement, reverse initial)

    Returns:
        selected_variant(Seq): randomly reversed sequence
        variant_description(str): string indicating which variant is selected
    """
    variants = ["initial sequence", "reverse complement", "complement", "reverse initial"]
    variant_sequences = [sequence, sequence.reverse_complement(), sequence.complement(), sequence[::-1]]
    
    selected_variant_index = random.choices(range(len(variants)), weights=weights, k=1)[0]
    selected_variant = variant_sequences[selected_variant_index]
    variant_description = variants[selected_variant_index]

    return selected_variant, variant_description



def break_sequence_with_probability(sequence, break_prob_function, nbreaks = 2):
    """
    Simulates breakages in the sequence, which could happen using the MinIon tool

    Arguments:
        sequence(Seq): sequence of interest
        break_prob_funtion(function): probability function used as weight for the breakages

    Returns:
        final_sequence (Seq): final sequence after breakages
        break_info (dict): dictionary containing information about the breakages
    """
    broken_sequence = []
    break_info = {'number_of_breaks': 0, 'part_taken': ''}
    
    for i, base in enumerate(sequence):
        # Computes the breakage probability
        break_prob = break_prob_function(i, len(sequence)) / 10

        # Checks if sequence breaks at this iteration
        if random.random() <= break_prob:
            if nbreaks > 0:
                broken_sequence.append('N')
                nbreaks -= 1   # Breakage marker is 'N'
                break_info['number_of_breaks'] += 1
        else:
            broken_sequence.append(base)

    if broken_sequence.count('N') == 2:  # If 2 breakages, take the middle part
        start_index = broken_sequence.index('N')
        end_index = len(broken_sequence) - broken_sequence[::-1].index('N') - 1
        final_sequence = broken_sequence[start_index + 1:end_index]
        break_info['part_taken'] = 'middle'
    elif broken_sequence.count('N') == 1:  # If 1 breakage, take the biggest part
        n_index = broken_sequence.index('N')
        if n_index < len(broken_sequence) - n_index - 1:
            final_sequence = broken_sequence[n_index + 1:]
            break_info['part_taken'] = 'end'
        else:
            final_sequence = broken_sequence[:n_index]
            break_info['part_taken'] = 'start'
    else:  # If 0 breakages, take full sequence
        final_sequence = broken_sequence
        break_info['part_taken'] = 'full sequence'

    return final_sequence, break_info


def break_prob_function(position, sequence_length):
    """
    Example of probability function usable for the previous function

    Arguments:
        position(int): position in sequence
        sequence_length(int)
    Returns:
        probability at a specific position in sequence(float)
    """
    max_prob = 0.5  # Max probability at start and end of sequence
    return max_prob * (position / sequence_length)

def mutate(base):
    """
    Simulates substitution mutation

    Arguments:
        base(char)

    Returns:
        muted base (char)
    """
    # Bases list
    bases = ['A', 'T', 'C', 'G']
    bases.remove(base)
    # Randomly selected new base
    return random.choice(bases)


def assign_quality_scores(sequence, mutation_probability=0.1,mutation_mean = 16,mutation_sd = 3,basic_mean = 48,basic_sd = 3):
    """
    Assigns fake scores for each base, following a different normal law if base is mutated or not
    
    Arguments:
        sequence(Seq): sequence of interest
        mutation_probability(float)
        mutation_mean(int): mean of normal law associated with mutation
        mutation_sd(int): standard deviation of normal law associated with mutation
        basic_mean(int): mean of normal law associated with no-change bases
        mutation_mean(int): standard deviation of normal law associated with no-change bases
    
    Returns:
        quality score of sequence (in base 64 characters -> ASCII)
        nb_mutations(int): number of mutations
        mutations_positions(array of ints): positions of mutations
    """
    quality_scores = []
    base64_table = "!°#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~"
    nb_mutations=0
    mutations_positions= []
    for i, base in enumerate(sequence):
        # Computes score for each base depending on mutation
        if random.random() < mutation_probability:
            # If mutation needed
            base = mutate(base)
            nb_mutations+=1
            mutations_positions.append(i+1)
            quality_score = np.random.normal(mutation_mean, mutation_sd)
        else:
            # If no mutation needed
            quality_score = np.random.normal(basic_mean, basic_sd)
        # Limits the score between 0 and 93 (ASCII go from 33 to 126)
        quality_score = max(min(quality_score, 93), 0)
        base64_score = base64_table[int(round(quality_score))]
        quality_scores.append(base64_score)
    return ''.join(quality_scores), nb_mutations, mutations_positions


def multiplex(folder_paths: list, dst: str):
    """
    merge multiple fastq files into one, while writing in the description the origin of each read

    Parameters
    ----------
    folder_paths: list of folder paths to look for fastq files for merging
    dst: path of the output merged fastq file

    """
    parsers=dict()
    for folder_path in folder_paths:
        for root,_,files in os.walk(folder_path):
            for file in files:
                if file.endswith(".fastq"):
                    file_path = ospath.join(root,file)
                    parsers[file_path] = list(SeqIO.parse(file_path, "fastq"))
    
    with open(dst, "a") as writer:
        while parsers:
            parent_file_path = random.choice(list(parsers.keys()))
            current_list=parsers[parent_file_path]
            index= random.randint(0, len(current_list)-1)
            seqrecord = current_list[index]
            seqrecord.description = seqrecord.description + " source=" + parent_file_path
            SeqIO.write(seqrecord, writer, "fastq")
            current_list.pop(index)
            if not current_list:
                parsers.pop(parent_file_path)






