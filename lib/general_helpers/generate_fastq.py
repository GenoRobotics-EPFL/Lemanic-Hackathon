from Bio import SeqIO
import random
from Bio.Seq import Seq
import numpy as np
import os
import os.path as ospath
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt


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
    variants = ["initial_sequence", "reverse_complement", "complement", "reverse_initial"]
    variant_sequences = [sequence, sequence.reverse_complement(), sequence.complement(), sequence[::-1]]
    
    selected_variant_index = random.choices(range(len(variants)), weights=weights, k=1)[0]
    selected_variant = variant_sequences[selected_variant_index]
    variant_description = variants[selected_variant_index]

    return selected_variant, variant_description



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

def bimodal_distribution_sampler(sequence_length):
    """
    draw random position from a biomodal probability function with peaks at a third and two thirds of the sequence length
    
    Parameters
    ----------
    sequence_length(int): length of the sequence

    Returns
    ----------
    position
    """
    peaks = [sequence_length/3, 2*sequence_length/3]
    choice = np.random.random()
    if choice < 0.5:
        pos = np.random.normal(peaks[0],sequence_length/15)
    else:
        pos = np.random.normal(peaks[1],sequence_length/15)
    return max(min(round(pos),sequence_length-1),0)


def bimodal_distribution_pdf(position,sequence_length):
    """
    Probability density function of a bimodal normal distribution, highest chance is centered around positions at one or two thirds of sequence length

    Parameters
    ----------
    position: position in the sequence for which probability of event is calculated
    sequence_length: length of sequence

    Returns
    ---------
    probability of event, between 0 and 1
    """
    peaks = [sequence_length/3, 2*sequence_length/3]
    return 0.5* scipy.stats.norm(peaks[0],sequence_length/15).pdf(position) + 0.5*scipy.stats.norm(peaks[1],sequence_length/15).pdf(position)


def break_sequence_with_probability(sequence, break_position_sampler = bimodal_distribution_sampler, nbreaks = 2, take_longest= False):
    """
    Simulates breakages in the sequence, which could happen using the MinIon tool. IF YOU CHANGE THE NUMBER OF BREAKS; YOU MIGHT WANT TO CHANGE THE BREAK POSITION SAMPLER

    Arguments:
        sequence(Seq): sequence of interest
        break_prob_funtion(function): probability function sampler for determining position of breaks based on the sequence length

    Returns:
        final_sequence (Seq): final sequence after breakages
        break_info (dict): dictionary containing information about the breakages
    """
    new_sequence=sequence
    break_info = {'number_of_breaks': nbreaks, 'part_taken': ''}
    
    for i in range(nbreaks):
        position = break_position_sampler(len(sequence))
        new_sequence = new_sequence[:position+i] +"N" +new_sequence[position+i:]

    broken_parts = new_sequence.split("N")

    if take_longest:
        length_index = lambda index: len(broken_parts[index])
        final_seq_index = max(range(len(broken_parts)), key=length_index)
    else:
        final_seq_index = random.choice(range(len(broken_parts)))
    
    final_seq=broken_parts[final_seq_index]
    break_info["part_taken"] = final_seq_index
    
    return final_seq, break_info


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
    if base not in bases:
        raise ValueError("Invalid base provided. Base must be one of 'A', 'T', 'C', or 'G'.")
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
    base64_table = """!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~"""
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

def generate_fastq(input_seq, n_reads, dst, variants_weights= [0.25, 0.25,0.25,0.25], break_position_sampler=bimodal_distribution_sampler, mutation_prob =0.1, quality_score_params= [10,2,50,5], nbreaks=2, take_longest=False):

    with open(dst, "a") as handle: # Operations on newly created or updated file
        for i in range(n_reads) :
            seq, orientation = generate_weighted_sequence_variant(input_seq,variants_weights)
            broken_seq, break_info = break_sequence_with_probability(seq,bimodal_distribution_sampler, nbreaks,take_longest)
            quality_string, n_mutations,mutations_positions = assign_quality_scores(broken_seq,mutation_prob,quality_score_params[0],quality_score_params[1],quality_score_params[2],quality_score_params[3])

            variant_str = "Orientation_of_sequence="+orientation+" "
            breakage_str = "Number_of_breakages="+str(break_info['number_of_breaks'])+" Part_taken="+str(break_info['part_taken'])+" "
            mutations_str = "Number_of_mutations="+str(n_mutations)+" Positions="+str(mutations_positions)

            handle.write(f"@read{i} "+variant_str+breakage_str+mutations_str+'\n') #writes the sequence id
            handle.write(str(broken_seq)+ '\n') #writes the sequence generated by transforming it into string
            handle.write("+\n") #writes the + separator
            handle.write(quality_string + '\n') #writes the quality score

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






