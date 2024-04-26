from logging import getLogger
from Bio import SeqIO
from Bio.Align import sam
from lib.general_helpers.run_command import run_command
import pandas as pd

# matK = SeqIO.read("Sequences/matK.fasta", "fasta").seq
# rbcL = SeqIO.read("Sequences/rbcL.fasta", "fasta").seq
# psbA_trnH = SeqIO.read("Sequences/psbA-trnH.fasta", "fasta").seq
# ITS = SeqIO.read("Sequences/ITS.fasta", "fasta").seq

count = 0
for record in SeqIO.parse("../../data/fake_multiplexed/Ficus_religiosa_mutiplexed.fastq", "fastq"):
    count += 1
print(count)

logger = getLogger('ciao')

genes = ['matK', 'rbcL', 'psbA-trnH', 'ITS']

alignment_paths = []

for gene in genes:
    reference_path = f"Sequences/Ficus_religiosa/{gene}.fasta"
    read_path = "../../data/fake_multiplexed/Ficus_religiosa_mutiplexed.fastq"
    alignment_path = f"Alignment/Ficus_religiosa/{gene}.sam"
    minimap2_command = f"minimap2 -ax map-ont {reference_path} {read_path} > {alignment_path}"
    _, minimap2_time = run_command(minimap2_command, logger, windows=False)
    alignment_paths.append(alignment_path)


data = {'matK_seq_id': [], 'matK_flag': [], 'matK_mapq': [], 'rbcL_seq_id': [], 'rbcL_flag': [], 'rbcL_mapq': [],
        'psbA-trnH_seq_id': [], 'psbA-trnH_flag': [], 'psbA-trnH_mapq': [],
        'ITS_seq_id': [], 'ITS_flag': [], 'ITS_mapq': []}

c = 0
for alignment_path in alignment_paths:
    print(alignment_path)
    alignment_iterator = sam.AlignmentIterator(alignment_path)
    l = 0
    for alignment in alignment_iterator:
        if alignment.sequences[1].id in data[f'{genes[c]}_seq_id']:
            if bin(alignment.flag)[2:].zfill(12)[0] == 1 and bin(alignment.flag)[2:].zfill(12)[3] == 1:
                pass
        else:
            data[f'{genes[c]}_seq_id'].append(alignment.sequences[1].id)
            data[f'{genes[c]}_flag'].append(alignment.flag)
            data[f'{genes[c]}_mapq'].append(alignment.mapq)
    print(l + count)
    c += 1

ciao = pd.DataFrame.from_dict(data)
print(ciao.head())


