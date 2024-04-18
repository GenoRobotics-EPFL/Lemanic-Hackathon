# Lemanic-Hackathon

This is working repository for Genorobotic's participation to the Lemanic Life Sciences Hackathon of April 2024. 

## Genorobotics in brief

Genorobotics aims to develop field tools for plant biodiversity identification based on DNA barcoding. This eliminates the need to transport biological samples to labs for testing and aims to support biodiversity conservation efforts.

Four genes are used to identify the species of a sample: matk, trnH-psbA, ITS and rbcL. These genes are amplified through Polymerase Chain Reaction (PCR) and sequenced with an Oxford Nanopore Portable Sequencer.

The goal of the bioinformatics team is to interpret the data generated by the sequencer. First, the raw DNA reads from a fastq file need to be aligned to generate a consensus sequence for each gene. Then, this sequence is compared to NCBI's database of genetic sequences, GenBank with BLASTn. The results of the BLASTn queries for the four genes are put together to predict the species.

<p align="center">
    <img src="images/Photomontage_jungle.jpg" alt="montage" width="500"/>
</p>

## The datasets

The datasets used are the sequencing results of two Genorobotics expeditions:
- [An expedition to Lausanne's botanical garden](data/expedition_jardin_botanique/) consisting of 27 samples of 14 various species and the 4 barcoding genes. The fastq reads are not multiplexed by genes and the DNA was collected using Genorobotics' lower fidelity microneedle patches.
- [A summer expedition to ?](data/summer_expedition/) consisting of 12 samples of 12 different plants. The reads are multiplexed by genes and the DNA was collected with a higher fidelity Qiagen DNeasy kit.

### data organization
- Each sub-folder represents a sample and spells out the sample's species, genes sequenced and barcode used.
- Each subfolder contains a `fastq` file containing the raw reads and a `fasta` file containing the species reference sequences for the genes amplified from the GenBank database
- Three csv files `general_info.csv`, `primer_info.csv` and `sample_info.csv` contain information about the expedition, the primers used and the species/genes for each sample respectively.

### data generation
Yair's code

## The project

The project revolves around DNA sequence clustering. This can be used for 2 purposes:
- **Demultiplexing without additional barcodes**: In the summer expedition, multiple genes are sequenced for every plant, all within one sample. This reduces the number of pipetting steps on site and tubes used. It also simplifies the protocole, reducing manipulation errors. However, This means the reads from multiple different genes end up in the same file. This is problematic for consensus sequence generation, as no one sequence will emerge.
- **Increasing coverage depth by directionalizing reads**: Even when one gene is sequenced, there are four kinds of reads in one file. This is because PCR amplifies both coding and non-coding strands, and the direction in which the sequence is inserted in the sequencer is random. The four orientiations are then: coding 5', coding 3', non-coding 5', non-coding 3'. 