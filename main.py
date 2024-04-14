

test_name= "photinia_davidiana_test_all"
input_folder = ospath.join("data","expedition_jardin_botanique","Photinia_davidiana_trnH-psbA_barcode82")
print(ospath.exists(input_folder))
total_time, time_minimap2, time_racon= run_consensus(test_name, ospath.join(input_folder,"Photinia_davidianatrnH-psbA_barcode82.fastq"), consensus_method="straightforward_all_sequences", output_dir=ospath.join(os.getcwd(), "output", test_name))
ref_seq = SeqIO.read(ospath.join(input_folder,"Photinia_davidiana_reference_seq.fasta"), "fasta")
consensus_sequences = list(SeqIO.parse(ospath.join("output", test_name, "photinia_davidiana_test_all_final_consensus.fasta"), "fasta"))
print(compute_alignment_score(consensus_sequences[0],ref_seq))
