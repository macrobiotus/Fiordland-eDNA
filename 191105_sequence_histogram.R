library("Biostrings")

fasta_input <- readDNAStringSet("/Users/paul/Documents/OU_eDNA/191105_reference_sequences/191105_Miya2015_Tab1_reference_sequences.fasta")

str(fasta_input)

lengths <- width(fasta_input)
summary(lengths) 
