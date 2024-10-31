# Viral Sequence Analysis Pipeline in Python

# Import necessary libraries
from scripts.utils import fetch_viral_sequences, save_sequences_to_file, run_cd_hit, run_blast, create_phylogenetic_tree_mafft, create_phylogenetic_tree_muscle, visualize_tree, plot_sequence_length_distribution, run_functional_annotation
import pandas as pd
import matplotlib.pyplot as plt

# Pipeline Execution

# Set parameters
search_terms = ["unclassified virus", "unknown virus", "environmental virus"]
email = "sebastien.delandtsheer@example.com"
retmax = 500

# Fetch sequences
sequences = fetch_viral_sequences(search_terms, email, retmax)

# Save sequences to file
save_sequences_to_file(sequences, "data/raw/viral_sequences.fasta")

# Run CD-HIT for clustering
run_cd_hit("data/raw/viral_sequences.fasta", "data/processed/filtered_sequences.fasta")

# Run BLAST for classification
run_blast("data/processed/filtered_sequences.fasta", "results/blast/blast_results.out")

# Create a phylogenetic tree using MAFFT and FastTree
create_phylogenetic_tree_mafft("data/processed/filtered_sequences.fasta", "results/alignments/aligned_sequences_mafft.fasta", "results/trees/tree_file_mafft.nwk")

# Alternatively, create a phylogenetic tree using MUSCLE and FastTree
#create_phylogenetic_tree_muscle("filtered_sequences.fasta", "aligned_sequences_muscle.fasta", "tree_file_muscle.nwk")

# Visualize the tree
visualize_tree("results/trees/tree_file_mafft.nwk")

# Assuming metadata file exists (can be generated manually based on sequence data)
metadata_file = "sequence_metadata.csv"
plot_sequence_length_distribution(metadata_file)

# Run functional annotation
run_functional_annotation("data/processed/filtered_sequences.fasta")

# End of the pipeline
