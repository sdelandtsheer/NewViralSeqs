# Viral Sequence Analysis Pipeline in Python

# Import necessary libraries
import subprocess
from Bio import Entrez, SeqIO, Phylo
import pandas as pd
import matplotlib.pyplot as plt


# Function Definitions

# 1. Fetch viral sequences from GenBank
def fetch_viral_sequences(search_terms, email, retmax=100):
    """
    Fetch viral sequences from GenBank.

    Parameters:
    search_terms (list): List of terms to search in GenBank (e.g., ["unclassified virus", "unknown virus"]).
    email (str): Email for NCBI Entrez identification.
    retmax (int): Maximum number of records to retrieve per term.

    Returns:
    list: List of SeqRecord objects representing viral sequences.
    """
    Entrez.email = email
    sequences = []
    for term in search_terms:
        handle = Entrez.esearch(db="nucleotide", term=term, retmax=retmax)
        record = Entrez.read(handle)
        ids = record['IdList']
        handle = Entrez.efetch(db="nucleotide", id=ids, rettype="gb", retmode="text")
        sequences.extend(list(SeqIO.parse(handle, "genbank")))
    return sequences


# 2. Save sequences to file
def save_sequences_to_file(sequences, output_file):
    """
    Save sequences to a FASTA file.

    Parameters:
    sequences (list): List of SeqRecord objects representing viral sequences.
    output_file (str): Path to the output FASTA file.
    """
    SeqIO.write(sequences, output_file, "fasta")


# 3. Run CD-HIT for Quality Check
def run_cd_hit(input_file, output_file, identity_cutoff=0.9):
    """
    Run CD-HIT to cluster and filter sequences based on identity.

    Parameters:
    input_file (str): Path to input FASTA file.
    output_file (str): Path to output FASTA file after clustering.
    identity_cutoff (float): Sequence identity cutoff for clustering.
    """
    subprocess.run(["cd-hit", "-i", input_file, "-o", output_file, "-c", str(identity_cutoff)])


# 4. Run BLAST for Taxonomic Classification
def run_blast(input_file, output_file, database="nt", evalue=0.001, num_threads=4):
    """
    Run BLAST for taxonomic classification of viral sequences.

    Parameters:
    input_file (str): Path to input FASTA file.
    output_file (str): Path to output file for BLAST results.
    database (str): BLAST database to use.
    evalue (float): E-value threshold for BLAST.
    num_threads (int): Number of threads to use for BLAST.
    """
    subprocess.run(
        ["blastn", "-query", input_file, "-db", database, "-out", output_file, "-evalue", str(evalue), "-num_threads",
         str(num_threads)])


# 5. Create Phylogenetic Tree (Multiple Options)
def create_phylogenetic_tree_mafft(input_file, alignment_file, tree_output_file):
    """
    Create a phylogenetic tree using MAFFT for alignment and FastTree for tree construction.

    Parameters:
    input_file (str): Path to input FASTA file.
    alignment_file (str): Path to output aligned sequences file.
    tree_output_file (str): Path to output Newick tree file.
    """
    subprocess.run(["mafft", input_file, "-o", alignment_file])
    subprocess.run(["FastTree", alignment_file, "-out", tree_output_file])


def create_phylogenetic_tree_muscle(input_file, alignment_file, tree_output_file):
    """
    Create a phylogenetic tree using MUSCLE for alignment and FastTree for tree construction.

    Parameters:
    input_file (str): Path to input FASTA file.
    alignment_file (str): Path to output aligned sequences file.
    tree_output_file (str): Path to output Newick tree file.
    """
    subprocess.run(["muscle", "-in", input_file, "-out", alignment_file])
    subprocess.run(["FastTree", alignment_file, "-out", tree_output_file])


# 6. Visualize Phylogenetic Tree
def visualize_tree(tree_file):
    """
    Visualize a phylogenetic tree using Bio.Phylo.

    Parameters:
    tree_file (str): Path to Newick tree file.
    """
    tree = Phylo.read(tree_file, "newick")
    Phylo.draw(tree)


# 7. Data Analysis and Visualization
def plot_sequence_length_distribution(metadata_file):
    """
    Plot the sequence length distribution from a metadata CSV file.

    Parameters:
    metadata_file (str): Path to the CSV file containing sequence metadata.
    """
    df = pd.read_csv(metadata_file)
    plt.hist(df['sequence_length'], bins=50)
    plt.xlabel('Sequence Length')
    plt.ylabel('Frequency')
    plt.title('Sequence Length Distribution')
    plt.show()


# 8. Functional Annotation Step
def run_functional_annotation(input_file):
    """
    Run functional annotation using InterProScan or similar tools.

    Parameters:
    input_file (str): Path to input FASTA file.
    """
    subprocess.run(["interproscan.sh", "-i", input_file, "-o", "functional_annotation.out", "-f", "tsv"])


# Pipeline Execution

# Set parameters
search_terms = ["unclassified virus", "unknown virus", "environmental virus"]
email = "your_email@example.com"
retmax = 100

# Fetch sequences
sequences = fetch_viral_sequences(search_terms, email, retmax)

# Save sequences to file
save_sequences_to_file(sequences, "viral_sequences.fasta")

# Run CD-HIT for clustering
run_cd_hit("viral_sequences.fasta", "filtered_sequences.fasta")

# Run BLAST for classification
run_blast("filtered_sequences.fasta", "blast_results.out")

# Create a phylogenetic tree using MAFFT and FastTree
create_phylogenetic_tree_mafft("filtered_sequences.fasta", "aligned_sequences_mafft.fasta", "tree_file_mafft.nwk")

# Alternatively, create a phylogenetic tree using MUSCLE and FastTree
# create_phylogenetic_tree_muscle("filtered_sequences.fasta", "aligned_sequences_muscle.fasta", "tree_file_muscle.nwk")

# Visualize the tree
visualize_tree("tree_file_mafft.nwk")

# Assuming metadata file exists (can be generated manually based on sequence data)
metadata_file = "sequence_metadata.csv"
plot_sequence_length_distribution(metadata_file)

# Run functional annotation
run_functional_annotation("filtered_sequences.fasta")

# End of the pipeline
