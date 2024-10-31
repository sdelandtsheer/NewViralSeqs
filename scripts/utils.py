# utils.py - Utility Functions for Viral Sequence Analysis Pipeline

import subprocess
from Bio import Entrez, SeqIO, Phylo
import requests
import os


def fetch_viral_sequences(search_terms, email, retmax=100):
    """
    Fetch viral sequences from GenBank.

    Parameters:
    search_terms (list): List of terms to search in GenBank (e.g., ["unclassified virus", "unknown virus"]).
    email (str): Email for NCBI Entrez identification.
    retmax (int): Maximum number of records to retrieve per term (default: 100).

    Returns:
    list: List of SeqRecord objects representing viral sequences.
    """
    Entrez.email = email
    sequences = []
    for term in search_terms:
        try:
            handle = Entrez.esearch(db="nucleotide", term=term, retmax=retmax)
            record = Entrez.read(handle)
            ids = record['IdList']
            handle = Entrez.efetch(db="nucleotide", id=ids, rettype="gb", retmode="text")
            sequences.extend(list(SeqIO.parse(handle, "genbank")))
        except Exception as e:
            print(f"Error fetching sequences for term '{term}': {e}")
    return sequences


def save_sequences_to_file(sequences, output_file):
    """
    Save sequences to a FASTA file.

    Parameters:
    sequences (list): List of SeqRecord objects representing viral sequences.
    output_file (str): Path to the output FASTA file.
    """
    try:
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        SeqIO.write(sequences, output_file, "fasta")
    except Exception as e:
        print(f"Error saving sequences to file '{output_file}': {e}")


def run_cd_hit(input_file, output_file, identity_cutoff=0.9):
    """
    Run CD-HIT to cluster and filter sequences based on identity.

    Parameters:
    input_file (str): Path to input FASTA file.
    output_file (str): Path to output FASTA file after clustering.
    identity_cutoff (float): Sequence identity cutoff for clustering (default: 0.9).
    """
    try:
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        subprocess.run(["cd-hit", "-i", input_file, "-o", output_file, "-c", str(identity_cutoff)], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running CD-HIT: {e}")


def run_blast(input_file, output_file, database="nt", evalue=0.001, num_threads=4):
    """
    Run BLAST for taxonomic classification of viral sequences.

    Parameters:
    input_file (str): Path to input FASTA file.
    output_file (str): Path to output file for BLAST results.
    database (str): BLAST database to use (default: "nt").
    evalue (float): E-value threshold for BLAST (default: 0.001).
    num_threads (int): Number of threads to use for BLAST (default: 4).
    """
    try:
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
        with open(input_file, 'r') as f:
            fasta_data = f.read()
        params = {
            'CMD': 'Put',
            'DATABASE': database,
            'PROGRAM': 'blastn',
            'MEGABLAST': 'on',
            'QUERY': fasta_data
        }
        response = requests.post(url, data=params)
        if response.status_code == 200:
            request_id = response.text.split("RID = ")[1].split()[0]
            # Check status and retrieve results (simplified)
            params = {'CMD': 'Get', 'RID': request_id, 'FORMAT_TYPE': 'Text'}
            response = requests.get(url, params=params)
            if response.status_code == 200:
                with open(output_file, 'w') as f:
                    f.write(response.text)
            else:
                raise Exception("BLAST result retrieval failed")
        else:
            raise Exception("BLAST server request failed")
    except Exception as e:
        print(f"Error running BLAST: {e}")


def create_phylogenetic_tree_mafft(input_file, alignment_file, tree_output_file):
    """
    Create a phylogenetic tree using MAFFT for alignment and FastTree for tree construction.

    Parameters:
    input_file (str): Path to input FASTA file.
    alignment_file (str): Path to output aligned sequences file.
    tree_output_file (str): Path to output Newick tree file.
    """
    try:
        os.makedirs(os.path.dirname(alignment_file), exist_ok=True)
        os.makedirs(os.path.dirname(tree_output_file), exist_ok=True)
        url = "https://mafft.cbrc.jp/alignment/server/alignments/submit/"
        files = {'file': open(input_file, 'rb')}
        response = requests.post(url, files=files)
        if response.status_code == 200:
            with open(alignment_file, 'w') as f:
                f.write(response.text)
        else:
            raise Exception("MAFFT server request failed")
        subprocess.run(["FastTree", alignment_file, "-out", tree_output_file], check=True)
    except Exception as e:
        print(f"Error creating phylogenetic tree: {e}")


def create_phylogenetic_tree_muscle(input_file, alignment_file, tree_output_file):
    """
    Create a phylogenetic tree using MUSCLE for alignment and FastTree for tree construction.

    Parameters:
    input_file (str): Path to input FASTA file.
    alignment_file (str): Path to output aligned sequences file.
    tree_output_file (str): Path to output Newick tree file.
    """
    try:
        os.makedirs(os.path.dirname(alignment_file), exist_ok=True)
        os.makedirs(os.path.dirname(tree_output_file), exist_ok=True)
        subprocess.run(["muscle", "-in", input_file, "-out", alignment_file], check=True)
        subprocess.run(["FastTree", alignment_file, "-out", tree_output_file], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error creating phylogenetic tree with MUSCLE: {e}")


def visualize_tree(tree_file):
    """
    Visualize a phylogenetic tree using Bio.Phylo.

    Parameters:
    tree_file (str): Path to Newick tree file.
    """
    try:
        tree = Phylo.read(tree_file, "newick")
        Phylo.draw(tree)
    except Exception as e:
        print(f"Error visualizing tree '{tree_file}': {e}")


def plot_sequence_length_distribution(metadata_file):
    """
    Plot the sequence length distribution from a metadata CSV file.

    Parameters:
    metadata_file (str): Path to the CSV file containing sequence metadata.
    """
    try:
        df = pd.read_csv(metadata_file)
        plt.hist(df['sequence_length'], bins=50)
        plt.xlabel('Sequence Length')
        plt.ylabel('Frequency')
        plt.title('Sequence Length Distribution')
        plt.show()
    except Exception as e:
        print(f"Error plotting sequence length distribution from '{metadata_file}': {e}")


def run_functional_annotation(input_file):
    """
    Run functional annotation using InterProScan or similar tools.

    Parameters:
    input_file (str): Path to input FASTA file.
    """
    try:
        subprocess.run(["interproscan.sh", "-i", input_file, "-o", "functional_annotation.out", "-f", "tsv"],
                       check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running functional annotation: {e}")
