# Viral Sequence Analysis Pipeline

## Overview
This project provides a Python-based pipeline for analyzing viral sequences retrieved from GenBank. The pipeline includes steps for sequence retrieval, quality control, taxonomic classification, phylogenetic analysis, data visualization, and functional annotation. It is designed to streamline the analysis of cryptic or unclassified viral sequences.

## Features
- **Fetch Viral Sequences from GenBank**: Retrieve viral sequences using specific keywords.
- **Quality Control and Clustering**: Use CD-HIT to cluster sequences and filter out redundancies.
- **Taxonomic Classification**: Run BLAST to determine the potential classification of viral sequences.
- **Phylogenetic Analysis**: Create phylogenetic trees using MAFFT or MUSCLE for alignment and FastTree for tree construction.
- **Visualization**: Visualize phylogenetic trees and sequence length distributions.
- **Functional Annotation**: Annotate sequences using InterProScan or similar tools.

## Prerequisites
To run this pipeline, you need the following software and Python libraries:
- Python 3.x
- [Biopython](https://biopython.org/)
- [Pandas](https://pandas.pydata.org/)
- [Matplotlib](https://matplotlib.org/)
- CD-HIT
- BLAST
- MAFFT or MUSCLE
- FastTree
- InterProScan (optional for functional annotation)

## Installation
1. Install Python and necessary Python libraries:
   ```bash
   pip install biopython pandas matplotlib
   ```
2. Install bioinformatics tools:
   - [CD-HIT](http://weizhongli-lab.org/cd-hit/)
   - [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
   - [MAFFT](https://mafft.cbrc.jp/alignment/software/) or [MUSCLE](https://www.drive5.com/muscle/)
   - [FastTree](http://www.microbesonline.org/fasttree/)
   - [InterProScan](https://www.ebi.ac.uk/interpro/interproscan.html) (optional)

## Usage
1. **Set Parameters**: Edit the parameters such as search terms, email, and maximum number of records to retrieve.
   ```python
   search_terms = ["unclassified virus", "unknown virus", "environmental virus"]
   email = "your_email@example.com"
   retmax = 100
   ```

2. **Run the Pipeline**: Execute the script to perform the following tasks:
   - Fetch sequences from GenBank.
   - Save sequences to a file.
   - Run CD-HIT for quality control.
   - Classify sequences using BLAST.
   - Create a phylogenetic tree using MAFFT or MUSCLE and FastTree.
   - Visualize the phylogenetic tree and sequence length distribution.
   - Run functional annotation (optional).

3. **Example Command**:
   ```bash
   python viral_sequence_analysis_pipeline.py
   ```

## Pipeline Steps
1. **Fetch Viral Sequences**: Retrieve sequences from GenBank using keywords like "unclassified virus."
2. **Save Sequences to File**: Save retrieved sequences to a FASTA file for downstream analysis.
3. **Quality Control**: Use CD-HIT to cluster sequences and remove redundancies.
4. **Taxonomic Classification**: Run BLAST against the nucleotide database to determine the classification of the sequences.
5. **Phylogenetic Analysis**: Align sequences using MAFFT or MUSCLE and generate a phylogenetic tree with FastTree.
6. **Visualization**: Visualize the phylogenetic tree and plot sequence length distribution using Matplotlib.
7. **Functional Annotation**: Optionally, run functional annotation using InterProScan to predict protein functions.

## Example Output
- **Phylogenetic Tree**: A Newick format tree file that can be visualized.
- **BLAST Results**: A file containing BLAST hits and taxonomic classification.
- **Sequence Length Distribution**: A histogram showing the distribution of sequence lengths.

## Notes
- Ensure all required bioinformatics tools are installed and accessible in your system PATH.
- Some steps, like BLAST and functional annotation, may require significant computational resources, depending on the number of sequences.

## License
This project is licensed under the MIT License - see the LICENSE file for details.

## Contact
For any questions or issues, please contact [your_email@example.com].

