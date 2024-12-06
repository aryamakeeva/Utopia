# Utopian tools
Learning how to code softly and gently. This repository is the result of half-year Python course. It provite 2 toolskits: `Utopian Seqs` and `Utopian Files`.

## Installation
```bash
git clone git@github.com:aryamakeeva/Utopia.git
```
## Utopian Seqs

### Description
Is a Python toolkit designed for operations with DNA/RNA sequences and FASTQ sequence filtration based on GC content, sequence length and quality threshold. The toolkit provides two main functionalities:

1. `run_dna_rna_tools` : Manipulates with DNA/RNA sequences, including functions to generate transcribed, reversed, complemented, and reverse-complemented sequences.
2. `filter_fastq`: Filter FASTQ sequences based on defined thresholds for GC content, length, and quality.


### `run_dna_rna_tools`

This function allows users to perform DNA/RNA sequence operations, including:

- transcribe: Convert DNA to RNA.
- reverse: Reverse a sequence.
- complement: Generate the complement of a sequence.
- reverse_complement: Generate the reverse complement.

Input:

One or more n.a. sequences followed by one of 4 operations, where the last argument is the function name.

Returns: A single sequence (if one is provided) or a list of sequences.

~Example~

```python
from bioinf_seqs import run_dna_rna_tools

# For a single sequence
result = run_dna_rna_tools("ATG", "transcribe")
# Returns: "AUG"

# For multiple sequences
results = run_dna_rna_tools("ttG", "AT", "ATc", "complement")
# Returns: "aaC", "TA", "TAg" 
```


### `filter_fastq`

This function filters Filters FASTQ sequences by GC content, sequence length, and quality threshhold.

Input:

- input_fastq : pass to file from current directory or name of the FASTQ file in current directory
- output_fastq : name of the filtered FASTQ file.
- gc_bounds: GC content interval (%) for filtering. Defaults to (0, 100).

  If a single number is provided, it is considered the upper bound, and the lower bound is set to 0.

- length_bounds: sequence length interval for filtering. Defaults to (0, 2**32).


- quality_threshold: quality threshold for filtering. Defaults to 0.

Example

```python
from utopian_seqs import filter_fastq

# Input
filter_fastq(
    input_fastq="pass_to_file/file_name.fastq",
    output_fastq="output_file_name.fastq",
    gc_bounds=(30, 80),
    length_bounds=(10, 50),
    quality_threshold=20
)
```
## Utopian Files

### Description

This Python toolkit is designed for processing bioinformatics files, offering functions to manipulate and extract data from various bioinformatics file formats. The toolkit includes the following functionalities:

1. `convert_multiline_fasta_to_oneline`: Takes an input FASTA file where sequences (DNA/RNA/protein, etc.) are splited across multiple lines and makes a new FASTA file where each sequence is contained in a single line. 
  
2. `parse_blast_output`: Extracts the best matches from BLAST results for each sequence and saves the names of the best matching proteins in a new file, sorted alphabetically. 

3. `select_genes_from_gbk_to_fasta`: This function extracts specific genes, their aminoacid sequence and surrounding genes from a GenBank file. The extracted protein sequences (translations) are saved in a FASTA file, which can be directly used in BLAST analysis for further exploration.

### `convert_multiline_fasta_to_oneline`

- input_fasta: the path to the input FASTA file containing sequences split across multiple lines.
    
- output_fasta: the path to the output FASTA file where the single-line sequences will be saved.

If not provided, the output will be saved to a default file.

### `parse_blast_output`

This function parses a BLAST output file and extracts the first description line for each QUERY, saves the list of proteins to a new file sorted alphabetically.

- input_file: the path to the input BLAST output file (txt format).
        
- output_file: the path to the output file where the extracted protein descriptions will be saved.

```python
from utopian_files import parse_blast_output

# Input
parse_blast_output('./filename.txt', 'output_filename.txt')
```

### `select_genes_from_gbk_to_fasta`

This function extracts genes and their neighboring genes from a GBK file, retrieves their amino acid sequences, and saves them to a FASTA file.

```python
from utopian_files import select_genes_from_gbk_to_fasta

# Input
select_genes_from_gbk_to_fasta(
    input_gbk="./example/example_gbk.gbk",
    genes=["dtp"],
    n_before = 1, #default
    n_after = 1, #default 
    output_fasta="output.fasta"
)

#Output

Matching genes:
1: dtpD
2: dtpC
3: dtpB
4: dtpA
Enter the numbers of desired genes separated by commas (e.g., 1,3):  1
Gene(s) and neighbours are now in output.fasta.
```


## Installation
To use the Utopian_bioinf tool, follow steps bellow:
```bash
git clone git@github.com:aryamakeeva/Utopia.git
cd Utopia
python Utopian_bioinf.py
```

## Contact

You can send your feedback to aryamakeeva@gmail.com
