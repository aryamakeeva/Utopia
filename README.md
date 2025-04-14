# Utopian tools
Learning how to code softly and gently. 

This repository is the result of a half-year Python course and provides two toolkits: `Utopian Seqs` and `Utopian Files`. 

- `Utopian Seqs` is designed for operations with DNA/RNA/AminoAcid sequences and FASTQ files filtration based on GC content, sequence length, and quality threshold.  
- `Utopian Files` is designed for processing bioinformatics files:
  - Converts a multi-line FASTA file into a new format where each sequence is on a single line.
  - Extracts and alphabetically sorts the best matches from BLAST results, saving the names of top matching proteins.
  - Extracts specific genes and their protein sequences from GenBank files, saving them in FASTA format for use in BLAST analysis.

## Content
- [Installation](#installation)
- [Usage](#usage)
  - [Utopian Seqs](#utopian-seqs)
    - [dna_rna_tool](#dna_rna_tool)
    - [fastq_filter_tool](#fastq_filter_tool)
  - [Utopian Files](#utopian-files)
    - [convert_multiline_fasta_to_oneline](#convert_multiline_fasta_to_oneline)
    - [parse_blast_output](#parse_blast_output)
    - [select_genes_from_gbk_to_fasta](#select_genes_from_gbk_to_fasta)
- [Contact](#contact)

## Installation 
```bash
git clone git@github.com:aryamakeeva/Utopia.git
```
## Usage

### Utopian Seqs
The toolkit provides several classes and functionalities for manipulating biological sequences:

1. `BiologicalSequence`
This is the **base class** for all biological sequences. It includes common functionality such as:

- **Length Calculation**
- **Indexing**
- **Concatenation**
- **String Representation**

2. `NucleicAcidSequence` (Subclass of `BiologicalSequence`)
This class is for **DNA/RNA sequences**. It includes methods that:

- Validate DNA or RNA sequences
- Generate reversed, complemented, and reverse-complemented sequences

**Subclasses**:

- **`DNASequence`**: Works with DNA sequences with additional functionality for **transcription** to RNA.
- **`RNASequence`**

3. `AminoAcidSequence` (Subclass of `BiologicalSequence`)
This class is for **amino acid sequences**. It includes methods that:

- Validate amino acid sequences
- Calculate the **transcript length** 
  
*Input:*

One or more n.a. sequences followed by one of 4 operations, where the last argument is the function name.

Returns: A single sequence (if one is provided) or a list of sequences.

*Example*

```python
from utopian_seqs import DNASequence, RNASequence

# For a single DNA sequence
dna = DNASequence("ATGCT")
dna.complement()
# Returns: 'TACGA'
dna.reverse()
# Returns: 'TCGTA'
dna.transcribe()
# Returns: 'AUGCU'
DNASequence("ATGC") + DNASequence("AATC")
# Returns: 'ATGCAATC'
```

2. `fastq_filter_tool`: Filters FASTQ sequences based on defined thresholds for GC content, length, and quality threshhold.

*Input:*

- input_fastq : pass to file from current directory or name of the FASTQ file in current directory
- output_fastq : name of the filtered FASTQ file.
- gc_bounds: GC content interval (%) for filtering. Defaults to (0, 100).

  If a single number is provided, it is considered the upper bound, and the lower bound is set to 0.

- length_bounds: sequence length interval for filtering. Defaults to (0, 2**32).

- quality_threshold: quality threshold for filtering. Defaults to 0.

*Example*

```python
from utopian_seqs import fastq_filter_tool

# Input
fastq_filter_tool(
    input_fastq="pass_to_file/file_name.fastq",
    output_fastq="output_file_name.fastq",
    gc_bounds=(30, 80),
    length_bounds=(10, 50),
    quality_threshold=20
)
```
### Utopian Files

The toolkit includes the following functionalities:

1. `convert_multiline_fasta_to_oneline`: Takes an input FASTA file where sequences (DNA/RNA/protein, etc.) are splited across multiple lines and makes a new FASTA file where each sequence is contained in a single line.

- input_fasta: the path to the input FASTA file containing sequences split across multiple lines.
    
- output_fasta: the path to the output FASTA file where the single-line sequences will be saved.

If not provided, the output will be saved to a default file.
  
2. `parse_blast_output`: Extracts the best matches from BLAST results for each sequence and saves the names of the best matching proteins in a new file, sorted alphabetically.

This function parses a BLAST output file and extracts the first description line for each QUERY, saves the list of proteins to a new file sorted alphabetically.

- input_file: the path to the input BLAST output file (txt format).
        
- output_file: the path to the output file where the extracted protein descriptions will be saved.

```python
from utopian_files import parse_blast_output

# Input
parse_blast_output('./filename.txt', 'output_filename.txt')
```

3. `select_genes_from_gbk_to_fasta`: This function extracts specific genes, their aminoacid sequence and surrounding genes from a GenBank file. The extracted protein sequences (translations) are saved in a FASTA file, which can be directly used in BLAST analysis for further exploration.

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

## Contact

You can send your feedback to aryamakeeva@gmail.com
