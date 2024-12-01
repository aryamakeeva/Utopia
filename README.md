# Utopian bioinf
Learning how to code softly and gently

## Description
It is a Python toolkit designed for DNA/RNA sequence manipulations and FASTQ sequence filtering based on GC content, sequence length and quality threshold. The toolkit provides two main functionalities:

1. `run_dna_rna_tools` : Manipulate DNA/RNA sequences (generate transcribed, reversed, complemented and  reverse-complemented sequences).
2. `filter_fastq`: Filter FASTQ sequences from dictionary (for now) based on defined thresholds for GC content, length, and quality.

## Functions Overview
1. `run_dna_rna_tools`

This function allows users to perform DNA/RNA sequence manipulations, including:

- transcribe: Convert DNA to RNA.
- reverse: Reverse a sequence.
- complement: Generate the complement of a sequence.
- reverse_complement: Generate the reverse complement.

Arguments: *args: The sequences followed by the procedure type ("transcribe", "reverse", "complement", "reverse_complement").

Returns: A single sequence (if one is provided) or a list of sequences.

2. `filter_fastq`

This function filters fastq sequences based on the following criteria: GC content, sequence length and quality.

Arguments:

- seq: dict - a dictionary {name: (sequence, quality)}
- output_fastq: filtered sequences (saved in the filtered folder).
- gc_bounds: The GC content bounds for filtering. Default is (0, 100).
- length_bounds: The sequence length bounds for filtering. Default is (0, 2**32).
- quality_threshold: The threshold for filtering by sequence quality. Default is 0.

Returns: filtered_seqs: filtered dictionary {name: (sequence, quality)}
