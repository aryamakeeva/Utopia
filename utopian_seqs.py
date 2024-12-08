from modules.modules_filter_fastq import (
    read_fastq,
    write_filtered_fastq,
    gc_filter,
    quality_filter,
    length_filter,
)

from modules.modules_dna_rna import (
    reverse,
    transcribe,
    reverse_complement,
    complement,
    validation,
)


def filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds: tuple = (0, 100),
    length_bounds: tuple = (0, 2 ** 32),
    quality_threshold: int = 0,
):
    """
    Filters FASTQ sequences by GC content,
    sequence length, and quality threshhold.
    -------------------------------
    Parameters:

    input_fastq : str
        - pass to file from current directory or
          name of the FASTQ file in current directory
            

    output_fastq : str
        - name of the filtered FASTQ file.

    gc_bounds: tuple
        - GC content interval (%) for filtering.
        If a single number is provided, it is considered the upper bound,
        and the lower bound is set to 0.
        Defaults to (0, 100).

    length_bounds: tuple
        - sequence length interval for filtering.
        Defaults to (0, 2**32).

    quality_threshold: int
        - quality threshold for filtering.
        Defaults to 0.
    -------------------------------
    ValueErrors:
    If the name of output file already exists.
    """

    filtered_seqs = {}

    for header, (sequence, quality) in read_fastq(input_fastq):
        if (
            gc_filter(sequence, gc_bounds)
            and length_filter(sequence, length_bounds)
            and quality_filter(quality, quality_threshold)
        ):
            filtered_seqs[header] = (sequence, quality)

    if filtered_seqs:  # Проверяем, есть ли отфильтрованные данные
        write_filtered_fastq(output_fastq, filtered_seqs)
    else:
        print("The reads were filtered out... Survival rate: 0%.")


def run_dna_rna_tools(*args):
    """
    Performs base molecular operations with DNA/RNA sequences
    -------------------------------
    Parameters:

    *args:  accepts one or more n.a. sequences and performs 4
    operations, where the last argument is the function name

    operations:
    "transcribe"- transcription of DNA into RNA

    "reverse" - returns the reverse sequence

    "complement" - returns the complementary sequence

    "reverse_complement" - returns the reverse complementary sequence

    The function excludes the input sequence for the simultaneous
    presence of "U", "u"  and "T", "t" in the composition

    Function takes into account the register of nucleotides in the sequence
    -------------------------------
    ValueErrors
    - If sequence is not DNA or RNA
    - If user tries to transcribe RNA-sequence
    - If molucular operation is not listed in description

    """
    *seqs, action = args
    seq_types = validation(*seqs)

    for seq, seq_type in zip(seqs, seq_types):
        if seq_type == "NotNA":
            raise ValueError(
                "Invalid sequence: \
                the sequence must contain only DNA or RNA nucleotides"
            )

        if action == "transcribe" and seq_type != "DNA":
            raise ValueError(
                "Invalid sequence: \
            a DNA sequence is required for transcription"
            )

    if action == "transcribe":
        dna_seqs = [seq for seq, t in zip(seqs, seq_types) if t == "DNA"]
        result = transcribe(*dna_seqs)
    elif action == "reverse":
        reversed_seqs = reverse(*seqs)
        result = reversed_seqs[0] if len(reversed_seqs) == 1 else reversed_seqs
    elif action == "complement":
        result = complement(*seqs)
    elif action == "reverse_complement":
        result = reverse_complement(*seqs)
    else:
        raise ValueError("Action is not in the list of options")

    return result
