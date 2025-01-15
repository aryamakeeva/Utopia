from modules.modules_dna_rna import (
    complement,
    reverse,
    reverse_complement,
    transcribe,
    validation,
)

from modules.modules_fastq import (
    gc_filter,
    length_filter,
    quality_filter,
    read_fastq,
    write_filtered_fastq,
)


def fastq_filter_tool(
    input_fastq: str,
    output_fastq: str,
    gc_bounds: tuple[float, float] = (0, 100),
    length_bounds: tuple[int, int] = (0, 2 ** 32),
    quality_threshold: int = 0,
):
    """
    Filters FASTQ sequences by GC content,
    sequence length, and quality threshhold.
    -------------------------------
    Parameters:

    :param input_fastq: str
    - pass to file from current directory or
      name of the FASTQ file in current directory
    :type input_fastq: str

    :param output_fastq: str
    - name of the filtered FASTQ file.
    :type output_fastq: str

    :param gc_bounds: tuple
    - GC content interval (%) for filtering.
      If a single number is provided, it is considered the upper bound,
      and the lower bound is set to 0.
      Defaults to (0, 100).
    :type gc_bounds: tuple[float, float]

    :param length_bounds: tuple
    - sequence length interval for filtering.
      Defaults to (0, 2**32).
    :type length_bounds: tuple[int, int]

    :param quality_threshold: int
    - quality threshold for filtering.
      Defaults to 0.
    :type quality_threshold: int
    -------------------------------
    ValueErrors:
    If the name of output file already exists.
    """
    filtered_seqs = {}

    for header, (sequence, comment, quality) in read_fastq(input_fastq):
        if (
            gc_filter(sequence, gc_bounds)
            and length_filter(sequence, length_bounds)
            and quality_filter(quality, quality_threshold)
        ):
            filtered_seqs[header] = (sequence, comment, quality)

    if filtered_seqs:  
        write_filtered_fastq(output_fastq, filtered_seqs)
    else:
        print("The reads were filtered out... Survival rate: 0%.")


def dna_rna_tool (*args: str) -> str | list[str]:
    """
    Performs base molecular operations with DNA/RNA sequences
    -------------------------------
    Parameters:

    :param *args: 
    - accepts one or more nucleotide sequences and performs 4 operations, 
      where the last argument is the function name.
    :type *args: str

    Operations:
    "transcribe" - transcription of DNA into RNA
    "reverse" - returns the reverse sequence
    "complement" - returns the complementary sequence
    "reverse_complement" - returns the reverse complementary sequence

    The function excludes the input sequence for the simultaneous presence of "U", "u" 
    and "T", "t" in the composition.

    Function takes into account the register of nucleotides in the sequence.
    -------------------------------
    ValueErrors:
    - If sequence is not DNA or RNA
    - If user tries to transcribe an RNA sequence
    - If molecular operation is not listed in description
    """
    *seqs, action = args
    seq_types = validation(*seqs)

    for seq, seq_type in zip(seqs, seq_types):
        if seq_type == "NotNA":
            raise ValueError(
                "Invalid sequence: the sequence must contain only DNA or RNA nucleotides"
            )

        if action == "transcribe" and seq_type != "DNA":
            raise ValueError(
                "Invalid sequence: a DNA sequence is required for transcription"
            )

    action_map = {
        "transcribe": transcribe,
        "reverse": reverse,
        "complement": complement,
        "reverse_complement": reverse_complement,
    }

    if action not in action_map:
        raise ValueError("Action is not in the list of options")

    func = action_map[action]

    if action == "transcribe":
        dna_seqs = [seq for seq, t in zip(seqs, seq_types) if t == "DNA"]
        result = func(*dna_seqs)
    else:
        result = func(*seqs)

    if len(result) == 1:
        return result[0]
    else:
        return ", ".join(result)
