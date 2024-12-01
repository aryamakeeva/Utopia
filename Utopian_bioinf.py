from modules.modules_filter_fastq import (
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
    seqs: dict,
    gc_bounds: tuple = (0, 100),
    length_bounds: tuple = (0, 2 ** 32),
    quality_threshold: int = 0,
):
    """
    Filters FASTQ sequences in a dict using
    the connected validation functions.

    seqs: dict - a dictionary {name: (sequence, quality)}.
    gc_bounds: tuple - bounds for GC content (%).
    length_bounds: tuple - bounds for sequence length.
    quality_threshold: int - minimum quality value.
    """

    # Чтобы не перебирать словарь при вызове каждой проверочной функции,
    # обратимся к нему прямо здесь,
    # а функциям будем передавать конкретные строки из него
    # Новый словарь для отфильтрованных последовательностей
    filtered_seqs = {}

    # Перебираем элементы словаря и применяем фильтры
    for key, value in seqs.items():
        # Применяем все фильтры
        # если все фильтры проходят, добавляем в новый словарь
        if (
            gc_filter(value[0], gc_bounds) != "EXCLUDE"
            and length_filter(value[0], length_bounds) != "EXCLUDE"
            and quality_filter(value[1], quality_threshold) != "EXCLUDE"
        ):
            filtered_seqs[key] = value  # Добавляем в новый словарь

    return filtered_seqs


def run_dna_rna_tools(*args):
    """
    Function for base molecular operations with DNA/RNA sequences

    :args:  accepts one or more n.a. sequences and performs 4
    operations

    :operations:
        "transcribe"- transcription of DNA into RNA
        "reverse" - returns the reverse sequence
        "complement" - returns the complementary sequence
        "reverse_complement" - returns the reverse complementary sequence
        The function excludes the input sequence for the simultaneous
        presence of "U", "u"  and "T", "t" in the composition
        Function takes into account the register of nucleotides in the sequence

    The funcion uses additional module dna_rna_tools
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
