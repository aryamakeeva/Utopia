"""
Module for DNA and RNA sequence manipulations, including validation, transcription, 
and generation of complement and reverse complement sequences.
"""

def validation(*seqs: str) -> list[str]:
    """
    Validates DNA or RNA sequences.
    """
    dna_bases = {"A", "a", "T", "t", "C", "c", "G", "g"}
    rna_bases = {"A", "a", "U", "u", "C", "c", "G", "g"}

    return [
        "DNA" if set(seq).issubset(dna_bases) else "RNA" if set(seq).issubset(rna_bases) else "NotNA"
        for seq in seqs
    ]

def reverse(*seqs: str) -> list[str]:
    """
    Reverses the input DNA/RNA sequences.
    """
    return [seq[::-1] for seq in seqs]


def transcribe(*seqs: str) -> list[str]: 
    """
    Transcribes DNA sequences to RNA by replacing 'T' with 'U' and 't' with 'u'.
    """
    return [seq.replace("T", "U").replace("t", "u") for seq in seqs]


def complement(*seqs: str) -> list[str]:
    """
    Returns the complementary DNA/RNA sequence for each input sequence.
    """
    complement_dict = {
        "A": "T", "C": "G", "G": "C", "T": "A",
        "a": "t", "c": "g", "g": "c", "t": "a",
        "U": "A", "u": "a"
    }

    return [
        "".join([complement_dict.get(nucl, "") for nucl in seq]) for seq in seqs
    ]


def reverse_complement(*seqs: str) -> list[str]:
    """
    Returns the reverse complement of the input DNA/RNA sequences.
    """
    return complement(*reverse(*seqs))
