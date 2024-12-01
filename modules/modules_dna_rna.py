def validation(*seqs):
    dna_bases = {"A", "a", "T", "t", "C", "c", "G", "g"}
    rna_bases = {"A", "a", "U", "u", "C", "c", "G", "g"}

    results = []
    for seq in seqs:
        sequence_set = set(seq)

        if sequence_set.issubset(dna_bases):
            results.append("DNA")
        elif sequence_set.issubset(rna_bases):
            results.append("RNA")
        else:
            results.append("NotNA")

    return results


def reverse(*seqs):
    return [seq[::-1] for seq in seqs]


def transcribe(*seqs):  # мы не транскрибируем  из РНК ДНК
    transcribe_d = {
        "T": "U",
        "A": "A",
        "G": "G",
        "C": "C",
        "t": "u",
        "a": "a",
        "g": "g",
        "c": "c",
    }

    out_seqs = []

    for seq in seqs:
        out_seq = []
        for nucl in seq:
            out_seq.append(transcribe_d.get(nucl, ""))
        out_seqs.append("".join(out_seq))

    return out_seqs[0] if len(out_seqs) == 1 else out_seqs


def complement(*seqs):
    complement = {
        "A": "T",
        "C": "G",
        "G": "C",
        "T": "A",
        "a": "t",
        "c": "g",
        "g": "c",
        "t": "a",
        "U": "A",
        "u": "a",
    }

    out_seqs = []

    for seq in seqs:
        out_seq = []
        for nucl in seq:
            if nucl in complement:
                out_seq.append(complement[nucl])

        out_seqs.append("".join(out_seq))

    return out_seqs[0] if len(out_seqs) == 1 else out_seqs


def reverse_complement(*seqs):
    return complement(*reverse(*seqs))
