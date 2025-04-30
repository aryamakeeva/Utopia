from abc import ABC, abstractmethod


class BiologicalSequence(ABC):
    def __init__(self, seq: str):
        self.seq = seq
        self._validate_alphabet()

    def __len__(self) -> int:
        return len(self.seq)

    def __getitem__(self, index: int) -> str:
        return self.seq[index]

    def __repr__(self) -> str:
        return self.seq

    def __str__(self):
        return self.seq

    def __add__(self, other: "BiologicalSequence") -> "BiologicalSequence":
        """
        Concatenates two sequences.
        :param other: Another BiologicalSequence instance.
        :return: A new instance with concatenated sequences.
        """
        return self.__class__(self.seq + other.seq)

    @abstractmethod
    def _validate_alphabet(self) -> None:
        pass


class NucleicAcidSequence(BiologicalSequence):
    complements: dict[str, str] = {}
    valid_bases: set[str] = set()

    def __init__(self, seq: str):
        if type(self) is NucleicAcidSequence:
            raise NotImplementedError("Is it DNA or RNA?")
        super().__init__(seq)

    def _validate_alphabet(self) -> None:
        if not set(self.seq).issubset(self.valid_bases):
            raise ValueError(
                "Invalid sequence: must contain DNA or RNA nucleotides only"
            )

    def reverse(self) -> "NucleicAcidSequence":
        return self.__class__(self.seq[::-1])

    def complement(self) -> "NucleicAcidSequence":
        complement_seq = self.seq.translate(self.complements)
        return self.__class__(complement_seq)

    def reverse_complement(self) -> "NucleicAcidSequence":
        return self.complement().reverse()


class DNASequence(NucleicAcidSequence):

    complements = str.maketrans("ATGCatgc", "TACGtacg")
    valid_bases = set("ATGCatgc")

    def transcribe(self) -> "RNASequence":
        return RNASequence(self.seq.replace("T", "U").replace("t", "u"))


class RNASequence(NucleicAcidSequence):

    complements = str.maketrans("AUGCaugc", "UACGuacg")
    valid_bases = set("AUGCaugc")


class AminoAcidSequence(BiologicalSequence):

    valid_bases = set("ACDEFGHIKLMNPQRSTVWY")

    def _validate_alphabet(self):
        if not set(self.seq).issubset(self.valid_bases):
            raise ValueError("Invalid sequence: must contain only amino-acids")

    def transcript_len(self) -> int:
        return len(self.seq) * 3
