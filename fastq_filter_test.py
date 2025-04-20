import os
import pytest
from Bio import SeqIO
from utopian_fastqfilter import fastq_filter_tool


@pytest.fixture
def tmp_fastq_files():
    fastq_data = """@seq1
TAAATATCG
+
IIIIIIIII
@seq2
ATGCATGCATGC
+
IIIIIIIIIIII
@seq3
GCGC
+
IIII
"""
    input_path = "tmp_input.fastq"
    output_path = "tmp_output.fastq"
    with open(input_path, "w") as f:
        f.write(fastq_data)
    yield input_path, output_path
    for path in [input_path, output_path]:
        if os.path.exists(path):
            os.remove(path)


@pytest.fixture
def wrong_fasta_file():
    wrong_data = """>seq1
TAAAFATCG
>seq2
ATGCATGCATGC
>seq3
GCGC
"""
    wrong_input_path = "wrong_input.fastq"
    with open(wrong_input_path, "w") as f:
        f.write(wrong_data)
    yield wrong_input_path
    if os.path.exists(wrong_input_path):
        os.remove(wrong_input_path)


class TestGCFilter:
    def test_gc_lower_bound(self, tmp_fastq_files):
        input_fastq, output_fastq = tmp_fastq_files
        fastq_filter_tool(
            input_fastq=input_fastq,
            output_fastq=output_fastq,
            gc_bounds=(80, 100),
            length_bounds=(0, 100),
            quality_threshold=0,
        )
        with open(output_fastq, "r") as f:
            records = list(SeqIO.parse(f, "fastq"))
            assert len(records) == 1

    def test_gc_upper_bound(self, tmp_fastq_files):
        input_fastq, output_fastq = tmp_fastq_files
        fastq_filter_tool(
            input_fastq=input_fastq,
            output_fastq=output_fastq,
            gc_bounds=(0, 30),
            length_bounds=(0, 100),
            quality_threshold=0,
        )
        with open(output_fastq, "r") as f:
            records = list(SeqIO.parse(f, "fastq"))
            assert len(records) == 1


class TestLengthFilter:
    def test_length_lower_bound(self, tmp_fastq_files):
        input_fastq, output_fastq = tmp_fastq_files
        fastq_filter_tool(
            input_fastq=input_fastq,
            output_fastq=output_fastq,
            gc_bounds=(0, 100),
            length_bounds=(10, 100),
            quality_threshold=0,
        )
        with open(output_fastq, "r") as f:
            records = list(SeqIO.parse(f, "fastq"))
            assert len(records) == 1

    def test_length_upper_bound(self, tmp_fastq_files):
        input_fastq, output_fastq = tmp_fastq_files
        fastq_filter_tool(
            input_fastq=input_fastq,
            output_fastq=output_fastq,
            gc_bounds=(0, 100),
            length_bounds=(0, 5),
            quality_threshold=0,
        )
        with open(output_fastq, "r") as f:
            records = list(SeqIO.parse(f, "fastq"))
            assert len(records) == 1


class TestQualityFilter:
    def test_quality_threshold(self, tmp_fastq_files):
        input_fastq, output_fastq = tmp_fastq_files
        fastq_filter_tool(
            input_fastq=input_fastq,
            output_fastq=output_fastq,
            gc_bounds=(0, 100),
            length_bounds=(0, 100),
            quality_threshold=41,
        )
        with open(output_fastq, "r") as f:
            records = list(SeqIO.parse(f, "fastq"))
            assert len(records) == 0


class TestErrorHandling:
    def test_invalid_file(self):
        with pytest.raises(FileNotFoundError):
            fastq_filter_tool(
                input_fastq="non_existing.fastq",
                output_fastq="output.fastq",
                gc_bounds=(0, 100),
                length_bounds=(0, 100),
                quality_threshold=0,
            )

    def test_invalid_format(self, wrong_fasta_file):
        output_fastq = "tmp_output_invalid.fastq"
        try:
            with pytest.raises(ValueError):
                fastq_filter_tool(
                    input_fastq=wrong_fasta_file,
                    output_fastq=output_fastq,
                    gc_bounds=(0, 100),
                    length_bounds=(0, 100),
                    quality_threshold=0,
                )
        finally:
            if os.path.exists(output_fastq):
                os.remove(output_fastq)


class TestIO:
    def test_output_creation(self, tmp_fastq_files):
        input_fastq, output_fastq = tmp_fastq_files
        fastq_filter_tool(
            input_fastq=input_fastq,
            output_fastq=output_fastq,
            gc_bounds=(0, 100),
            length_bounds=(0, 100),
            quality_threshold=0,
        )
        assert os.path.exists(output_fastq)
