import os


def read_fastq(input_fastq: str):

    current_directory = os.getcwd()
    file_path = os.path.join(current_directory, input_fastq)
    
    with open(file_path, "r") as f:
        while True:

            header = f.readline().strip()
            if not header:
                break
            sequence = f.readline().strip()
            f.readline()
            quality = f.readline().strip()

            yield header, (sequence, quality)


def write_filtered_fastq(output_fastq: str, filtered_seqs):

    if os.path.exists(output_fastq):
        raise FileExistsError(
            f"File {output_fastq} exists. Please, choose another filename."
        )

    with open(output_fastq, "w") as f:
        for header, (sequence, quality) in filtered_seqs.items():
            f.write(f"{header}\n")
            f.write(f"{sequence}\n")
            f.write("+\n")
            f.write(f"{quality}\n")


def gc_filter(seq: str, gc_bounds: tuple):

    if isinstance(gc_bounds, (int, float)):
        gc_bounds = (0, gc_bounds)

    gc_content = (seq.count("C") + seq.count("G")) / len(seq) * 100

    return gc_bounds[0] <= gc_content <= gc_bounds[1]


def length_filter(seq: str, length_bounds: tuple):

    if isinstance(length_bounds, (int, float)):
        length_bounds = (
            0,
            length_bounds,
        )
    return length_bounds[0] <= len(seq) <= length_bounds[1]


def quality_filter(quality_string: str, threshold: int):
    unique_values = set(quality_string)
    total_quality = 0

    for i in unique_values:
        total_quality += quality_string.count(i) * (ord(i) - 33)

    total_quality = total_quality / len(quality_string)
    return total_quality >= threshold
