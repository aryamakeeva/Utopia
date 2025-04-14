import re


def convert_multiline_fasta_to_oneline(
    input_fasta: str, output_fasta: str = "output.fasta"
):
    """
    Converts a multi-line FASTA sequence to a single-line format and saves it to a new file.
    -------------------------------
    Parametrs:

    input_fasta: str
        - the path to the input FASTA file containing sequences split across multiple lines.

    output_fasta: str (optional)
        - the path to the output FASTA file where the single-line sequences will be saved.
    If not provided, the output will be saved to a default file.
    """

    with open(input_fasta, "r") as infile, open(output_fasta, "w") as outfile:
        sequence = ""
        header = None

        for line in infile:
            line = line.strip()
            if line.startswith(">"):

                if header:
                    outfile.write(f"{header}\n{sequence}\n")

                header = line
                sequence = ""
            else:
                sequence += line

        if header:
            outfile.write(f"{header}\n{sequence}\n")


def parse_blast_output(input_file: str, output_file: str):
    """
    Parses a BLAST output file and extracts the first description line for each QUERY,
    saves the list of proteins to a new file sorted alphabetically.
    -------------------------------
    Parametrs:

    input_file: str
        - the path to the input BLAST output file (txt format).

    output_file: str
        - the path to the output file where the extracted protein descriptions will be saved.
    """

    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        lines = infile.readlines()
        result_lines = []

        for i, line in enumerate(lines):
            if "Sequences producing significant alignments:" in line:
                # Получаем строку через 3 после текущей
                target_line = lines[i + 3].strip()
                result_lines.append(target_line)

        sorted_lines = sorted(result_lines)
        outfile.write("\n".join(sorted_lines) + "\n")


def select_genes_from_gbk_to_fasta(
    input_gbk: str,
    genes: str,
    n_before: int = 1,
    n_after: int = 1,
    output_fasta: str = "output.fasta",
):
    """
    Extracts genes and their neighboring genes from a GBK file, retrieves their amino acid sequences,
    and saves them to a FASTA file.
    -------------------------------
    Parameters:

    input_gbk: str
        - the path to the input GBK file.

    genes: str
        - a list of genes of interest. The program finds all matches
        (e.g., abc, abc_1, abc_A) and asks the user to confirm which genes to select.

    n_before: int (default = 1)
        - the number of neighboring genes to retrieve before the gene of interest.
    n_after: int (default = 1)
        - the number of neighboring genes to retrieve after the gene of interest.

    output_fasta: str
        - the path to the output FASTA file where the amino acid sequences will be saved.
    """

    def write_fasta_entry(outfile, gene_name, translation):
        outfile.write(f">{gene_name}\n")
        outfile.write(f"{translation}\n")

    with open(input_gbk, "r") as infile, open(output_fasta, "w") as outfile:
        content = infile.read()

        # rexex for genes and translation
        gene_pattern = re.compile(
            r"\/gene=\"(?P<gene>[^\"]+)\".*?"
            r"\/translation=\"(?P<translation>[A-Z\s]+)\"",
            re.DOTALL,
        )

        # search
        matches = list(gene_pattern.finditer(content))
        genes_found = []

        # find some variants
        for match in matches:
            gene_name = match.group("gene")
            if any(re.search(g, gene_name) for g in genes):
                genes_found.append((gene_name, match))

        if not genes_found:
            print(
                "No matching genes were found. \
            Try another gene or check the spelling and capitalization."
            )
            return

        # re-chek with user
        print("Matching genes:")
        for i, (gene_name, _) in enumerate(genes_found, 1):
            print(f"{i}: {gene_name}")

        selected_indices = input(
            "Enter the numbers of desired genes separated by commas (e.g., 1,3): "
        )
        selected_indices = [int(idx.strip()) - 1 for idx in selected_indices.split(",")]

        selected_genes = [genes_found[idx] for idx in selected_indices]

        # now searching for neighbours and writing out info
        for gene_name, match in selected_genes:
            gene_index = matches.index(match)
            start = max(0, gene_index - n_before)
            end = min(len(matches), gene_index + n_after + 1)

            for i in range(start, end):
                neighbour_gene = matches[i].group("gene")
                translation = matches[i].group("translation")
                translation = translation.replace("\n", "").replace(" ", "")
                write_fasta_entry(outfile, neighbour_gene, translation)

        print(f"Gene(s) and neighbours are now in {output_fasta}.")
