import sys
import argparse
from Bio import SeqIO
from Bio.SeqUtils import GC
from loguru import logger


def setup_logger(log_level="INFO"):
    logs_format = "<green>{time:YYYY-MM-DD HH:mm:ss.SSS}</green> | <level>{level: <8}</level> | <cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan> - <level>{message}</level>"

    logger.remove()
    logger.add(sys.stdout, level=log_level, format=logs_format, colorize=True)
    logger.add("all_logs.log", level=0, format=logs_format, colorize=False)
    logger.add(
        "errors_logs.log",
        level=40,
        format=logs_format,
        colorize=False,
        rotation="500 KB",
    )


def fastq_filter_tool(
    input_fastq: str,
    output_fastq: str,
    gc_bounds: float | tuple[float, float] = (0, 100),
    length_bounds: int | tuple[int, int] = (0, 2 ** 32),
    quality_threshold: int = 0,
):
    logger.info(
        f"💅🏻Started to filter {input_fastq} with GC bounds {gc_bounds}, length bounds {length_bounds}, and quality threshold {quality_threshold}"
    )

    if isinstance(gc_bounds, (int, float)):
        gc_bounds = (0, gc_bounds)
    if isinstance(length_bounds, int):
        length_bounds = (0, length_bounds)

    count = 0
    try:
        with open(output_fastq, "w") as output_handle:
            for record in SeqIO.parse(input_fastq, "fastq"):
                sequence = record.seq
                quality = record.letter_annotations["phred_quality"]

                if (
                    gc_bounds[0] <= GC(sequence) <= gc_bounds[1]
                    and length_bounds[0] <= len(sequence) <= length_bounds[1]
                    and sum(quality) / len(quality) >= quality_threshold
                ):
                    SeqIO.write(record, output_handle, "fastq")
                    count += 1

            if count == 0:
                logger.warning("😩 The reads were filtered out... Survival rate: 0%.")
            else:
                logger.info(f"🌞Filtering completed. {count} reads passed the filter")

    except Exception as e:
        logger.error(f"💥An error occurred during filtering: {e}")
        raise


def parse_args():
    parser = argparse.ArgumentParser(
        prog="fastq_filter_tool",
        description="Filter FASTQ reads by GC content, length, and quality.",
        epilog="💚 Developed during IB2024-2025 Python Course by aryamakeeva",
    )

    parser.add_argument("input_fastq", help="Path to the input FASTQ file.")
    parser.add_argument("output_fastq", help="Path to save the filtered FASTQ file.")
    parser.add_argument(
        "--gc_bounds",
        nargs=2,
        type=float,
        default=(0, 100),
        help="Min and max GC content, e.g. --gc_bounds 30 70 (default: 0 100)",
    )
    parser.add_argument(
        "--length_bounds",
        nargs=2,
        type=int,
        default=(0, 2 ** 32),
        help="Min and max sequence length (default: all sequenses)",
    )
    parser.add_argument(
        "--quality_threshold",
        type=int,
        default=0,
        help="Minimum average quality score (default: 0)",
    )
    parser.add_argument(
        "--log_level",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        default="INFO",
        help="Logging level (default: INFO)",
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    setup_logger(args.log_level)

    fastq_filter_tool(
        input_fastq=args.input_fastq,
        output_fastq=args.output_fastq,
        gc_bounds=tuple(args.gc_bounds),
        length_bounds=tuple(args.length_bounds),
        quality_threshold=args.quality_threshold,
    )
