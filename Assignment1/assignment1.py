"""
This module provides functionality to calculate average PHRED scores from FastQ files.
The script processes multiple FastQ files simultaneously, efficiently handling large datasets using Python's multiprocessing capabilities.
The average PHRED scores are output either to the console or to CSV files for further analysis.
The computational workload scales dynamically based on the number of available CPU cores.

Features:
- Processes multiple FastQ files simultaneously.
- Efficiently handles large datasets using Python's multiprocessing capabilities.
- Outputs average PHRED scores either to the console or to CSV files for further analysis.
- Dynamically scales the computational workload based on the number of available CPU cores.

Dependencies:
- Python 3.11 or higher
- NumPy library

Usage:
The script is executed from the command line with the following syntax:

    python assignment1.py -n <number_of_cores> [-o] <fastq_file1> [<fastq_file2> ...]

Parameters:
- `-n` : Required. Specifies the number of CPU cores to use for processing.
- `-o` : Optional. If set, the program outputs the results to CSV files named after the input files with an added '.output.csv' suffix.
- `<fastq_files>` : Required. One or more paths to FastQ files that need to be processed.

Output:
- If `-o` is not specified, the average PHRED scores for each base position are printed to the console in the format "position,score".
- If `-o` is specified, each FastQ file's results are written to a corresponding CSV file with each line containing the base position and its average PHRED score.

Examples:
    python assignment1.py -n 4 sample1.fastq sample2.fastq
    python assignment1.py -n 8 -o large_dataset.fastq

Author: Mats Slik
Version: 1.
"""  # METADATA
# Imports
import argparse
import multiprocessing
import csv
import itertools
import numpy as np
from pathlib import Path
from typing import Generator, List, Tuple

__author__ = "Mats Slik"
__version__ = "1.2"


def compute_phred_scores(quality_str):
    """Convert a quality string into PHRED scores."""
    return np.array([ord(char) - 33 for char in quality_str], dtype=np.float64)


def quality_line_generator(start: int, stop: int, filepath: Path) -> Generator[str, None, None]:
    """A generator for FastQ file's quality lines."""
    with open(filepath, "rb") as fastq_file:
        fastq_file.seek(start)
        while fastq_file.tell() < stop:
            line = fastq_file.readline()
            if not line:
                break
            if line.startswith(b"@"):  # Assuming '@' is at the beginning of each read identifier
                fastq_file.readline()  # Skip the sequence line
                fastq_file.readline()  # Skip the '+' line
                yield fastq_file.readline().strip().decode()  # Yield the quality line


def process_chunk(chunk_tuple: Tuple[Path, int, int]) -> Tuple[Path, np.ndarray, np.ndarray]:
    """Process a chunk of a FastQ file to calculate cumulative PHRED scores."""
    filepath, start, stop = chunk_tuple

    count_iterator, parsing_iterator = itertools.tee(
        quality_line_generator(start, stop, filepath)
    )

    amount_of_lines = 0
    max_line_length = 0

    for line in count_iterator:
        amount_of_lines += 1
        if len(line) > max_line_length:
            max_line_length = len(line)

    all_phred_scores = np.full((amount_of_lines, max_line_length), np.nan, dtype=np.float64)

    for i, line in enumerate(parsing_iterator):
        for j, char in enumerate(line):
            all_phred_scores[i, j] = ord(char) - 33

    chunk_phred_sum = np.nansum(all_phred_scores, axis=0)
    chunk_phred_count = np.count_nonzero(~np.isnan(all_phred_scores), axis=0)

    return filepath, chunk_phred_sum, chunk_phred_count


def aggregate_results(results: List[Tuple[Path, np.ndarray, np.ndarray]], fastq_files: List[Path], output: bool):
    """Aggregate results from all chunks, calculating average PHRED scores."""
    for fastq_file in fastq_files:
        all_sum_arrays = [result[1] for result in results if result[0] == fastq_file]
        all_count_arrays = [result[2] for result in results if result[0] == fastq_file]

        total_phred_sums = np.sum(all_sum_arrays, axis=0)
        total_phred_counts = np.sum(all_count_arrays, axis=0)
        average_phred_scores = total_phred_sums / total_phred_counts

        if output:
            output_file_name = f"{fastq_file.stem}.output.csv"
            with open(output_file_name, "w", newline="") as csvfile:
                writer = csv.writer(csvfile)
                writer.writerows(enumerate(average_phred_scores))
            print(f"Output for {fastq_file} written to {output_file_name}")
        else:
            for index, score in enumerate(average_phred_scores):
                print(f"{index},{score}")


def get_chunks(file_paths: List[Path], number_of_chunks: int) -> List[Tuple[Path, int, int]]:
    """Divide files into chunks for parallel processing."""
    chunks = []
    chunks_per_file = number_of_chunks // len(file_paths)
    for file_path in file_paths:
        file_size = file_path.stat().st_size
        chunk_size = file_size // chunks_per_file
        start = 0
        stop = chunk_size
        while stop < file_size:
            chunks.append((file_path, start, stop))
            start = stop
            stop += chunk_size
        chunks.append((file_path, start, file_size))
    return chunks


def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description="Calculate average PHRED scores from FastQ files."
    )
    parser.add_argument("-n", required=True, type=int, help="Number of cores to use.")
    parser.add_argument(
        "-o", nargs='?', const=True, help="Output results to a CSV file or specify the output file name"
    )
    parser.add_argument("fastq_files", nargs="+", type=Path, help="FastQ files to process")
    args = parser.parse_args()

    # Validate the existence of each FastQ file
    for fastq_path in args.fastq_files:
        if not fastq_path.exists():
            print(f"Error: File {fastq_path} not found.")
            return

    # Create chunks for multiprocessing
    data_chunks = get_chunks(args.fastq_files, args.n)

    # Use multiprocessing pool to process chunks
    with multiprocessing.Pool(args.n) as pool:
        results = pool.map(process_chunk, data_chunks)

    # Aggregate and output results
    aggregate_results(results, args.fastq_files, output=args.o)


if __name__ == "__main__":
    main()
