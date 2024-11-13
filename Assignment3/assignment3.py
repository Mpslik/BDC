#!/usr/bin/env python3

"""
Calculate average PHRED scores per base position from FastQ files using GNU Parallel
"""

import argparse
import sys
import numpy as np
import csv
from pathlib import Path
from typing import List, Tuple


def compute_phred_scores(quality_str):
    """
    Convert a quality string into PHRED scores.
    Each character in the quality string represents a PHRED score
    (ASCII-encoded). This function converts it to a numerical array.
    """
    return np.array([ord(char) - 33 for char in quality_str], dtype=np.float64)


def get_chunks(file_paths: List[Path], number_of_chunks: int) -> List[Tuple[Path, int, int]]:
    """
    Divide each file into chunks for parallel processing.
    Each chunk is defined by a start and stop byte position within the file.
    This enables each chunk to be processed independently.
    """
    chunks = []
    chunks_per_file = number_of_chunks // len(file_paths)  # Determine number of chunks per file

    for file_path in file_paths:
        file_size = file_path.stat().st_size  # Get file size in bytes
        chunk_size = file_size // chunks_per_file  # Calculate chunk size in bytes
        start = 0
        stop = chunk_size

        # Divide file into chunks based on byte size
        while stop < file_size:
            chunks.append((file_path, start, stop))
            start = stop
            stop += chunk_size
        chunks.append((file_path, start, file_size))  # Final chunk up to end of file
    return chunks


def process_chunk(chunk_tuple: Tuple[Path, int, int]) -> Tuple[np.ndarray, np.ndarray]:
    """
    Process a single chunk of a FastQ file to calculate cumulative PHRED scores.
    Reads each quality line in the chunk, converts it to PHRED scores, and
    accumulates the scores and counts for each base position.
    """
    filepath, start, stop = chunk_tuple  # Unpack chunk information
    all_phred_scores = []  # List to hold PHRED scores for each line in chunk

    # Open file and move to the start position for this chunk
    with open(filepath, "r") as fastq_file:
        fastq_file.seek(start)

        # Read lines until reaching the stop position for this chunk
        for line in fastq_file:
            if line.startswith("@"):  # Detect quality line sequence start
                fastq_file.readline()  # Skip the sequence line
                fastq_file.readline()  # Skip the '+' line
                quality_line = fastq_file.readline().strip()  # Read quality line
                all_phred_scores.append(compute_phred_scores(quality_line))  # Convert to PHRED scores
            if fastq_file.tell() >= stop:  # Stop if we have reached the chunk's end
                break

    # Convert list of scores to a numpy array for further processing
    scores = np.array(all_phred_scores)
    return np.nansum(scores, axis=0), np.count_nonzero(~np.isnan(scores), axis=0)  # Sum and count non-NaN values


def aggregate_results(results: List[Tuple[np.ndarray, np.ndarray]]):
    """
    Aggregate results from all chunks to compute average PHRED scores per base position.
    Takes cumulative scores and counts from each chunk, combines them,
    and calculates the average score per position, outputting as CSV.
    """
    summed_scores = np.sum([res[0] for res in results], axis=0)  # Sum scores across all chunks
    summed_counts = np.sum([res[1] for res in results], axis=0)  # Sum counts across all chunks
    avg_scores = summed_scores / summed_counts  # Compute average score per position

    # Output results as CSV format to stdout
    writer = csv.writer(sys.stdout)
    writer.writerow(["Position", "Average PHRED Score"])
    writer.writerows(enumerate(avg_scores))


def main():
    """
    Main function to handle argument parsing and execute processing.
    Divides a FastQ file into chunks, processes each chunk, and aggregates results.
    """
    parser = argparse.ArgumentParser(description="GNU Parallel PHRED score calculator for FastQ files.")
    parser.add_argument("fastq_file", type=Path, help="Path to the FastQ file to process.")
    parser.add_argument("-n", type=int, required=True, help="Number of chunks to split the file into.")
    args = parser.parse_args()

    # Create chunks for parallel processing
    data_chunks = get_chunks([args.fastq_file], args.n)

    # Process each chunk individually
    results = [process_chunk(chunk) for chunk in data_chunks]

    # Aggregate and output results to stdout
    aggregate_results(results)


if __name__ == "__main__":
    main()