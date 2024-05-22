# METADATA
__author__ = "Mats Slik"
__version__ = "0.6"


# Imports
import argparse
import multiprocessing
import csv
import time
import os
import numpy as np

# Code

def compute_phred_scores(quality_str):
    """Convert a quality string into PHRED scores."""
    return np.array([ord(char) - 33 for char in quality_str])

def process_chunk(quality_lines):
    """Process a chunk of quality lines to calculate cumulative PHRED scores."""
    phred_arrays = [compute_phred_scores(line.strip()) for line in quality_lines if line]
    if not phred_arrays:
        return None
    # Use numpy to sum up PHRED scores for better performance
    combined_array = np.stack(phred_arrays)
    sum_scores = np.sum(combined_array, axis=0)
    return sum_scores, len(phred_arrays)

def aggregate_results(results):
    """Aggregate results from all chunks, calculating average PHRED scores."""
    total_sum = None
    total_count = 0
    for result in results:
        if result is None:
            continue
        sum_scores, count = result
        if total_sum is None:
            total_sum = sum_scores
        else:
            total_sum += sum_scores
        total_count += count
    return total_sum / total_count if total_count > 0 else None

def main():
    """Main function."""
    start_time = time.time()

    parser = argparse.ArgumentParser(description="Calculate average PHRED scores from FastQ files.")
    parser.add_argument("-n", required=True, type=int, help="Number of cores to use.")
    parser.add_argument("-o", action="store_true", help="Output results to individual CSV files")
    parser.add_argument("fastq_files", nargs="+", help="FastQ files to process")
    args = parser.parse_args()

    with multiprocessing.Pool(args.n) as pool:
        for fastq_path in args.fastq_files:
            with open(fastq_path, "r") as fastq_file:
                quality_lines = [line.strip() for i, line in enumerate(fastq_file) if i % 4 == 3]
                results = pool.map(process_chunk, [quality_lines[i::args.n] for i in range(args.n)])
                average_scores = aggregate_results(results)

                if average_scores is not None:
                    if args.o:
                        output_file_name = f"{os.path.splitext(fastq_path)[0]}.output.csv"
                        with open(output_file_name, "w", newline="") as csvfile:
                            writer = csv.writer(csvfile)
                            writer.writerows(enumerate(average_scores))
                        print(f"Output for {fastq_path} written to {output_file_name}")
                    else:
                        for index, score in enumerate(average_scores):
                            print(f"{index},{score}")

    print(f"Execution time: {time.time() - start_time:.2f} seconds")

if __name__ == "__main__":
    main()