# METADATA
__author__ = "Mats Slik"
__version__ = "0.6"

import argparse
import multiprocessing
import csv
import time
import os
import numpy as np


def compute_phred_scores(quality_str):
    """Convert a quality string into PHRED scores."""
    return np.array([ord(char) - 33 for char in quality_str])


def process_file(fastq_path, n_cores):
    """Process an individual FastQ file to calculate average PHRED scores."""
    with open(fastq_path, "r") as fastq_file:
        quality_lines = [line.strip() for i, line in enumerate(fastq_file) if i % 4 == 3]
    if not quality_lines:
        return fastq_path, None
    with multiprocessing.Pool(n_cores) as pool:
        results = pool.map(process_chunk, [quality_lines[i::n_cores] for i in range(n_cores)])
    return fastq_path, aggregate_results(results)


def process_chunk(quality_lines):
    """Process a chunk of quality lines to calculate cumulative PHRED scores."""
    phred_arrays = [compute_phred_scores(line.strip()) for line in quality_lines if line]
    if not phred_arrays:
        return None
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

    script_dir = os.path.dirname(os.path.abspath(__file__))

    results = []
    with multiprocessing.Pool(args.n) as pool:
        # Map each file to the process_file function
        result_objects = [pool.apply_async(process_file, args=(fastq_path, args.n)) for fastq_path in args.fastq_files]
        # Collect results as they come in
        for result in result_objects:
            results.append(result.get())

    for fastq_path, average_scores in results:
        if average_scores is not None:
            if args.o:
                output_file_name = os.path.join(script_dir, f"{os.path.basename(fastq_path)}.output.csv")
                with open(output_file_name, "w", newline="") as csvfile:
                    writer = csv.writer(csvfile)
                    writer.writerows(enumerate(average_scores))
                print(f"Output for {fastq_path} written to {output_file_name}")
            else:
                print(f"Results for {fastq_path}:")
                for index, score in enumerate(average_scores):
                    print(f"{index},{score}")
        else:
            print(f"No valid data to process in {fastq_path}")

    print(f"Execution time: {time.time() - start_time:.2f} seconds")


if __name__ == "__main__":
    main()
