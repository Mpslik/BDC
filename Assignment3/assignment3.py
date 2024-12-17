#!/usr/bin/env python3

"""
Calculate average PHRED scores per base position from FastQ files using GNU Parallel.
This script has two modes:
1. --chunk-parser: Reads a chunk of FASTQ data from stdin, computes PHRED sums and counts.
2. --combine-chunks: Reads sums and counts from stdin and combines them into final averages.
"""

# Metadata
__author__ = "Mats Slik"
__version__ = "1.0"

import argparse
import sys
import numpy as np
import csv
from pathlib import Path
from typing import List, Tuple


def argument_parser():
    parser = argparse.ArgumentParser(
        description="GNU Parallel PHRED score calculator for FastQ files."
    )
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument(
        "--chunk-parser",
        action="store_true",
        help="Parse a chunk of FASTQ data from stdin and output sum and count per position."
    )
    mode.add_argument(
        "--combine-chunks",
        action="store_true",
        help="Combine results from chunk parses (sum and count arrays) into final CSV output."
    )
    return parser.parse_args()


def compute_phred_scores(quality_str):
    """
    Convert a quality string into PHRED scores.
    Each character in the quality string represents a PHRED score
    (ASCII-encoded). This function converts it to a numerical array.
    """
    return np.array([ord(char) - 33 for char in quality_str], dtype=np.float64)

def read_fastq_quality_lines():
    """
    Reads FASTQ data from stdin and yields only the quality lines.
    """
    lines = sys.stdin.read().strip().split('\n')
    # Process lines in sets of four
    for i in range(0, len(lines), 4):
        if i+3 < len(lines):
            quality_line = lines[i+3].strip()
            yield quality_line

def aggregate_results_from_stdin():
    """
    Aggregate results from stdin to compute average PHRED scores per base position.
    Input lines should contain 'sum: [...]' and 'count: [...]' arrays from multiple chunks.
    This function parses these arrays, aggregates sums and counts, and then outputs averages.
    """
    sum_arrays = []
    count_arrays = []

    for line in sys.stdin:
        line = line.strip()
        if line.startswith("sum:"):
            arr_str = line[line.index('[')+1:line.rindex(']')]
            if arr_str.strip():
                sums = np.fromstring(arr_str, sep=",")
                sum_arrays.append(sums)
            else:
                sum_arrays.append(np.array([]))
        elif line.startswith("count:"):
            arr_str = line[line.index('[')+1:line.rindex(']')]
            if arr_str.strip():
                counts = np.fromstring(arr_str, sep=",")
                count_arrays.append(counts)
            else:
                count_arrays.append(np.array([]))

    # If no data, just print a header
    writer = csv.writer(sys.stdout)
    writer.writerow(["Position", "Average PHRED Score"])

    if not sum_arrays or not count_arrays:
        return

    # Pad arrays to the max length
    max_len = 0
    for arr in sum_arrays + count_arrays:
        if len(arr) > max_len:
            max_len = len(arr)

    def pad_arrays(arr_list, length):
        padded = np.zeros((len(arr_list), length))
        for i, arr in enumerate(arr_list):
            padded[i, :len(arr)] = arr
        return padded

    total_sums = np.sum(pad_arrays(sum_arrays, max_len), axis=0)
    total_counts = np.sum(pad_arrays(count_arrays, max_len), axis=0)
    with np.errstate(divide='ignore', invalid='ignore'):
        averages = total_sums / total_counts
        averages[np.isnan(averages)] = 0.0

    for i, val in enumerate(averages):
        writer.writerow([i, val])


def chunk_parser_mode():
    """
    In chunk-parser mode, read a chunk of FASTQ data from stdin,
    compute cumulative PHRED sums and counts per base position,
    and print them to stdout in a parseable format.
    """
    all_scores = []
    for q_line in read_fastq_quality_lines():
        if q_line:
            scores = compute_phred_scores(q_line)
            all_scores.append(scores)

    if len(all_scores) == 0:
        print("sum: []")
        print("count: []")
        return

    max_len = max(len(arr) for arr in all_scores)
    padded = np.zeros((len(all_scores), max_len), dtype=np.float64)
    for i, arr in enumerate(all_scores):
        padded[i, :len(arr)] = arr

    pos_sum = np.sum(padded, axis=0)
    pos_count = np.count_nonzero(padded, axis=0)

    print("sum:", list(pos_sum))
    print("count:", list(pos_count))


def main():
    """
    Main function:
    Based on the chosen mode, either parse a chunk of FASTQ data for sums and counts,
    or combine these chunks to produce a final CSV of average PHRED scores.
    """
    args = argument_parser()
    if args.chunk_parser:
        chunk_parser_mode()
    elif args.combine_chunks:
        aggregate_results_from_stdin()

if __name__ == "__main__":
    main()
