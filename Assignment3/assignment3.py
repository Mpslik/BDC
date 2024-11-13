#!/usr/bin/env python3

"""
Calculate average PHRED scores per base position from FastQ files using GNU Parallel
"""

import argparse
import numpy as np
import csv
from pathlib import Path
from typing import List, Tuple


def compute_phred_scores(quality_str):
    """Convert a quality string into PHRED scores."""
    return np.array([ord(char) - 33 for char in quality_str], dtype=np.float64)


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


def process_chunk(chunk_tuple: Tuple[Path, int, int]) -> Tuple[np.ndarray, np.ndarray]:
    filepath, start, stop = chunk_tuple
    all_phred_scores = []
    with open(filepath, "r") as fastq_file:
        fastq_file.seek(start)
        for line in fastq_file:
            if line.startswith("@"):
                fastq_file.readline()
                fastq_file.readline()
                quality_line = fastq_file.readline().strip()
                all_phred_scores.append(compute_phred_scores(quality_line))
            if fastq_file.tell() >= stop:
                break
    scores = np.array(all_phred_scores)
    return np.nansum(scores, axis=0), np.count_nonzero(~np.isnan(scores), axis=0)
