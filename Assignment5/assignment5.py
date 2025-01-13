#!/usr/bin/env python3

"""
BDC Assignment 5 - MapReduce & PySpark

This script:
1. Creates a Spark Session (either local or on a cluster).
2. Reads a "gbff" file as text and infers a logical schema to extract genomic features.
3. Performs several analyses:
   - Average number of features per Archaea genome.
   - Ratio of coding vs. non-coding features.
   - Minimum and maximum number of proteins (coding sequences) across all organisms.
   - Removes all non-coding (RNA) features and writes them to a separate DataFrame/disk.
   - Computes the average length of a feature.
"""

# Metadata
__author__ = "Mats Slik"
__version__ = "0.1"


def parse_args():
    pass


def extract_records():
    pass


def parse_features():
    pass


def det_coding():
    pass


def calculate_statistics():
    pass


def main():
    """
    Main
    """


if __name__ == "__main__":
    main()

