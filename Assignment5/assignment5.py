#!/usr/bin/env python3

"""
BDC Assignment 5 - MapReduce & PySpark

This script:
1) Uses PySpark to read and parse a .gbff file for Archaea.
2) Identifies cryptic genes (non-coding) vs. coding genes.
3) Answers 5 questions:
   Q1: Average number of features per Archaea genome?
   Q2: Ratio coding vs. non-coding features?
   Q3: Minimum and maximum number of proteins?
   Q4: Remove non-coding RNA features and store the remainder separately.
   Q5: Average feature length?
"""

# Metadata
__author__ = "Mats Slik"
__version__ = "0.1"

import os
import shutil

from io import StringIO
from pyspark.sql import SparkSession
from pyspark.sql import functions as F
from pyspark.sql.types import (
    StructType, StructField, StringType, BooleanType, IntegerType
)
from Bio import SeqIO
from Bio.SeqFeature import CompoundLocation


FILENAME = "archaea.2.genomic.gbff"
FILEPATH = "/data/datasets/NCBI/refseq/ftp.ncbi.nlm.nih.gov/refseq/release/archaea/" + FILENAME

FEATURES_TO_KEEP = ["ncRNA", "rRNA", "gene", "propeptide", "CDS"]

# Define Spark schema
schema_for_features = StructType([
    StructField("accession_id",   StringType(),  True),
    StructField("feature_type",   StringType(),  True),
    StructField("start_pos",      IntegerType(), True),
    StructField("end_pos",        IntegerType(), True),
    StructField("organism_name",  StringType(),  True),
    StructField("protein_flag",   BooleanType(), True),
])

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

