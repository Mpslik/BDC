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

import os
import shutil

from io import StringIO
from Bio import SeqIO
from Bio.SeqFeature import CompoundLocation

from pyspark.sql import SparkSession
from pyspark.sql.types import (
    StructType,
    StructField,
    StringType,
    BooleanType,
    IntegerType,
)
from pyspark.sql.functions import (
    col,
    avg,
    count,
    expr,
    min as min_ps,
    max as max_ps,
)


FILENAME = "archaea.2.genomic.gbff"
FILE = "/data/datasets/NCBI/refseq/ftp.ncbi.nlm.nih.gov/refseq/release/archaea/" + FILENAME

FEATURES_OF_INTEREST = ["ncRNA", "rRNA", "gene", "propeptide", "CDS"]

# Spark session
spark = (
    SparkSession.builder.master("local[16]")
    .config("spark.executor.memory", "64g")
    .config("spark.driver.memory", "64g")
    .getOrCreate()
)
spark.conf.set("spark.task.maxBroadcastSize", "2m")
sc = spark.sparkContext
sc.setLogLevel("OFF")

# Schema for the final DataFrame
feature_schema = StructType(
    [
        StructField("accession", StringType(), True),
        StructField("type", StringType(), True),
        StructField("location_start", IntegerType(), True),
        StructField("location_end", IntegerType(), True),
        StructField("organism", StringType(), True),
        StructField("is_protein", BooleanType(), True),
    ]
)

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

