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
    StructField("accession_id", StringType(), True),
    StructField("feature_type", StringType(), True),
    StructField("start_pos", IntegerType(), True),
    StructField("end_pos", IntegerType(), True),
    StructField("organism_name", StringType(), True),
    StructField("protein_flag", BooleanType(), True),
])


def parse_single_record(record_text: str):
    """
    Parse a GenBank record chunk with Biopython.

    Returns: a list of dicts of relevant features.
    """
    parsed_features = []
    for biorec in SeqIO.parse(StringIO(record_text), "genbank"):
        org = biorec.annotations.get("organism", "")
        for feat in biorec.features:
            if feat.type not in FEATURES_TO_KEEP:
                continue

            start = int(feat.location.start)
            end = int(feat.location.end)
            if isinstance(feat.location, CompoundLocation):
                start = int(feat.location.parts[0].start)
                end = int(feat.location.parts[-1].end)

            is_protein = bool(feat.type == "CDS" and "protein_id" in feat.qualifiers)

            parsed_features.append({
                "accession_id": biorec.id,
                "feature_type": feat.type,
                "start_pos": start,
                "end_pos": end,
                "organism_name": org,
                "protein_flag": is_protein,
            })
    return parsed_features


def make_spark_session():
    """
    Creates and returns a Spark session with desired config.
    """
    spark_session = (
        SparkSession.builder
        .master("local[16]")
        .config("spark.executor.memory", "64g")
        .config("spark.driver.memory", "64g")
        .getOrCreate()
    )
    # Optionally reduce logging
    spark_session.sparkContext.setLogLevel("OFF")
    spark_session.conf.set("spark.task.maxBroadcastSize", "2m")
    return spark_session


def parse_gbff_to_df(spark_session, gbff_path):
    """
    1. Read entire GBFF file from `gbff_path`.
    2. Split into records by "//"
    3. Parallelize, parse each record with Biopython => flatten => Spark DF
    4. Filter out ambiguous (< or >) coords if necessary
    """
    with open(gbff_path, mode="r", encoding="utf-8") as fhandle:
        content = fhandle.read()

    # Each record ends with "//"
    raw_records = content.split("//\n")
    # Re-attach "//" so Biopython recognizes it
    records_cleaned = [r.strip() + "\n//" for r in raw_records if r.strip()]

    # Parallel parse
    rdd_parsed = spark_session.sparkContext.parallelize(records_cleaned, 16).flatMap(parse_single_record)
    df_full = spark_session.createDataFrame(rdd_parsed, schema=schema_for_features)

    df_ok = df_full.filter(F.col("start_pos") >= 0).filter(F.col("end_pos") >= 0)

    df_ok = df_ok.filter(F.col("feature_type").isin(FEATURES_TO_KEEP))

    return df_ok


def separate_genes(features_df):
    """
    Identify 'gene' features that do NOT match any 'CDS'
    on (accession, start, end).

    Returns a tuple: (df_with_cryptic_removed, df_cryptic_genes)

    """
    # Genes
    df_genes = features_df.filter(F.col("feature_type") == "gene").alias("g")
    # CDS
    df_cds = features_df.filter(F.col("feature_type") == "CDS").alias("c")

    # Genes that match a CDS by accession + start/end
    matched_genes = df_genes.join(
        df_cds,
        on=[
            df_genes.accession_id == df_cds.accession_id,
            df_genes.start_pos == df_cds.start_pos,
            df_genes.end_pos == df_cds.end_pos
        ],
        how="inner"
    ).select("g.*")

    # Cryptic = all genes minus matched
    cryptic_genes = df_genes.join(
        matched_genes,
        on=[
            df_genes.accession_id == matched_genes.accession_id,
            df_genes.start_pos == matched_genes.start_pos,
            df_genes.end_pos == matched_genes.end_pos,
        ],
        how="left_anti"
    )

    # Remove these cryptic genes from the main DF
    df_without_coding_genes = features_df.join(
        matched_genes,
        on=[
            features_df.feature_type == matched_genes.feature_type,
            features_df.accession_id == matched_genes.accession_id,
            features_df.start_pos == matched_genes.start_pos,
            features_df.end_pos == matched_genes.end_pos,
        ],
        how="left_anti"
    )

    return df_without_coding_genes, cryptic_genes

# Methods for awnsering the questions of assignment 5

def question_1(archaea_features):
    """
        Q1: How many features does an Archaea genome have on average?
        """

    avg_feats = (
        archaea_features.groupBy("accession_id")
        .count()
        .agg(F.avg("count"))
        .collect()[0][0]
    )
    print("Question 1: How many features does an Archaea genome have on average?")
    print(f"An Archaea genome has {round(avg_feats)} features on average")

def question_2(df_no_cds_genes, cryptic_df):
    """
        Q2: Ratio coding vs. non-coding features (coding / non-coding).

        """

    main_no_cryptic = df_no_cds_genes.join(cryptic_df, on="accession_id", how="left_anti")

    # non-coding => rRNA, ncRNA plus cryptic
    rna_cnt = main_no_cryptic.filter(F.col("feature_type").isin(["ncRNA", "rRNA"])).count()
    cryptic_cnt = cryptic_df.count()
    total_non_coding = rna_cnt + cryptic_cnt

    # coding => propeptide, gene, CDS in main_no_cryptic
    coding_cnt = main_no_cryptic.filter(
        F.col("feature_type").isin(["gene", "CDS", "propeptide"])
    ).count()

    ratio = coding_cnt / total_non_coding if total_non_coding else 0

    print("Question2: What is the ratio between coding and non-coding features? (coding / non-coding totals)")
    print(f"The ratio of coding/ non-coding = {round(ratio, 3)}.")

def question_3():

def question_4():

def question_5():


def main():
    """
    Main
    """


if __name__ == "__main__":
    main()
