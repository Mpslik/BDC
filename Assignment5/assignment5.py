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
__version__ = "0.2"

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

FILENAME = "archaea.3.genomic.gbff"
FILEPATH = "/data/datasets/NCBI/refseq/ftp.ncbi.nlm.nih.gov/refseq/release/archaea/" + FILENAME
FEATURES_TO_KEEP = ["ncRNA", "rRNA", "gene", "propeptide", "CDS"]

# Define Spark schema
schema_for_features_str = StructType([
    StructField("accession", StringType(), True),
    StructField("type", StringType(), True),
    StructField("location_start", IntegerType(), True),
    StructField("location_end", IntegerType(), True),
    StructField("organism", StringType(), True),
    StructField("is_protein", BooleanType(), True),
])


def make_spark_session():
    """
    Creates and returns a Spark session with desired config.
    """
    spark = (
        SparkSession.builder
        .appName("assignment5_mats")
        .master("local[16]")
        .config("spark.driver.memory", "64g")
        .config("spark.executor.memory", "64g")
        .getOrCreate())
    spark.conf.set("spark.task.maxBroadcastSize", "16m")
    spark.sparkContext.setLogLevel("OFF")
    return spark


def parse_partition_of_lines(lines_iter):
    """
    Partition-level function that accumulates lines until we see '//'
    and parses that record chunk, yielding feature dicts.
    """
    buffer = []
    for line in lines_iter:
        buffer.append(line)
        if line.strip() == "//":
            record_text = "\n".join(buffer)
            buffer.clear()
            for feat_dict in parse_record_text(record_text):
                yield feat_dict


def parse_record_text(record_chunk):
    """
    Generator function that parses a single record chunk using Biopython
    and yields feature dictionaries.
    """
    for biorec in SeqIO.parse(StringIO(record_chunk), "genbank"):
        organism_val = biorec.annotations.get("organism", "")
        for feat in biorec.features:
            if feat.type not in FEATURES_TO_KEEP:
                continue

            if isinstance(feat.location, CompoundLocation):
                start = int(feat.location.parts[0].start)
                end = int(feat.location.parts[-1].end)
            else:
                start = int(feat.location.start)
                end = int(feat.location.end)

            is_protein = (feat.type == "CDS") and ("protein_id" in feat.qualifiers)
            yield {
                "accession": biorec.id,
                "type": feat.type,
                "location_start": start,
                "location_end": end,
                "organism": organism_val,
                "is_protein": is_protein,
            }


def parse_gbff_to_df(spark_session, gbff_path):
    """
    Reads the GBFF file line-by-line, accumulates lines per record,
    parses with Biopython, and returns a Spark DF.
    """
    lines_rdd = spark_session.sparkContext.textFile(gbff_path)
    parsed_rdd = lines_rdd.mapPartitions(parse_partition_of_lines)
    df_raw = spark_session.createDataFrame(parsed_rdd, schema=schema_for_features_str)

    df_filtered = (
        df_raw
        .filter(~F.col("location_start").startswith("<"))
        .filter(~F.col("location_end").startswith(">"))
        .filter(F.col("type").isin(FEATURES_TO_KEEP))
    )
    return df_filtered


def filter_coding_and_cryptic(df_features):
    """
    Filters the DataFrame to separate cryptic genes (non-coding) from coding genes.
    """
    # All genes and CDS
    all_genes = df_features.filter(F.col("type") == "gene").alias("genes")
    all_cds = df_features.filter(F.col("type") == "CDS").alias("cds")

    # Coding genes: genes that have matching CDS
    coding_genes = all_genes.join(
        all_cds,
        on=[
            all_genes.accession == all_cds.accession,
            all_genes.location_start == all_cds.location_start,
            all_genes.location_end == all_cds.location_end,
        ]
    ).select(all_genes["*"])

    # Cryptic genes: genes that are not coding
    cryptic_genes = all_genes.join(
        coding_genes,
        on=[
            all_genes.accession == coding_genes.accession,
            all_genes.location_start == coding_genes.location_start,
            all_genes.location_end == coding_genes.location_end,
        ],
        how="left_anti"
    )

    # Features excluding cryptic genes
    df_without_cryptic = df_features.join(
        cryptic_genes,
        on=[
            df_features.accession == cryptic_genes.accession,
            df_features.location_start == cryptic_genes.location_start,
            df_features.location_end == cryptic_genes.location_end,
        ],
        how="left_anti"
    )

    return df_without_cryptic, cryptic_genes


# Define methods for answering questions
def question_1(df):
    avg_features = (
        df.groupBy("accession").count()
        .agg(F.avg("count").alias("avg_features"))
        .first()["avg_features"]
    )
    print(f"Q1: Average features per genome: {round(avg_features)}")


def question_2(df, cryptic_df):
    non_coding_count = (
            df.filter(F.col("type").isin(["ncRNA", "rRNA"])).count() + cryptic_df.count()
    )
    coding_count = df.filter(F.col("type").isin(["gene", "CDS", "propeptide"])).count()
    ratio = coding_count / non_coding_count if non_coding_count else 0
    print(f"Q2: Coding/Non-Coding Ratio: {round(ratio, 3)}")


def question_3(df):
    proteins = df.filter((F.col("type") == "CDS") & (F.col("is_protein") == True))
    org_counts = proteins.groupBy("organism").count()
    min_max = org_counts.agg(F.min("count").alias("min"), F.max("count").alias("max")).first()
    print(f"Q3: Min proteins: {min_max['min']}, Max proteins: {min_max['max']}")


def question_4(df, cryptic_df):
    coding_only = (
        df.filter(~F.col("type").isin(["rRNA", "ncRNA"]))
        .join(cryptic_df, on="accession", how="left_anti")
    )
    table_name = "coding_spark_frame"
    table_path = f"spark-warehouse/{table_name}"

    # Remove directory if it exists
    if os.path.exists(table_path):
        shutil.rmtree(table_path)

    coding_only.write.saveAsTable("coding_spark_frame")
    print("Q4: Saved coding features as a separate table.")


def question_5(df):
    avg_length = (
        df.withColumn("length", F.col("location_end") - F.col("location_start"))
        .agg(F.avg("length").alias("avg_length"))
        .first()["avg_length"]
    )
    print(f"Q5: Average feature length: {round(avg_length)}")


# Main execution block
def main():
    spark_sess = make_spark_session()
    df_features = parse_gbff_to_df(spark_sess, FILEPATH)
    df_without_cryptic, cryptic_genes = filter_coding_and_cryptic(df_features)

    question_1(df_without_cryptic)
    question_2(df_without_cryptic, cryptic_genes)
    question_3(df_without_cryptic)
    question_4(df_without_cryptic, cryptic_genes)
    question_5(df_without_cryptic)

    spark_sess.stop()


if __name__ == "__main__":
    main()
