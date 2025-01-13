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

from pyspark.shell import spark
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
schema_for_features = StructType([
    StructField("accession_id", StringType(), True),
    StructField("feature_type", StringType(), True),
    StructField("start_pos", IntegerType(), True),
    StructField("end_pos", IntegerType(), True),
    StructField("organism_name", StringType(), True),
    StructField("protein_flag", BooleanType(), True),
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
    spark.sparkContext.setLogLevel("WARN")
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
            yield from parse_record_text(record_text)


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
                start_val = int(feat.location.parts[0].start)
                end_val = int(feat.location.parts[-1].end)
            else:
                start_val = int(feat.location.start)
                end_val = int(feat.location.end)

            protein_bool = (feat.type == "CDS") and ("protein_id" in feat.qualifiers)

            yield {
                "accession_id": biorec.id,
                "feature_type": feat.type,
                "start_pos": start_val,
                "end_pos": end_val,
                "organism_name": organism_val,
                "protein_flag": protein_bool,
            }


def parse_gbff_to_df(spark_session, gbff_path):
    """
    Reads the GBFF file line-by-line, accumulates lines per record,
    parses with Biopython, and returns a Spark DF.
    """
    # Force more partitions if the file is large
    lines_rdd = spark_session.sparkContext.textFile(gbff_path, minPartitions=100)

    # For each partition, parse the lines
    parsed_rdd = lines_rdd.mapPartitions(parse_partition_of_lines)

    # Convert to Spark DataFrame
    df_full = spark_session.createDataFrame(parsed_rdd, schema=schema_for_features)

    # Filter out invalid coords
    df_ok = df_full.filter((F.col("start_pos") >= 0) & (F.col("end_pos") >= 0))

    # Keep relevant features
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

    # Remove cryptic genes from the main DF
    df_without_coding_genes = features_df.join(
        matched_genes,
        on=[
            features_df.feature_type == matched_genes.feature_type,
            features_df.accession_id == matched_genes.accession_id,
            features_df.start_pos == matched_genes.start_pos,
            features_df.end_pos == matched_genes.end_pos
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
        archaea_features.groupBy("accession_id").count()
        .agg(F.avg("count"))
        .first()[0]
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


def question_3(features_df):
    """
    Q3: Minimum and maximum number of proteins (CDS with protein_flag=True).
    """
    only_prot_cds = features_df.filter(
        (F.col("protein_flag") == True) & (F.col("feature_type") == "CDS")
    )
    org_counts = only_prot_cds.groupBy("organism_name").agg(F.count("*").alias("count"))
    min_max = org_counts.agg(
        F.min("count").alias("minimum"),
        F.max("count").alias("maximum")
    ).first()
    print("Question3: What are the minimum and maximum number of proteins of all organisms in the file?")
    print(f"Minimum = {min_max['minimum']} and maximum = {min_max['maximum']}")


def question_4(all_features, cryptic_df):
    """
    Q4: Remove all non-coding (RNA) features => (rRNA, ncRNA) plus cryptic genes.
    """
    keep_code = all_features.filter(~(F.col("feature_type").isin(["rRNA", "ncRNA"])))
    keep_code = keep_code.join(cryptic_df, on="accession_id", how="left_anti")

    table_name = "coding_spark_frame"
    spark.catalog.clearCache()

    # Remove old table (dir) if it exists
    if os.path.exists(f"spark-warehouse/{table_name}"):
        shutil.rmtree(f"spark-warehouse/{table_name}")

    # Save as table
    keep_code.write.saveAsTable(table_name)
    print("Question4: Remove all non-coding (RNA) features and write this as a separate DataFrame (Spark format)")
    print(f"Wrote table -> {table_name}")


def question_5(df_any):
    row_avg = (
        df_any.withColumn("length", (F.col("end_pos") - F.col("start_pos")))
              .agg(F.avg("length").alias("avg_length"))
              .first()["avg_length"]
    )
    print("Question5: What is the average length of a feature?")
    print(f"Average length is: {round(row_avg)}")


def main():
    spark_sess = make_spark_session()

    # 1) Read + parse
    df_parsed = parse_gbff_to_df(spark_sess, FILEPATH)

    # 2) Identify cryptic genes => remove from main DF
    df_no_coding_genes, cryptic_genes = separate_genes(df_parsed)

    # 3) Answer questions
    question_1(df_no_coding_genes)
    print("\n")
    question_2(df_no_coding_genes, cryptic_genes)
    print("\n")
    question_3(df_no_coding_genes)
    print("\n")
    question_4(df_no_coding_genes, cryptic_genes)
    print("\n")
    question_5(df_no_coding_genes)

    spark_sess.stop()

    # Done
    spark_sess.stop()


if __name__ == "__main__":
    main()
