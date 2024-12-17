#!/bin/bash

# This script uses GNU parallel to process a FASTQ file in chunks.
# It relies on the Python script in chunk-parser mode to parse chunks,
# and then in combine-chunks mode to aggregate results.

# Input FASTQ file
FASTQ_FILE="/commons/Themas/Thema12/HPC/rnaseq.fastq"
NUM_JOBS=4
BLOCK_SIZE=1M

# Explanation:
# --jobs $NUM_JOBS: number of parallel jobs
# --pipepart: process a part of a file in parallel
# --block $BLOCK_SIZE: size of each chunk
# --sshlogin: multiple hosts can be specified like host1,host2
# --recstart '@': each record starts at a line starting with '@'
# We first run the python script in --chunk-parser mode on each chunk,
# then pipe all results to the same python script in --combine-chunks mode
# to produce a final CSV.

parallel --jobs $NUM_JOBS --pipepart --block $BLOCK_SIZE \
    --recstart '@' \
    python3 assignment3.py --chunk-parser :::: $FASTQ_FILE \
    | python3 assignment3.py --combine-chunks > output.csv