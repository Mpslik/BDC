#!/bin/bash

# Define the FastQ file to process

FASTQ_FILE=""

# Define the number of jobs and block size for parallel processing
NUM_JOBS=4  # Number of parallel jobs (adjust as needed)
BLOCK_SIZE=1M  # Block size (adjust based on file size and performance)

parallel --jobs $NUM_JOBS --pipepart --block $BLOCK_SIZE --sshlogin assemblix2019,nuc109 \
    --recstart '@' \
    python3 assignment3.py --chunk-parser | \
    python3 assignment3.py --combine-chunks > output.csv