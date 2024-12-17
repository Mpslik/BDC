#!/bin/bash

# Define the FastQ file to process

FASTQ_FILE="/commons/Themas/Thema12/HPC/rnaseq.fastq"

NUM_JOBS=4  # Number of parallel jobs
BLOCK_SIZE=1M  # Block size

parallel --jobs $NUM_JOBS --pipepart --block $BLOCK_SIZE --sshlogin assemblix2019,nuc109 \
    --recstart '@' \
    python3 assignment3.py --chunk-parser | \
    python3 assignment3.py --combine-chunks > output.csv