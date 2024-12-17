#!/usr/bin/env python3

"""
Calculate mean PHRED scores per base position from FastQ files using distributed computing.

Examples:
    $ python assignment2.py -s --host localhost --port 25715 --chunks 4 -o output.csv fastq_file1.fastq fastq_file2.fastq
    $ python assignment2.py -c --host localhost --port 25715 -n 4
"""

# Metadata
__author__ = "Mats Slik"
__version__ = "1.0"

import argparse
from collections import defaultdict
import multiprocessing as mp
from multiprocessing.managers import BaseManager
import numpy as np
import os
import queue
import sys
import time
import csv

# Constants
AUTHKEY = b"secretkey"
POISON_PILL = "TERMINATE"


# Utility functions
def compute_phred_scores(line):
    """Convert FastQ quality scores to numerical PHRED scores."""
    return np.array([ord(char) - 33 for char in line.strip()])


def parse_cli_args():
    """Parse command-line arguments to configure the script execution mode and parameters."""
    argparser = argparse.ArgumentParser(
        description="Distributed PHRED score calculator for FastQ files."
    )
    mode = argparser.add_mutually_exclusive_group(required=True)
    mode.add_argument(
        "-s",
        action="store_true",
        help="Run the program in Server mode; see extra options needed below",
    )
    mode.add_argument(
        "-c",
        action="store_true",
        help="Run the program in Client mode; see extra options needed below",
    )
    server_args = argparser.add_argument_group(title="Arguments when run in server mode")
    server_args.add_argument(
        "-o",
        action="store",
        dest="csvfile",
        type=argparse.FileType('w', encoding='UTF-8'),
        required=False,
        help="CSV file to save the output to. Default is output to terminal STDOUT",
    )
    server_args.add_argument(
        "fastq_files",
        action="store",
        type=argparse.FileType('r'),
        nargs="*",
        help="At least 1 Illumina Fastq Format file to process",
    )
    server_args.add_argument("--chunks", action="store", type=int, required=False, default=4)

    client_args = argparser.add_argument_group(title="Arguments when run in client mode")
    client_args.add_argument(
        "-n",
        action="store",
        dest="n",
        required=False,
        type=int,
        help="Number of cores to use per host.",
    )
    argparser.add_argument(
        "--host", action="store", type=str, help="The hostname where the Server is listening"
    )
    argparser.add_argument(
        "--port", action="store", type=int, help="The port on which the Server is listening"
    )

    return argparser.parse_args()


# Processing classes
class JobManager(BaseManager):
    """A custom manager to manage server/client communication."""


JobManager.register("get_job_queue")
JobManager.register("get_result_queue")


class PhredScoreCalculator:
    """Handles the computation of PHRED scores from chunks of FastQ data."""

    def __init__(self, file_paths, output_path=None, chunks=4):
        self.file_paths = file_paths
        self.output_path = output_path
        self.chunks = chunks

    def split_into_chunks(self, file_path):
        """Yield chunks of lines from a FastQ file."""
        with open(file_path, "r", encoding="utf-8") as file:
            lines = file.readlines()
            chunk_size = len(lines) // self.chunks + (len(lines) % self.chunks > 0)
            for i in range(0, len(lines), chunk_size):
                yield lines[i:i + chunk_size]

    @staticmethod
    def calculate_scores(chunk):
        """
        Calculate PHRED score sums and counts for each position from a list of FastQ quality lines.
        We store sums and counts so final averages can be computed after aggregating all chunks.
        """
        score_data = defaultdict(lambda: [0, 0])  # pos -> [sum, count]

        for line in chunk:
            if line.startswith("+") or line.startswith("@"):
                continue
            # Check if line is not a sequence line (A,T,C,G,N)
            scores = compute_phred_scores(line)
            for idx, score in enumerate(scores):
                score_data[idx][0] += score
                score_data[idx][1] += 1
        return dict(score_data)


# Server and Client implementations
class Server(mp.Process):
    """Server class for distributing jobs to clients."""

    def __init__(self, host, port, files, output, chunks):
        super().__init__()
        self.host = host
        self.port = port
        self.files = files
        self.output = output
        self.chunks = chunks

    def run(self):
        """Manage job distribution and result collection."""
        manager = JobManager(address=(self.host, self.port), authkey=AUTHKEY)
        manager.start()
        job_queue = manager.get_job_queue()
        result_queue = manager.get_result_queue()

        # Enqueue jobs
        for file_path in self.files:
            calculator = PhredScoreCalculator([file_path], chunks=self.chunks)
            for chunk in calculator.split_into_chunks(file_path):
                job_queue.put((PhredScoreCalculator.calculate_scores, chunk))

        # Collect results
        results = []
        expected_results = self.chunks * len(self.files)
        while len(results) < expected_results:
            result = result_queue.get()
            results.append(result)
            print("Collected a result.")

        # Output results
        if self.output:
            with open(self.output, "w", encoding="utf-8") as file:
                for result in results:
                    file.write(f"{result}\n")
        else:
            for result in results:
                print(result)

        # Signal termination to clients
        job_queue.put(POISON_PILL)
        manager.shutdown()


class Client(mp.Process):
    """Client class for processing jobs from the server."""

    def __init__(self, host, port, num_cores):
        super().__init__()
        self.host = host
        self.port = port
        self.num_cores = num_cores

    def run(self):
        """Process jobs from the server using available cores."""
        manager = JobManager(address=(self.host, self.port), authkey=AUTHKEY)
        manager.connect()
        job_queue = manager.get_job_queue()
        result_queue = manager.get_result_queue()

        while True:
            job = job_queue.get()
            if job == POISON_PILL:
                job_queue.put(POISON_PILL)
                break
            func, data = job
            result = func(data)
            result_queue.put(result)


# Main execution logic
def main():
    """Main function."""
    args = parse_cli_args()
    if args.server:
        server = Server(args.host, args.port, args.fastq_files, args.output, args.chunks)
        server.start()
        server.join()
    elif args.client:
        client = Client(args.host, args.port, args.num_cores)
        client.start()
        client.join()


if __name__ == "__main__":
    main()
