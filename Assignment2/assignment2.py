#!/usr/bin/env python3

"""
Calculate mean PHRED scores per base position from FastQ files using distributed computing.

Examples:
    $ python assignment2.py -s --host localhost --port 25715 --chunks 4 -o output.csv fastq_file1.fastq fastq_file2.fastq
    $ python assignment2.py -c --host localhost --port 25715 -n 4
"""

# Metadata
__author__ = "Mats Slik"
__version__ = "2.0"

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

g_job_queue = mp.Queue()
g_result_queue = mp.Queue()


# Utility functions
def compute_phred_scores(line):
    """Convert FastQ quality scores to numerical PHRED scores."""
    return np.array([ord(char) - 33 for char in line.strip()])


def read_all_quality_lines(fastq_path: str) -> list[str]:
    """
    Read all quality lines from the FastQ file, ensuring 4-line blocks:
      1) @read_header
      2) DNA sequence
      3) plus line
      4) quality line
    Return a list of *just* the quality lines (no partial reads).
    """
    quality_lines = []
    with open(fastq_path, "r", encoding="utf-8") as f:
        while True:
            header = f.readline()
            if not header:
                break
            seq_line = f.readline()
            plus_line = f.readline()
            qual_line = f.readline()
            if not qual_line:
                break
            # Only the quality line
            quality_lines.append(qual_line.rstrip("\n"))
    return quality_lines


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
    pass


JobManager.register("get_job_queue", callable=lambda: g_job_queue)
JobManager.register("get_result_queue", callable=lambda: g_result_queue)


class PhredScoreCalculator:
    """Handles the computation of PHRED scores from chunks of FastQ data."""

    def __init__(self, file_paths, output_path=None, chunks=4):
        self.file_paths = file_paths
        self.output_path = output_path
        self.chunks = chunks

    def split_into_chunks(self, file_handle):
        """
        read-boundary chunking.
        """
        file_path = file_handle.name
        quality_lines = read_all_quality_lines(file_path)

        n_lines = len(quality_lines)
        chunk_size = n_lines // self.chunks
        remainder = n_lines % self.chunks

        start_idx = 0
        for i in range(self.chunks):

            size = chunk_size + (1 if i < remainder else 0)
            if size == 0:
                break
            end_idx = start_idx + size
            yield quality_lines[start_idx:end_idx]
            start_idx = end_idx

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
        total_jobs = 0
        for file_handle in self.files:
            calculator = PhredScoreCalculator([file_handle], chunks=self.chunks)
            for chunk in calculator.split_into_chunks(file_handle):
                job_queue.put((PhredScoreCalculator.calculate_scores, chunk))
                total_jobs += 1

        # Collect results
        results = []
        while len(results) < total_jobs:
            result = result_queue.get()
            results.append(result)
            print("Collected a result.", file=sys.stderr)

        # Aggregate all results
        final_scores = defaultdict(lambda: [0, 0])
        for res in results:
            for pos, (s, c) in res.items():
                final_scores[pos][0] += s
                final_scores[pos][1] += c

        final_means = {
            pos: final_scores[pos][0] / final_scores[pos][1]
            for pos in final_scores
        }

        # Output results in CSV format
        if self.output:
            writer = csv.writer(self.output)
        else:
            writer = csv.writer(sys.stdout)
        writer.writerow(["line_number", "average_phredscore"])
        for pos in sorted(final_means.keys()):
            writer.writerow([pos, final_means[pos]])

        # Send termination signals to clients
        num_clients = mp.cpu_count()
        for _ in range(num_clients):
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
    if args.s:
        # Server mode
        if args.csvfile:
            output_file = args.csvfile
        else:
            output_file = None
        server = Server(args.host, args.port, args.fastq_files, output_file, args.chunks)
        server.start()
        server.join()
    elif args.c:
        # Client mode
        num_cores = args.n if args.n else mp.cpu_count()
        client = Client(args.host, args.port, num_cores)
        client.start()
        client.join()


if __name__ == "__main__":
    main()
