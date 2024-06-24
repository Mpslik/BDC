#!/usr/bin/env python3

"""Calculate mean PHRED scores per base position from FastQ files using distributed computing.

Examples:
    $ python phred_score_calculator.py -s --host localhost --port 25715 [--chunks 4] [-o output.csv]
    fastq_file1.fastq [fastq_file2.fastq ...]
    $ python phred_score_calculator.py -c --host localhost --port 25715 [-n 4]
"""

# Metadata
__author__ = "Mats Slik"
__version__ = "0.1"

import argparse
from collections import defaultdict
import multiprocessing as mp
from multiprocessing.managers import BaseManager
import numpy as np
from pathlib import Path

# Constants
POISON_PILL = "TERMINATE"
AUTHKEY = b"secretkey"


# Utility functions
def compute_phred_scores(line):
    """Convert FastQ quality scores to numerical PHRED scores."""
    return np.array([ord(char) - 33 for char in line.strip()])


def parse_cli_args():
    """Parse command-line arguments to configure the script execution mode and parameters."""
    parser = argparse.ArgumentParser(
        description="Distributed PHRED score calculator for FastQ files."
    )
    parser.add_argument("--host", required=True, help="Hostname for the server.")
    parser.add_argument("--port", required=True, type=int, help="Port for server communication.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-s", "--server", action="store_true", help="Run in server mode.")
    group.add_argument("-c", "--client", action="store_true", help="Run in client mode.")
    parser.add_argument("--chunks", type=int, default=4, help="Number of chunks to split the input files.")
    parser.add_argument("-o", "--output", type=Path, help="Output file path (default is STDOUT).")
    parser.add_argument("-n", "--num_cores", type=int, default=mp.cpu_count(),
                        help="Number of cores to use in client mode.")
    parser.add_argument("fastq_files", nargs="*", type=Path, help="Paths to FastQ files to process.")
    return parser.parse_args()


# Processing classes
class JobManager(BaseManager):
    """A custom manager to manage server/client communication."""
    pass


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
        """Calculate average PHRED scores from a list of FastQ quality lines."""
        score_sums = defaultdict(list)
        for line in chunk:
            if line.startswith("+"):
                continue  # Skip header lines
            scores = compute_phred_scores(line)
            for idx, score in enumerate(scores):
                score_sums[idx].append(score)
        return {pos: np.mean(scores) for pos, scores in score_sums.items()}


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
                job_queue.put((calculator.calculate_scores, chunk))

        # Collect results
        results = []
        while len(results) < len(self.files) * self.chunks:
            results.append(result_queue.get())
            print("Collected a result.")

        # Output results
        if self.output:
            with open(self.output, "w", encoding="utf-8") as file:
                for result in results:
                    file.write(f"{result}\n")
        else:
            for result in results:
                print(result)

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
    elif args.client:
        client = Client(args.host, args.port, args.num_cores)
        client.start()
        client.join()


if __name__ == "__main__":
    main()
