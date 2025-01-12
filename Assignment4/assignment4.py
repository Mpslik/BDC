#!/usr/bin/env python3

"""
Calculate average PHRED scores per base position from FastQ files using MPI (mpi4py).

This script:
1. Uses MPI with one controller (rank 0) and multiple worker ranks.
2. Distributes chunks of FASTQ files among workers.
3. Collects results and writes them to a CSV or to stdout.
4. Tracks runtime and allows specification of run index and worker count
   for performance and numerical stability experiments.
"""

# Metadata
__author__ = "Mats Slik"
__version__ = "0.1"

import argparse
import csv
import sys
from collections import defaultdict
from mpi4py import MPI
from pathlib import Path
import numpy as np


def zero_pair():
    """Return [0.0, 0.0] to avoid unpickleable lambdas."""
    return [0.0, 0.0]


def parse_args():
    """
    Parse command-line arguments to configure the script's execution parameters.
    """
    parser = argparse.ArgumentParser(
        description="Calculate average PHRED scores using MPI, naive byte chunking."
    )
    parser.add_argument(
        "--files",
        nargs="+",
        required=True,
        help="List of FastQ files to process."
    )
    parser.add_argument(
        "-n",
        type=int,
        default=4,
        help="Number of chunks per file ."
    )
    parser.add_argument(
        "--run-index",
        type=int,
        default=1,
        help="Which run index ."
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=0,
        help="Reported worker count."
    )
    parser.add_argument(
        "--output-prefix",
        type=str,
        default=None,
        help="If set, output results to CSVs with this prefix. Otherwise, to STDOUT."
    )
    return parser.parse_args()


class PhredscoreCalculator:
    """
    A class responsible for reading FastQ file slices (chunked by bytes),
    extracting PHRED scores, and computing partial sums.

    Attributes:
        n_chunks (int): Number of chunks to split each file into (naive byte-based).
    """

    def __init__(self, n_chunks):
        self.n_chunks = n_chunks

    @staticmethod
    def compute_phred_scores(quality_str):
        """
        Convert a single quality string into PHRED scores.
        """
        return np.array([ord(ch) - 33 for ch in quality_str.strip()], dtype=np.float64)

    @staticmethod
    def read_chunk_bytes(filepath: Path, start: int, end: int) -> bytes:
        """
        Read raw bytes from 'start' to 'end' in the file.
        """
        with open(filepath, "rb") as f:
            f.seek(start)
            return f.read(end - start)

    def chunk_file(self, file_path: Path):
        """
        Split 'file_path' into self.n_chunks by naive byte-splitting.

        Yields (file_path, start, end).
        """
        fsize = file_path.stat().st_size
        chunk_size = fsize // self.n_chunks
        remainder = fsize % self.n_chunks

        offset = 0
        for i in range(self.n_chunks):
            this_size = chunk_size + (1 if i < remainder else 0)
            start = offset
            end = offset + this_size
            yield file_path, start, end
            offset = end

    def process_chunk(self, file_chunk_info):
        """
        Given (file_path, start, end), read those bytes, split them into lines,
        and treat every 4th line as a quality line.
        Return partial sums as a dict: pos->[sum, count].
        """
        file_path, start, end = file_chunk_info
        raw_data = self.read_chunk_bytes(file_path, start, end)

        lines = raw_data.split(b"\n")
        partial_sums = defaultdict(lambda: [0.0, 0.0])

        i = 0
        # the 4th line in each block is a quality line
        while i < len(lines) - 3:
            # lines[i+3] is the quality line
            quality_line = lines[i + 3].decode("utf-8", errors="ignore")
            scores = self.compute_phred_scores(quality_line)
            for idx, val in enumerate(scores):
                partial_sums[idx][0] += val
                partial_sums[idx][1] += 1
            i += 4

        return partial_sums


def output_phred_results(all_averages, output_prefix, run_index, workers):
    """
    Output the final average PHRED scores. If output_prefix is given,
    write each file's results to <prefix>_fileN_runX_wY.csv. Otherwise print to STDOUT.

    Args:
        all_averages (dict): { filename: { pos: average_score } }
        output_prefix (str): If None, print to STDOUT. Otherwise CSV to files.
        run_index (int): The run # for performance.
        workers (int): The reported worker count (for numeric stability).
    """
    if output_prefix:
        # One CSV per file
        for idx, (fname, posdict) in enumerate(all_averages.items(), start=1):
            out_name = f"{output_prefix}_file{idx}_run{run_index}_w{workers}.csv"
            with open(out_name, "w", newline="") as csvfile:
                writer = csv.writer(csvfile)
                writer.writerow(["Position", "AveragePHRED"])
                for pos in sorted(posdict.keys()):
                    writer.writerow([pos, posdict[pos]])
        print(f"# [Rank0] Wrote {len(all_averages)} CSV(s) with prefix='{output_prefix}'")
    else:
        # Print to STDOUT
        print("Position,AveragePHRED,Filename")
        for fname, posdict in all_averages.items():
            for pos in sorted(posdict.keys()):
                print(f"{pos},{posdict[pos]},{fname}")


def main():
    """
    Main function:
    - Parse arguments
    - MPI environment: rank0 = controller, rank>0 = workers
    - rank0 splits each file into 'chunks'
    - workers do partial sums, rank0 merges results => final averages
    """
    args = parse_args()
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    start_time = MPI.Wtime()

    calculator = PhredscoreCalculator(n_chunks=args.n)

    if rank == 0:
        # Controller

        tasks = []
        file_names = []
        for fpath in args.files:
            path = Path(fpath)
            for chunk_info in calculator.chunk_file(path):
                tasks.append(chunk_info)
                file_names.append(str(path))

        num_workers = size - 1
        if num_workers < 1:
            print("Error: Need at least 2 ranks (1 controller + 1 worker).")
            sys.exit(1)

        #  Distribute tasks to workers
        dist_data = [[] for _ in range(num_workers)]
        for i, t in enumerate(tasks):
            w_index = i % num_workers
            dist_data[w_index].append((t, file_names[i]))

        for w in range(1, size):
            comm.send(dist_data[w - 1], dest=w, tag=77)

        # Receive partial results
        combined = {}

        for _ in range(num_workers):
            worker_partials = comm.recv(source=MPI.ANY_SOURCE, tag=88)
            # Merge
            for fname, partial_dict in worker_partials.items():
                if fname not in combined:
                    combined[fname] = defaultdict(zero_pair)
                for pos, (s, c) in partial_dict.items():
                    combined[fname][pos][0] += s
                    combined[fname][pos][1] += c

        # Compute final averages
        final_averages = {}
        for fname, sc_map in combined.items():
            final_averages[fname] = {}
            for pos, (s, c) in sc_map.items():
                if c == 0:
                    final_averages[fname][pos] = 0.0
                else:
                    final_averages[fname][pos] = s / c

        # Output
        output_phred_results(
            final_averages,
            args.output_prefix,
            run_index=args.run_index,
            workers=args.workers
        )

        end_time = MPI.Wtime()
        elapsed = end_time - start_time
        print(f"# [Rank0] Finished in {elapsed:.4f} seconds", file=sys.stderr)

    else:
        # Worker
        my_tasks = comm.recv(source=0, tag=77)

        # Process partial sums
        partial_results = {}
        for (chunkinfo, fname) in my_tasks:
            sums_map = calculator.process_chunk(chunkinfo)
            if fname not in partial_results:
                partial_results[fname] = defaultdict(zero_pair)
            for pos, (s, c) in sums_map.items():
                partial_results[fname][pos][0] += s
                partial_results[fname][pos][1] += c

        # Send results
        comm.send(partial_results, dest=0, tag=88)


if __name__ == "__main__":
    main()
