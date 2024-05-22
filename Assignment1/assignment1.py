import argparse
import multiprocessing
import csv
import time
import os


def compute_phred_scores(quality_str):
    """Compute the PHRED scores from a quality string."""
    return [ord(char) - 33 for char in quality_str]


def process_chunk(lines):
    """Process a chunk of FastQ file lines to calculate cumulative PHRED scores."""
    results = []
    for line in lines:
        if line:  # Check for non-empty line
            phred_scores = compute_phred_scores(line.strip())
            if len(results) < len(phred_scores):
                results.extend([0] * (len(phred_scores) - len(results)))
            results = [x + y for x, y in zip(results, phred_scores)]
    return results


def main():
    """Main"""
    start_time = time.time()

    parser = argparse.ArgumentParser(
        description="Calculate average PHRED scores from FastQ files."
    )
    parser.add_argument("-n", required=True, type=int, help="Number of cores to use.")
    parser.add_argument(
        "-o", action="store_true", help="Output results to individual CSV files"
    )
    parser.add_argument("fastq_files", nargs="+", help="FastQ files to process")
    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))

    with multiprocessing.Pool(args.n) as pool:
        for fastq_path in args.fastq_files:
            with open(fastq_path, "r") as fastq_file:
                read_data = [line.strip() for line in fastq_file if line.strip()]
                chunks = [read_data[i : i + 4] for i in range(0, len(read_data), 4)]
                results = pool.map(process_chunk, chunks)
                aggregated_scores = [sum(scores) for scores in zip(*results)]
                average_scores = [x / len(results) for x in aggregated_scores]

                formatted_scores = [(i, score) for i, score in enumerate(average_scores)]

                if args.o:
                    output_file_name = os.path.join(script_dir, f"{os.path.basename(fastq_path)}.output.csv")
                    with open(output_file_name, "w", newline="") as csvfile:
                        writer = csv.writer(csvfile)
                        writer.writerows(formatted_scores)
                    print(f"Output for {fastq_path} written to {output_file_name}")
                else:
                    for index, score in formatted_scores:
                        print(f"{index},{score}")

    end_time = time.time()  # End timing
    print(f"Execution time: {end_time - start_time:.2f} seconds")  # Print the execution time


if __name__ == "__main__":
    main()