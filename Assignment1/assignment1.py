import argparse
import multiprocessing
import csv


def compute_phred_scores(quality_str):
    return [ord(char) - 33 for char in quality_str]


def process_chunk(chunk):
    results = []
    for i, line in enumerate(chunk):
        if i % 4 == 3:  # Quality line
            phred_scores = compute_phred_scores(line.strip())
            if len(results) < len(phred_scores):
                results.extend([0] * (len(phred_scores) - len(results)))
            results = [x + y for x, y in zip(results, phred_scores)]
    return results


def main():
    parser = argparse.ArgumentParser(description="Script for Assignment 1 of Big Data Computing")
    parser.add_argument("-n", required=True, type=int, help="Number of cores to use.")
    parser.add_argument("-o", type=argparse.FileType('w', encoding='UTF-8'),
                        help="CSV file to save the output. Default is output to terminal STDOUT")
    parser.add_argument("fastq_files", type=argparse.FileType('r'), nargs='+',
                        help="At least one Illumina Fastq Format file to process")
    args = parser.parse_args()

    pool = multiprocessing.Pool(args.n)

    for fastq_file in args.fastq_files:
        chunk_size = 4  # Process four lines at a time
        file_results = []
        while True:
            chunk = [fastq_file.readline() for _ in range(chunk_size)]
            if not any(chunk):  # End of file
                break
            results = pool.apply_async(process_chunk, (chunk,))
            file_results.append(results.get())

        # Aggregate and calculate average PHRED scores across all chunks
        aggregated_scores = [sum(x) for x in zip(*file_results)]
        average_scores = [x / len(file_results) for x in aggregated_scores]

        # Output results
        if args.o:
            csv_writer = csv.writer(args.o)
            csv_writer.writerow(average_scores)
        else:
            print(','.join(map(str, average_scores)))


if __name__ == "__main__":
    main()

