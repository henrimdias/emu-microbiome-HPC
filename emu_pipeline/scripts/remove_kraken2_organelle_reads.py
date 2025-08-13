#!/usr/bin/env python3

import gzip
import argparse
from Bio import SeqIO

def parse_kraken_output(kraken_file):
    unclassified_ids = set()
    with open(kraken_file, "r") as f:
        for line in f:
            cols = line.strip().split("\t")
            if cols[0] == "U":  # Unclassified
                unclassified_ids.add(cols[1])
    return unclassified_ids

def filter_fastq(input_fastq, kraken_output, output_fastq):
    unclassified_reads = parse_kraken_output(kraken_output)
    total = 0
    kept = 0

    with gzip.open(input_fastq, "rt") as infile, open(output_fastq, "w") as outfile:
        for record in SeqIO.parse(infile, "fastq"):
            total += 1
            if record.id in unclassified_reads:
                kept += 1
                SeqIO.write(record, outfile, "fastq")

    print(f"Total reads: {total}")
    print(f"Kept {kept} unclassified reads in {output_fastq}")
    print(f"Removed {total - kept} classified reads")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Keep only unclassified reads from FASTQ based on Kraken2 output")
    parser.add_argument("--input", required=True, help="Input .fastq.gz file")
    parser.add_argument("--kraken_output", required=True, help="Kraken2 output file")
    parser.add_argument("--output", required=True, help="Output filtered FASTQ file")

    args = parser.parse_args()
    filter_fastq(args.input, args.kraken_output, args.output)
