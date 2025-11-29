#!/usr/bin/env python3
# =============================================================================
# Article:
# Dias, H.M., et al. Reproducible Emu-based workflow for high-fidelity soil and
# plant microbiome profiling on HPC clusters. Bio-protocol. 2025.
#
# Script:
# Filter FASTQ reads to keep only Kraken2-unclassified reads (non-organelle)
# prior to EMU-based taxonomic profiling.
#
# Author (script):
# Henrique M. Dias
#
# Affiliation:
# South Dakota State University
#
# Date:
# 2025
#
# Description:
# This script uses a Kraken2 classification output file to identify reads
# labeled as "Unclassified" (flag "U") and retains only those reads from a
# gzip-compressed FASTQ file. It is intended for removing organelle-matching
# reads (classified) and keeping non-organelle reads for downstream EMU
# analysis.
#
# The script performs:
#   - Parses Kraken2 output to collect IDs of unclassified reads.
#   - Iterates through a .fastq.gz file with Biopython SeqIO.
#   - Writes only unclassified reads to an output FASTQ file.
#   - Prints a simple summary of total, kept, and removed reads.
#
# Assumptions:
#   - Kraken2 output is in standard tab-delimited format with "U" or "C" in
#     the first column and read IDs in the second column.
#   - Input reads are gzip-compressed FASTQ (.fastq.gz).
#   - Biopython is installed and available (SeqIO).
#
# Inputs (command-line arguments):
#   --input          : input .fastq.gz file
#   --kraken_output  : Kraken2 output file (per-read classification)
#   --output         : output FASTQ file with only unclassified reads
#
# Outputs:
#   - A FASTQ file containing only unclassified reads, suitable for EMU
#     full-length 16S profiling after organelle contamination removal.
#
# Usage:
#   python remove_kraken2_organelle_reads.py \
#       --input sample.fastq.gz \
#       --kraken_output sample_kraken_output.txt \
#       --output sample_filtered.fastq
#
# For full reproducibility, the versions of Python, Biopython, and Kraken2
# are documented in the manuscript / accompanying documentation.
# =============================================================================

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
