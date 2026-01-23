#!/usr/bin/env python3
# =============================================================================
# Article:
# Dias, H.M., et al. Reproducible Emu-based workflow for high-fidelity soil and
# plant microbiome profiling on HPC clusters. Bio-protocol. 2025.
#
# Script:
# Collect EMU mapping summary metrics from multiple .out files into one table.
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
# This script scans a directory containing EMU (or pipeline) log files (*.out),
# extracts key read-count metrics, and writes a single tab-separated summary file:
#
#   summary_counts.tsv
#
# Extracted metrics (one row per .out file):
#   - Barcode (e.g., barcode01, barcode12, etc.)
#   - Unmapped read count
#   - Mapped read count
#   - Unclassified mapped read count
#
# Assumptions:
#   - The directory contains one or more files ending in ".out"
#   - Each .out file includes lines like:
#       "Unmapped read count: <INT>"
#       "Mapped read count: <INT>"
#       "Unclassified mapped read count: <INT>"
#   - A barcode identifier appears somewhere in the file as "barcodeXX"
#     (e.g., barcode01, barcode12).
#
# Inputs (command-line arguments):
#   1) <directory_with_out_files> : path to folder containing *.out files
#
# Outputs:
#   - summary_counts.tsv : tab-separated table written in the current directory
#
# Usage:
#   python collect_counts.py /path/to/out_files/
#
# Notes:
#   - If a metric is not found in a given file, it is reported as "NA".
#   - If no barcode is found in a given file, Barcode is reported as "NA".
# =============================================================================

import glob
import re
import os
import sys


def main():
    # === Step 1: Parse and validate command-line arguments ===
    # Expect exactly one argument: the directory containing the .out files.
    if len(sys.argv) != 2:
        print(f"Uso: {sys.argv[0]} {{diretorio_com_out}}")
        sys.exit(1)

    input_directory = sys.argv[1]

    # Confirm the user provided a real directory path.
    if not os.path.isdir(input_directory):
        print(f"Erro: '{input_directory}' não é um diretório válido.")
        sys.exit(1)

    # === Step 2: Define output file name ===
    # Output is written to the current working directory (where the script is run).
    output_file = "summary_counts.tsv"

    # === Step 3: Compile regex patterns for metrics and barcode ===
    # These patterns search the log text and capture the integer values.
    unmapped_re = re.compile(r"Unmapped read count:\s*(\d+)")
    mapped_re = re.compile(r"Mapped read count:\s*(\d+)")
    unclassified_re = re.compile(r"Unclassified mapped read count:\s*(\d+)")

    # Barcode pattern captures only the digits after "barcode".
    # Example: "barcode01" -> captures "01"
    barcode_re = re.compile(r"barcode(\d+)")

    # === Step 4: Open output file and write header ===
    # Using a context manager ensures the file is properly closed even if errors occur.
    with open(output_file, "w") as out_f:
        out_f.write(
            "Barcode\tUnmapped_Read_Count\tMapped_Read_Count\tUnclassified_Mapped_Read_Count\n"
        )

        # === Step 5: Iterate over all *.out files in the input directory ===
        pattern = os.path.join(input_directory, "*.out")

        for filepath in glob.glob(pattern):
            # Read entire file content (these logs are typically small enough for this)
            with open(filepath, "r") as f:
                content = f.read()

            # === Step 6: Extract metrics using regex ===
            # Each search returns a match object or None if not found.
            unmapped = unmapped_re.search(content)
            mapped = mapped_re.search(content)
            unclassified = unclassified_re.search(content)
            barcode_match = barcode_re.search(content)

            # If match exists, take the captured group; otherwise label as NA.
            barcode = barcode_match.group(1) if barcode_match else "NA"
            unmapped_count = unmapped.group(1) if unmapped else "NA"
            mapped_count = mapped.group(1) if mapped else "NA"
            unclassified_count = unclassified.group(1) if unclassified else "NA"

            # === Step 7: Write one row per file ===
            out_f.write(
                f"{barcode}\t{unmapped_count}\t{mapped_count}\t{unclassified_count}\n"
            )

    # === Step 8: Print completion message ===
    print(f"Resumo gerado em: {output_file}")


if __name__ == "__main__":
    main()
