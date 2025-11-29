#!/usr/bin/env python3
# =============================================================================
# Article:
# Dias, H.M., et al. Reproducible Emu-based workflow for high-fidelity soil and
# plant microbiome profiling on HPC clusters. Bio-protocol. 2025.
#
# Script:
# Summarize Kraken2-based organelle contamination (chloroplast, mitochondria)
# and bacterial/unclassified fractions across samples.
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
# This script parses multiple Kraken2 report files (one per sample) and
# summarizes the percentage of:
#   - chloroplast-assigned reads
#   - mitochondrion-assigned reads
#   - unclassified reads
#   - bacterial reads
# into a single TSV table.
#
# The script:
#   - Scans an input directory for files matching *_report.txt.
#   - Extracts percentage values from Kraken2 report columns.
#   - Aggregates percentages for categories based on taxon names.
#   - Produces a tab-delimited summary table sorted by sample name.
#
# Assumptions:
#   - Kraken2 report format is standard, with:
#       column 1 = percentage of reads
#       last column = taxon name (e.g., "Bacteria", "unclassified",
#                     "... chloroplast", "... mitochondrion").
#   - Organelle contamination is identified by the presence of the substrings
#     "chloroplast", "mitochondrion"/"mitochondria" in the taxon name.
#
# Inputs (command-line arguments):
#   --input_dir    : directory containing *_report.txt Kraken2 reports
#   --output_file  : path to write the summary TSV table
#
# Outputs:
#   - TSV file with columns:
#       Sample
#       Chloroplast (%)
#       Mitochondria (%)
#       Unclassified (%)
#       Bacterial (%)
#
# Usage:
#   python summarize_kraken2_reports.py \
#       --input_dir /path/to/reports \
#       --output_file kraken_organelle_summary.tsv
#
# For full reproducibility, the versions of Python, pandas, and Kraken2
# are documented in the manuscript / accompanying documentation.
# =============================================================================

import os
import argparse
import glob
import pandas as pd

def parse_report(report_path):
    chloroplast = 0.0
    mitochondria = 0.0
    unclassified = 0.0
    bacterial = 0.0

    with open(report_path, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 6:
                continue
            perc = float(parts[0])
            name = parts[-1].lower()

            if "chloroplast" in name:
                chloroplast += perc
            elif "mitochondrion" in name or "mitochondria" in name:
                mitochondria += perc
            elif name == "unclassified":
                unclassified += perc
            elif "bacteria" in name:
                bacterial += perc

    return chloroplast, mitochondria, unclassified, bacterial

def summarize_reports(input_dir, output_file):
    summary = []

    for report in glob.glob(os.path.join(input_dir, "*_report.txt")):
        sample = os.path.basename(report).replace("_report.txt", "")
        chl, mito, unclass, bact = parse_report(report)
        summary.append({
            "Sample": sample,
            "Chloroplast (%)": round(chl, 2),
            "Mitochondria (%)": round(mito, 2),
            "Unclassified (%)": round(unclass, 2),
            "Bacterial (%)": round(bact, 2)
        })

    df = pd.DataFrame(summary)
    df = df.sort_values(by="Sample")
    df.to_csv(output_file, sep="\t", index=False)
    print(f"✅ Summary saved to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Summarize Kraken2 organelle contamination from multiple report files.")
    parser.add_argument("--input_dir", required=True, help="Directory with *_report.txt files")
    parser.add_argument("--output_file", required=True, help="Output summary table (TSV)")

    args = parser.parse_args()
    summarize_reports(args.input_dir, args.output_file)
