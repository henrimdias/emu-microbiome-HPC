#!/usr/bin/env python3

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
