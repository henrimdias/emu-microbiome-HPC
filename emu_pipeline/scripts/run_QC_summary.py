#!/usr/bin/env python3
# =============================================================================
# Article:
# Dias, H.M., et al. Reproducible Emu-based workflow for high-fidelity soil and
# plant microbiome profiling on HPC clusters. Bio-protocol. 2025.
#
# Script:
# Summarize NanoStat reports (pre- and post-NanoFilt) into CSV tables and
# multi-page PDF reports with tables and QC plots.
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
# This script parses NanoStat text reports generated before and after read
# filtering (e.g., NanoFilt), aggregates key sequencing metrics per barcode,
# and exports:
#   - A master CSV with all barcodes × phases (pre/post).
#   - Per-barcode CSV summaries.
#   - A multi-page PDF with tabular summaries.
#   - A multi-page PDF with bar plots for each metric across barcodes.
#
# The script:
#   - Recursively scans two directories (pre- and post-filter) for NanoStat
#     report files.
#   - Extracts metrics using regular expressions (mean/median length, quality,
#     number of reads, N50, etc.).
#   - Filters out “unclassified” folders by name.
#   - Builds a tidy pandas DataFrame indexed by (barcode, phase).
#   - Generates A4-size PDFs for both tables and plots.
#
# Assumptions:
#   - NanoStat reports are plain text files matching patterns:
#       * pre_dir:  *NanoStats.txt
#       * post_dir: *filteredNanoStats.txt
#   - Folder names contain a barcode suffix like `_barcodeXX` to identify
#     barcodes; otherwise, the folder name is used as the barcode label.
#   - matplotlib and pandas are installed and available.
#
# Inputs (user configuration below):
#   pre_dir    : directory containing pre-filter NanoStat reports
#   post_dir   : directory containing post-filter NanoStat reports
#   output_dir : directory where CSV and PDF summaries will be written
#
# Outputs:
#   - sequencing_summary.csv         : master summary for all barcodes/phases
#   - {barcode}_summary.csv          : per-barcode summaries
#   - sequencing_table.pdf           : multi-page table summary
#   - sequencing_plots.pdf           : metric-wise bar plots (pre vs post)
#
# Usage:
#   python run_QC_summary.py
#
# For full reproducibility, the software versions (Python, pandas, matplotlib,
# NanoStat, NanoFilt) are documented in the manuscript / documentation.
# =============================================================================

import os
import re
import glob
import math
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import MaxNLocator

# 1. Adjust these paths to your directories:
pre_dir     = "/home/jacks.local/henrique.dias/scratch/emu_pipeline/nanoplot_reports"
post_dir    = "/home/jacks.local/henrique.dias/scratch/emu_pipeline/nanoplot_after_nanofilt_data"
output_dir  = "/home/jacks.local/henrique.dias/scratch/emu_pipeline/qc_summary_reports"
os.makedirs(output_dir, exist_ok=True)

# 2. Define the metrics and their regex patterns:
metrics = {
    "Mean read length":     r"Mean read length:\s*([\d,\.]+)",
    "Mean read quality":    r"Mean read quality:\s*([\d,\.]+)",
    "Median read length":   r"Median read length:\s*([\d,\.]+)",
    "Median read quality":  r"Median read quality:\s*([\d,\.]+)",
    "Number of reads":      r"Number of reads:\s*([\d,\.]+)",
    "Read length N50":      r"Read length N50:\s*([\d,\.]+)",
    "STDEV read length":    r"STDEV read length:\s*([\d,\.]+)",
    "Q10_count":            r">Q10:\s*([\d,]+)",
    "Q10_pct":              r">Q10:.*\(([\d\.]+)%\)",
}

def parse_report(filepath):
    """Read a NanoStat report and return a dict metric → value."""
    text = open(filepath, encoding="utf8").read()
    values = {}
    for name, pattern in metrics.items():
        m = re.search(pattern, text)
        values[name] = float(m.group(1).replace(",", "")) if m else float("nan")
    return values

records = []
# 3. Iterate over pre- and post-filter directories
for base_dir, phase, pattern in [
    (pre_dir,  "pre",  "*NanoStats.txt"),
    (post_dir, "post", "*filteredNanoStats.txt")
]:
    print(f"\nProcessing phase '{phase}' in {base_dir}")
    for folder_name in sorted(os.listdir(base_dir)):
        # skip any “unclassified” folders
        if "_unclassified" in folder_name.lower():
            print(f"  • Skipping unclassified folder: {folder_name}")
            continue

        folder_path = os.path.join(base_dir, folder_name)
        if not os.path.isdir(folder_path):
            continue

        matches = glob.glob(os.path.join(folder_path, pattern))
        if not matches:
            continue

        report_path = matches[0]
        # Extract barcode ID from folder name
        m = re.search(r"_barcode(\d+)", folder_name)
        barcode = m.group(1) if m else folder_name

        print(f"  ✓ Found report for barcode {barcode}: {os.path.basename(report_path)}")
        data = parse_report(report_path)
        data["barcode"] = barcode
        data["phase"]   = phase
        records.append(data)

if not records:
    raise RuntimeError("No valid reports found in any folder.")

# 4. Build DataFrame
df = pd.DataFrame(records).set_index(["barcode", "phase"]).sort_index()

# 5. Save master CSV and per-barcode CSVs
master_csv = os.path.join(output_dir, "sequencing_summary.csv")
df.to_csv(master_csv)
print(f"\n✔ Master summary CSV saved to: {master_csv}")

for bc, grp in df.groupby(level=0):
    csv_path = os.path.join(output_dir, f"{bc}_summary.csv")
    grp.to_csv(csv_path)
    print(f"  • Summary CSV for barcode {bc} saved to: {csv_path}")

# 6. Create a multi-page table PDF
table_pdf = os.path.join(output_dir, "sequencing_table.pdf")
rows_per_page = 25
table_df = df.reset_index()
num_pages = math.ceil(len(table_df) / rows_per_page)

with PdfPages(table_pdf) as pdf:
    for page in range(num_pages):
        start = page * rows_per_page
        end   = start + rows_per_page
        chunk = table_df.iloc[start:end]

        fig, ax = plt.subplots(figsize=(11.7, 8.3))  # A4 landscape in inches
        ax.axis("off")
        tbl = ax.table(
            cellText=chunk.values,
            colLabels=chunk.columns,
            loc="center",
            cellLoc="center"
        )
        tbl.auto_set_font_size(False)
        tbl.set_fontsize(8)
        tbl.scale(1, 1.2)
        ax.set_title(f"Metrics summary (rows {start+1}–{min(end, len(table_df))})", pad=20)
        pdf.savefig(fig)
        plt.close(fig)

print(f"\n✔ Table PDF saved to: {table_pdf} ({num_pages} pages)")

# 7. Create a separate PDF for all plots, each at A4 size with more y-axis breaks
plots_pdf = os.path.join(output_dir, "sequencing_plots.pdf")
with PdfPages(plots_pdf) as pdf:
    for metric_name in metrics.keys():
        pivot = df[metric_name].unstack("phase")
        fig, ax = plt.subplots(figsize=(11.7, 8.3))  # A4 landscape
        pivot.plot.bar(ax=ax, width=0.6)

        ax.set_title(metric_name)
        ax.set_ylabel(metric_name)
        ax.set_xticks(range(len(pivot.index)))
        ax.set_xticklabels(pivot.index, rotation=90, fontsize=8)
        ax.yaxis.set_major_locator(MaxNLocator(nbins=10))

        plt.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

print(f"✔ Plots PDF saved to: {plots_pdf}")
print(f"\nAll output files have been saved to: {output_dir}")
