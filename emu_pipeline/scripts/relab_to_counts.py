#!/usr/bin/env python3
# =============================================================================
# Article:
# Dias, H.M., et al. Reproducible Emu-based workflow for high-fidelity soil and
# plant microbiome profiling on HPC clusters. Bio-protocol. 2025.
#
# Script:
# Convert relative-abundance tables to raw count tables using per-barcode
# totals from EMU summary_counts.tsv.
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
# This script reconstructs integer raw counts from a relative-abundance table
# (e.g., EMU downstream output) by:
#   - Reading per-barcode totals from summary_counts.tsv.
#   - Matching sample columns to barcodes via numeric keys in their names.
#   - Multiplying relative abundances by the chosen per-barcode total
#     (mapped / mapped_classified / total).
#   - Optionally enforcing that integer counts per sample sum EXACTLY to the
#     chosen total using the Largest Remainder Method.
#
# Key features:
#   - Robust parsing of summary_counts.tsv (mapped, unmapped, unclassified).
#   - Heuristic detection of sample vs taxonomy/metadata columns.
#   - Automatic detection of relative scale (fractions 0–1 vs percentages 0–100).
#   - Strict-sum reconstruction to match per-sample totals if requested.
#
# Assumptions:
#   - summary_counts.tsv contains the columns:
#       Barcode, Unmapped_Read_Count, Mapped_Read_Count,
#       Unclassified_Mapped_Read_Count
#   - The relative-abundance table is TSV and has sample columns whose names
#     contain digits or start with "barcode"/"sample".
#
# Inputs (command-line arguments):
#   summary_counts       : path to summary_counts.tsv
#   relative_abundance   : path to relative-abundance table (TSV)
#   output_counts        : path to write raw-count table (TSV)
#
# Optional flags:
#   --total-mode {mapped, mapped_classified, total}
#       mapped            = Mapped_Read_Count (default)
#       mapped_classified = Mapped_Read_Count - Unclassified_Mapped_Read_Count
#       total             = Mapped_Read_Count + Unmapped_Read_Count
#   --strict-sum
#       Enforce that reconstructed counts per barcode sum exactly to the
#       chosen total using the Largest Remainder Method.
#
# Outputs:
#   - A TSV raw-count table with the same structure as the input relative
#     table, but with integer counts instead of relative values.
#
# Usage:
#   python relab_to_counts.py summary_counts.tsv rel_abundance.tsv counts.tsv \
#       --total-mode mapped_classified --strict-sum
#
# For full reproducibility, the versions of Python and dependencies are
# documented in the manuscript / accompanying documentation.
# =============================================================================

import argparse
import csv
import math
import os
import re
import statistics
from typing import List, Tuple, Dict

def parse_args():
    p = argparse.ArgumentParser(
        description="Convert a relative-abundance table to raw counts using per-barcode totals from summary_counts.tsv."
    )
    p.add_argument("summary_counts", help="Path to summary_counts.tsv")
    p.add_argument("relative_abundance", help="Path to relative abundance table (TSV)")
    p.add_argument("output_counts", help="Path to write raw-count table (TSV)")
    p.add_argument(
        "--total-mode",
        choices=["mapped", "mapped_classified", "total"],
        default="mapped",
        help=(
            "Per-barcode total (denominator):\n"
            "  mapped            = Mapped_Read_Count (default)\n"
            "  mapped_classified = Mapped_Read_Count - Unclassified_Mapped_Read_Count\n"
            "  total             = Mapped_Read_Count + Unmapped_Read_Count"
        ),
    )
    p.add_argument(
        "--strict-sum",
        action="store_true",
        help="Make integer counts per barcode sum EXACTLY to the chosen total (Largest Remainder Method)."
    )
    return p.parse_args()

def canonical_barcode_key(x: str) -> str:
    m = re.search(r"(\d+)", x)
    if not m:
        return None
    return f"{int(m.group(1)):02d}"

def read_summary_counts(path: str, total_mode: str) -> Dict[str, int]:
    totals = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        required = {"Barcode", "Unmapped_Read_Count", "Mapped_Read_Count", "Unclassified_Mapped_Read_Count"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f"summary_counts.tsv missing required columns. Found: {reader.fieldnames}")
        for row in reader:
            key = canonical_barcode_key(row["Barcode"])
            if key is None:
                continue
            unmapped = 0 if row["Unmapped_Read_Count"] in ("", "NA") else int(row["Unmapped_Read_Count"])
            mapped = 0 if row["Mapped_Read_Count"] in ("", "NA") else int(row["Mapped_Read_Count"])
            unclassified = 0 if row["Unclassified_Mapped_Read_Count"] in ("", "NA") else int(row["Unclassified_Mapped_Read_Count"])
            if total_mode == "mapped":
                total = mapped
            elif total_mode == "mapped_classified":
                total = mapped - unclassified
            elif total_mode == "total":
                total = mapped + unmapped
            else:
                raise ValueError("Invalid total_mode")
            totals[key] = max(0, total)
    if not totals:
        raise ValueError("No barcode totals parsed from summary_counts.tsv")
    return totals

def read_tsv(path: str) -> Tuple[List[str], List[List[str]]]:
    with open(path, newline="") as f:
        rows = list(csv.reader(f, delimiter="\t"))
    if not rows:
        raise ValueError(f"Empty file: {path}")
    return rows[0], rows[1:]

def write_tsv(path: str, header: List[str], rows: List[List[str]]) -> None:
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(header)
        w.writerows(rows)

def detect_columns(header: List[str]) -> Tuple[List[str], List[str]]:
    """
    Heuristic: sample columns are those whose names contain digits OR start with 'barcode' or 'sample' (case-insensitive).
    Everything else is metadata/taxonomy columns.
    """
    sample_cols, tax_cols = [], []
    for col in header:
        low = col.lower()
        if low.startswith("barcode") or low.startswith("sample") or re.search(r"\d", col):
            sample_cols.append(col)
        else:
            tax_cols.append(col)
    if not sample_cols and len(header) > 1:
        # Fallback: first column taxonomy, remaining are samples
        tax_cols, sample_cols = [header[0]], header[1:]
    return tax_cols, sample_cols

def collect_positive_values(rel_rows: List[List[str]], header: List[str], sample_cols: List[str], limit: int = 2000) -> List[float]:
    idx_map = {c: header.index(c) for c in sample_cols}
    vals = []
    for r in rel_rows[:limit]:
        for c in sample_cols:
            s = r[idx_map[c]].strip()
            if s and s.upper() != "NA":
                try:
                    v = float(s)
                    if v > 0:
                        vals.append(v)
                except Exception:
                    pass
    return vals

def largest_remainder_adjust(int_counts: List[int], float_counts: List[float], target_sum: int) -> List[int]:
    current_sum = sum(int_counts)
    diff = target_sum - current_sum
    if diff == 0 or len(int_counts) == 0:
        return int_counts
    remainders = [(fc - math.floor(fc), idx) for idx, fc in enumerate(float_counts)]
    if diff > 0:
        remainders.sort(reverse=True)
        for k in range(min(diff, len(int_counts))):
            int_counts[remainders[k][1]] += 1
    else:
        remainders.sort()
        removed = 0
        for _, idx in remainders:
            if int_counts[idx] > 0:
                int_counts[idx] -= 1
                removed += 1
                if removed == -diff:
                    break
    return int_counts

def relab_to_counts(summary_counts_path: str, relab_path: str, out_path: str,
                    total_mode: str = "mapped", strict_sum: bool = True) -> None:
    totals = read_summary_counts(summary_counts_path, total_mode=total_mode)
    rel_header, rel_rows = read_tsv(relab_path)
    tax_cols, sample_cols = detect_columns(rel_header)

    # Map samples to barcode keys (by the digits in the column name)
    barcode_keys = {col: canonical_barcode_key(col) for col in sample_cols}
    missing = [col for col, key in barcode_keys.items() if key not in totals or key is None]
    if missing:
        raise ValueError("These sample columns are missing in summary_counts (by numeric key): " + ", ".join(missing))

    per_barcode_total = {col: totals[barcode_keys[col]] for col in sample_cols}

    # Detect if relative values are 0..1 or 0..100 (percentages)
    positives = collect_positive_values(rel_rows, rel_header, sample_cols, limit=2000)
    is_percent = (statistics.median(positives) > 1.0) if positives else False

    out_header = tax_cols + sample_cols
    out_rows: List[List[str]] = []

    idx = {c: rel_header.index(c) for c in rel_header}

    # For strict-sum adjustment
    cols_float = {c: [] for c in sample_cols}
    cols_int   = {c: [] for c in sample_cols}

    for r in rel_rows:
        tax_part = [r[idx[c]] for c in tax_cols]

        cur_float_counts, cur_int_counts = [], []
        for c in sample_cols:
            s = r[idx[c]].strip()
            rel_v = 0.0
            if s and s.upper() != "NA":
                try:
                    rel_v = float(s)
                except Exception:
                    rel_v = 0.0
            if is_percent:
                rel_v /= 100.0
            total = per_barcode_total[c]
            fc = rel_v * total
            cur_float_counts.append(fc)
            cur_int_counts.append(math.floor(fc))

        for j, c in enumerate(sample_cols):
            cols_float[c].append(cur_float_counts[j])
            cols_int[c].append(cur_int_counts[j])

        out_rows.append(tax_part + [str(v) for v in cur_int_counts])

    if strict_sum:
        for c in sample_cols:
            cols_int[c] = largest_remainder_adjust(cols_int[c][:], cols_float[c], per_barcode_total[c])
        for r_i in range(len(out_rows)):
            for c_i, c in enumerate(sample_cols):
                out_rows[r_i][len(tax_cols) + c_i] = str(cols_int[c][r_i])

    write_tsv(out_path, out_header, out_rows)

def main():
    args = parse_args()
    relab_to_counts(args.summary_counts, args.relative_abundance, args.output_counts,
                    total_mode=args.total_mode, strict_sum=args.strict_sum)
    print(f"Wrote raw-count table to: {args.output_counts}")
    print(f"Total mode: {args.total_mode}")

if __name__ == "__main__":
    main()

