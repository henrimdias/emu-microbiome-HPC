!/usr/bin/env python
# =============================================================================
# Article:
# Dias, H.M., et al. Reproducible Emu-based workflow for high-fidelity soil and
# plant microbiome profiling on HPC clusters. Bio-protocol. 2025.
#
# Script:
# Convert EMU taxonomic output into a FAPROTAX-compatible abundance table.
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
# This script takes an EMU taxonomic abundance table (TSV) and:
#   - Concatenates EMU taxonomic ranks (superkingdom to genus) into a single
#     semicolon-separated "taxonomy" string compatible with FAPROTAX.
#   - Retains only the taxonomy column plus sample abundance columns.
#   - Writes a FAPROTAX-ready TSV file.
#
# Assumptions:
#   - The input table is a TSV file produced by EMU (or post-processed) and
#     contains the columns:
#       superkingdom, phylum, class, order, family, genus
#   - Sample columns are named with the prefix "barcode" (e.g., barcode01, etc.).
#
# Inputs (command-line arguments):
#   1) input_file.tsv  : EMU-style taxonomic table (TSV)
#   2) output_file.tsv : path to write FAPROTAX-ready table (TSV)
#
# Outputs:
#   - A TSV file with:
#       taxonomy    (semicolon-separated lineage)
#       barcode...  (one column per sample)
#
# Usage:
#   python emu_to_faprotax.py emu_taxonomy.tsv faprotax_input.tsv
#
# For full reproducibility, the versions of Python, pandas, and EMU used
# are documented in the manuscript / accompanying documentation.
# =============================================================================

import pandas as pd
import sys

# === Step 1: Parse command-line arguments ===
if len(sys.argv) != 3:
    print("Usage: python emu_to_faprotax.py <input_file.tsv> <output_file.tsv>")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

# === Step 2: Load input file ===
try:
    df = pd.read_csv(input_file, sep="\t")
except Exception as e:
    print(f"❌ Error reading input file: {e}")
    sys.exit(1)

# === Step 3: Build FAPROTAX-style taxonomy ===
df["taxonomy"] = (
    df["superkingdom"].fillna("")
    + ";" + df["phylum"].fillna("")
    + ";" + df["class"].fillna("")
    + ";" + df["order"].fillna("")
    + ";" + df["family"].fillna("")
    + ";" + df["genus"].fillna("")
)

# === Step 4: Keep only taxonomy + sample columns ===
sample_cols = [col for col in df.columns if col.startswith("barcode")]
if not sample_cols:
    print("❌ No sample columns found starting with 'barcode'")
    sys.exit(1)

output_df = df[["taxonomy"] + sample_cols]

# === Step 5: Write output file ===
try:
    output_df.to_csv(output_file, sep="\t", index=False)
    print(f"✅ FAPROTAX-ready file saved: {output_file}")
except Exception as e:
    print(f"❌ Error writing output file: {e}")
    sys.exit(1)

