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

