#!/bin/bash
# =============================================================================
# Article:
# Dias, H.M., et al. Reproducible Emu-based workflow for high-fidelity soil and
# plant microbiome profiling on HPC clusters. Bio-protocol. 2025.
#
# Script:
# Download and prepare RefSeq plastid and mitochondrion genomes for building
# an organelle-focused Kraken2 database within the workflow described above.
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
# This helper script downloads RefSeq plastid and mitochondrion genome FASTA
# files from NCBI, stores them under a user-defined directory, and decompresses
# them so they can be used as input for Kraken2 database construction.
#
# The script performs:
#   - Creation of the output directory (if it does not exist)
#   - Download of plastid and mitochondrion RefSeq genome archives (.fna.gz)
#   - Basic integrity check that both files were downloaded
#   - Decompression of the .gz archives into plain FASTA (.fna) files
#
# Assumptions:
#   - wget and gunzip are installed and available in the PATH
#   - The user has write permissions to the chosen OUTDIR
#
# Inputs (user configuration below):
#   OUTDIR      : path to the directory where organelle FASTA files will be saved
#   PLASTID_URL : URL to RefSeq plastid genomes (.fna.gz)
#   MITO_URL    : URL to RefSeq mitochondrion genomes (.fna.gz)
#
# Outputs:
#   - plastid.3.1.genomic.fna
#   - mitochondrion.1.1.genomic.fna
#   saved under $OUTDIR, ready for inclusion in a Kraken2 database build.
#
# Usage:
#   bash download_plant_organelle_RefSeq_fastas.sh
#
# For full reproducibility, the software versions and URLs used here are
# documented in the manuscript / accompanying documentation.
# =============================================================================

# Define output directory
OUTDIR="/put/your/path/emu_pipeline/kraken2_db"
mkdir -p "$OUTDIR"

echo "📁 Organellar FASTA output will be saved to: $OUTDIR"
echo "🔽 Starting download of RefSeq plastid and mitochondrion genome files..."

# File URLs
PLASTID_URL="https://ftp.ncbi.nlm.nih.gov/genomes/refseq/plastid/plastid.3.1.genomic.fna.gz"
MITO_URL="https://ftp.ncbi.nlm.nih.gov/genomes/refseq/mitochondrion/mitochondrion.1.1.genomic.fna.gz"

# Download files
wget -nc -P "$OUTDIR" "$PLASTID_URL"
wget -nc -P "$OUTDIR" "$MITO_URL"

# Check if files were downloaded
if [[ -f "$OUTDIR/plastid.3.1.genomic.fna.gz" && -f "$OUTDIR/mitochondrion.1.1.genomic.fna.gz" ]]; then
  echo "✅ Files downloaded successfully."
else
  echo "❌ Error: One or both files were not downloaded correctly."
  exit 1
fi

# Decompress
echo "🧬 Decompressing FASTA files..."
gunzip -f "$OUTDIR/plastid.3.1.genomic.fna.gz"
gunzip -f "$OUTDIR/mitochondrion.1.1.genomic.fna.gz"

# Confirm decompression
if [[ -f "$OUTDIR/plastid.3.1.genomic.fna" && -f "$OUTDIR/mitochondrion.1.1.genomic.fna" ]]; then
  echo "✅ FASTA files successfully decompressed."
  echo "🎉 Organelle reference files ready for Kraken2 database build."
else
  echo "❌ Error: Decompression failed."
  exit 1
fi
