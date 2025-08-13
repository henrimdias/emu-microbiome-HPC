#!/bin/bash

# Define output directory
OUTDIR="/home/jacks.local/henrique.dias/scratch/emu_pipeline/kraken2_db"
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
