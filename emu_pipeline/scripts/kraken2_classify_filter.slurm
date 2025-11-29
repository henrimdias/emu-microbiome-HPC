#!/bin/bash
# =============================================================================
# Article:
# Dias, H.M., et al. Reproducible Emu-based workflow for high-fidelity soil and
# plant microbiome profiling on HPC clusters. Bio-protocol. 2025.
#
# Script:
# Classify reads with an organelle-focused Kraken2 database and remove
# organelle-matching reads from ONT FASTQ files prior to EMU profiling.
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
# This SLURM batch script:
#   - Runs Kraken2 classification on demultiplexed FASTQ files using an
#     organelle-targeted database (plastid + mitochondrion).
#   - Generates Kraken2 reports and classification outputs for each sample.
#   - Calls a companion Python script to remove organelle-classified reads,
#     producing organelle-depleted FASTQ files for downstream EMU analysis.
#
# The script performs:
#   - Iteration over all *.fastq.gz files in READ_DIR
#   - Kraken2 classification with user-defined confidence threshold
#   - Organelle-read filtering via remove_kraken2_organelle_reads.py
#
# Assumptions:
#   - Miniconda is installed and accessible at ~/miniconda3
#   - Conda environment `preemu_env` contains Kraken2 and Python dependencies
#   - DBDIR points to a built Kraken2 organelle database
#   - READ_DIR contains gzip-compressed FASTQ files (*.fastq.gz)
#   - The Python script remove_kraken2_organelle_reads.py is available and
#     executable at the path specified in SCRIPT.
#
# Inputs (user configuration below):
#   READ_DIR : directory with input FASTQ files
#   DBDIR    : path to Kraken2 organelle database
#   OUTDIR   : directory for reports, Kraken2 outputs, and filtered FASTQ files
#   SCRIPT   : path to remove_kraken2_organelle_reads.py
#
# Outputs (per sample):
#   - ${sample}_report.txt         : Kraken2 classification report
#   - ${sample}_kraken_output.txt  : Kraken2 per-read classification output
#   - ${sample}_filtered.fastq     : organelle-depleted reads for EMU
#
# Usage:
#   sbatch kraken2_classify_filte.slurm
#
# For full reproducibility, the versions of Kraken2, Python, and dependencies
# are documented in the manuscript / accompanying documentation.
# =============================================================================

#SBATCH --job-name=kraken2_filter_reads
#SBATCH --output=/put/your/path/emu_pipeline/no_organelle_filtered_data/kraken2_filter_%j.out
#SBATCH --error=/put/your/path/emu_pipeline/no_organelle_filtered_data/kraken2_filter_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --partition=compute
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=put-your-email@email.edu

# Load environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate preemu_env

# Define input/output
READ_DIR="/put/your/path/emu_pipeline/filtered_data"
DBDIR="/put/your/path/emu_pipeline/kraken2_db/db"
OUTDIR="/put/your/path/emu_pipeline/no_organelle_filtered_data"
SCRIPT="/put/your/path/emu_pipeline/scripts/remove_kraken2_organelle_reads.py"

mkdir -p "$OUTDIR"

# Loop over all FASTQ files
for fq in "$READ_DIR"/*.fastq.gz; do
  sample=$(basename "$fq" .fastq.gz)

  echo "🧬 Processing $sample..."

  kraken2 \
    --db "$DBDIR" \
    --threads 8 \
    --confidence 0.50 \
    --report "$OUTDIR/${sample}_report.txt" \
    --output "$OUTDIR/${sample}_kraken_output.txt" \
    --gzip-compressed \
    "$fq"

  echo "🧹 Filtering organelle reads for $sample..."

  python "$SCRIPT" \
    --input "$fq" \
    --kraken_output "$OUTDIR/${sample}_kraken_output.txt" \
    --output "$OUTDIR/${sample}_filtered.fastq"
done

echo "✅ All samples processed."
