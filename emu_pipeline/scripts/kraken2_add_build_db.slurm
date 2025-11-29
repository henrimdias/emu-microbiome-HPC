!/bin/bash
# =============================================================================
# Article:
# Dias, H.M., et al. Reproducible Emu-based workflow for high-fidelity soil and
# plant microbiome profiling on HPC clusters. Bio-protocol. 2025.
#
# Script:
# Final build of an organelle-focused Kraken2 database (plastid + mitochondrion)
# for use in the pre-EMU contamination filtering workflow.
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
# This SLURM batch script takes pre-downloaded plastid and mitochondrion RefSeq
# FASTA files, adds them to a Kraken2 database, downloads the NCBI taxonomy,
# and performs the final Kraken2 database build.
#
# The script performs:
#   - Environment setup via a Conda environment (preemu_env)
#   - Addition of plastid and mitochondrion FASTA files to the Kraken2 library
#   - Download of NCBI taxonomy required by Kraken2
#   - Final Kraken2 database build using multiple threads
#
# Assumptions:
#   - Miniconda is installed and accessible at ~/miniconda3
#   - The Conda environment `preemu_env` has Kraken2 installed
#   - DBDIR points to a valid location for the Kraken2 database
#   - FASTA_DIR contains decompressed plastid and mitochondrion .fna files
#
# Inputs (user configuration below):
#   DBDIR     : path to the Kraken2 database directory
#   FASTA_DIR : directory containing plastid.3.1.genomic.fna and
#               mitochondrion.1.1.genomic.fna
#
# Outputs:
#   - Fully built Kraken2 database under $DBDIR, ready for plastid/mitochondrial
#     read filtering prior to EMU-based microbiome profiling.
#
# Usage:
#   sbatch kraken2_build_final.sh
#
# For full reproducibility, the versions of Kraken2 and dependencies are listed
# in the manuscript / accompanying documentation.
# =============================================================================

#SBATCH --job-name=kraken2_build_final
#SBATCH --output=/put/your/path/kraken_organelle_db/kraken2_build_%j.out
#SBATCH --error=/put/your/path/henrique.dias/kraken_organelle_db/kraken2_build_%j.err
#SBATCH --time=2-00:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --partition=compute
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=put-your-email@email.edu

# Load Conda
source ~/miniconda3/etc/profile.d/conda.sh
conda activate preemu_env

# Set environment variables
export OMP_NUM_THREADS=8

# Define paths
DBDIR="/put/your/path/emu_pipeline/kraken2_db/db"
FASTA_DIR="/put/your/path/kraken_env/data/organelle_refseq"
PLASTID="$FASTA_DIR/plastid.3.1.genomic.fna"
MITO="$FASTA_DIR/mitochondrion.1.1.genomic.fna"

# Optional: clean partial build (uncomment if re-running)
# rm -rf "$DBDIR/library" "$DBDIR/hash.k2d" "$DBDIR/opts.k2d" "$DBDIR/taxo.k2d"

# Add FASTA files to Kraken2 library
kraken2-build --add-to-library "$PLASTID" --db "$DBDIR"
kraken2-build --add-to-library "$MITO" --db "$DBDIR"

# Download taxonomy (required before build)
kraken2-build --download-taxonomy --db "$DBDIR"

# Build the Kraken2 database
kraken2-build --build --threads $OMP_NUM_THREADS --db "$DBDIR"
