!/bin/bash
# =============================================================================
# Article:
# Dias, H.M., et al. Reproducible Emu-based workflow for high-fidelity soil and
# plant microbiome profiling on HPC clusters. Bio-protocol. 2025.
#
# Script:
# Build organelle-focused Kraken2 database for plastid/mitochondrial filtering
# within the workflow described in the article above.
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
# This SLURM batch script downloads NCBI taxonomy and organellar sequences
# (chloroplast and mitochondrion) and builds a Kraken2 database for removing
# plant organelle reads prior to EMU-based microbiome profiling.
#
# The script performs:
#   - Download of NCBI taxonomy for Kraken2
#   - Download of chloroplast and mitochondrion libraries
#   - Build of a Kraken2 database using the specified number of threads
#
# Assumptions:
#   - Kraken2 is installed and available in the conda environment `kraken_env`
#   - Miniconda is installed and accessible at ~/miniconda3
#   - The target database directory exists or can be created by Kraken2
#
# Inputs (user configuration below):
#   DBDIR : path to the Kraken2 database directory
#
# Outputs:
#   - Populated Kraken2 database under $DBDIR containing taxonomy and
#     chloroplast/mitochondrion libraries
#
# Usage:
#   sbatch download_build_kraken2_db.sh
#
# For full reproducibility, the versions of Kraken2 and dependencies are listed
# in the manuscript / accompanying documentation.
# =============================================================================
#SBATCH --job-name=kraken2_db
#SBATCH --output=/put/your/path/emu_pipeline/kraken2_db/download_kraken2_db_%j.outt
#SBATCH --error=/put/your/path/emu_pipeline/kraken2_db/download_kraken2_db_%j.err
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=compute
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=henrique.dias@sdstate.edu

# Load Conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate kraken_env

# Set DB path
DBDIR=/put/your/path/h/emu_pipeline/kraken2_db

# Step 1: Download taxonomy
kraken2-build --download-taxonomy --db $DBDIR

# Step 2: Download organellar libraries
kraken2-build --download-library chloroplast --db $DBDIR
kraken2-build --download-library mitochondrion --db $DBDIR

# Step 3: Build the database
kraken2-build --build --threads 8 --db $DBDIR
