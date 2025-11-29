#!/bin/bash
# =============================================================================
# Article:
# Dias, H.M., et al. Reproducible Emu-based workflow for high-fidelity soil and
# plant microbiome profiling on HPC clusters. Bio-protocol. 2025.
#
# Script:
# Run EMU full-length 16S abundance profiling as a SLURM job array on
# organelle-filtered FASTQ files.
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
# This SLURM array script iterates over organelle-depleted FASTQ files and
# runs EMU abundance profiling for each sample in parallel. For every input
# *_filtered.fastq file, EMU is executed with options to keep read-level
# counts and assignments, and outputs are stored in per-sample folders.
#
# The script:
#   - Activates the `emu_env` Conda environment.
#   - Uses a SLURM array index to select one FASTQ per task.
#   - Sets EMU_DATABASE_DIR to point to the pre-built EMU database.
#   - Runs `emu abundance` with `--keep-counts` and `--keep-read-assignments`.
#
# Assumptions:
#   - Miniconda is installed and accessible at ~/miniconda3.
#   - The Conda environment `emu_env` is created and contains EMU and its
#     dependencies.
#   - BASE_DIR contains organelle-filtered FASTQ files matching *_filtered.fastq.
#   - EMU_DB points to a valid EMU full-length 16S database.
#   - The SLURM array bounds match the number of FASTQ files in BASE_DIR.
#
# Inputs (user configuration below):
#   BASE_DIR   : directory containing *_filtered.fastq files
#   OUTPUT_DIR : directory where EMU result subfolders will be written
#   EMU_DB     : EMU database directory (exported as EMU_DATABASE_DIR)
#
# Outputs (per sample):
#   - EMU abundance tables and read assignments under:
#       $OUTPUT_DIR/${BARCODE_NAME}/
#
# Usage:
#   Adjust the --array range to match the number of input FASTQ files:
#     FASTA_FILES=($(...))
#     # number of files = N
#     #SBATCH --array=0-(N-1)%N
#   Then submit:
#     sbatch run_emu.slurm
#
# For full reproducibility, EMU, database, and dependency versions are
# documented in the manuscript / accompanying documentation.
# =============================================================================

#SBATCH --job-name=emu_array
#SBATCH --partition=compute
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=2-00:00:00
#SBATCH --array=0-62%62  # Adjust range based on number of input FASTQ files
#SBATCH --output=/put/your/path/emu_pipeline/emu_results/emu_nofilt%A_%a.out
#SBATCH --error=/put/your/path/emu_pipeline/emu_results/emu_nofilt%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=put-your-emails@email.edu

# -----------------------------------
# 1. Initialize Conda environment
# -----------------------------------
source ~/miniconda3/etc/profile.d/conda.sh
conda activate emu_env

# -----------------------------------
# 2. Define paths
# -----------------------------------
BASE_DIR="/put/your/path/emu_pipeline/no_organelle_filtered_data"
OUTPUT_DIR="/put/your/path/emu_pipeline/emu_results"
EMU_DB="/home/jacks.local/henrique.dias/scratch/emu_pipeline/emu_db"

export EMU_DATABASE_DIR=$EMU_DB

# -----------------------------------
# 3. Get list of all relevant FASTQ files
# -----------------------------------
FASTA_FILES=($(find "$BASE_DIR" -maxdepth 1 -type f -name "*_filtered.fastq"))

# -----------------------------------
# 4. Get the current FASTQ file based on SLURM array index
# -----------------------------------
FASTA_FILE=${FASTA_FILES[$SLURM_ARRAY_TASK_ID]}

# -----------------------------------
# 5. Define output subfolder based on input file name
# -----------------------------------
BARCODE_NAME=$(basename "$FASTA_FILE" .fastq)
BARCODE_OUTPUT="$OUTPUT_DIR/$BARCODE_NAME"
mkdir -p "$BARCODE_OUTPUT"

# -----------------------------------
# 6. Run EMU abundance profiling
# -----------------------------------
echo "Running EMU on: $FASTA_FILE"
emu abundance \
  --output-dir "$BARCODE_OUTPUT" \
  --keep-counts \
  --keep-read-assignments \
  "$FASTA_FILE"

# -----------------------------------
# 7. Clean up
# -----------------------------------
conda deactivate
