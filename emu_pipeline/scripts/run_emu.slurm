#!/bin/bash
#SBATCH --job-name=emu_array
#SBATCH --partition=compute
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=2-00:00:00
#SBATCH --array=0-62%62  # Adjust range based on number of input FASTQ files
#SBATCH --output=/home/jacks.local/henrique.dias/scratch/emu_pipeline/emu_results/emu_nofilt%A_%a.out
#SBATCH --error=/home/jacks.local/henrique.dias/scratch/emu_pipeline/emu_results/emu_nofilt%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=henrique.dias@sdstate.edu

# -----------------------------------
# 1. Initialize Conda environment
# -----------------------------------
source ~/miniconda3/etc/profile.d/conda.sh
conda activate emu_env

# -----------------------------------
# 2. Define paths
# -----------------------------------
BASE_DIR="/home/jacks.local/henrique.dias/scratch/emu_pipeline/no_organelle_filtered_data"
OUTPUT_DIR="/home/jacks.local/henrique.dias/scratch/emu_pipeline/emu_results"
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
