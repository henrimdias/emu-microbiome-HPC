#!/bin/bash
#SBATCH --job-name=kraken2_db
#SBATCH --output=/home/jacks.local/henrique.dias/scratch/emu_pipeline/kraken2_db/download_kraken2_db_%j.outt
#SBATCH --error/home/jacks.local/henrique.dias/scratch/emu_pipeline/kraken2_db/download_kraken2_db_%j.err
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
DBDIR=/home/jacks.local/henrique.dias/scratch/emu_pipeline/kraken2_db

# Step 1: Download taxonomy
kraken2-build --download-taxonomy --db $DBDIR

# Step 2: Download organellar libraries
kraken2-build --download-library chloroplast --db $DBDIR
kraken2-build --download-library mitochondrion --db $DBDIR

# Step 3: Build the database
kraken2-build --build --threads 8 --db $DBDIR
