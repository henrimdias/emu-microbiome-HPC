#!/bin/bash
#SBATCH --job-name=kraken2_build_final
#SBATCH --output=/home/jacks.local/henrique.dias/kraken_organelle_db/kraken2_build_%j.out
#SBATCH --error=/home/jacks.local/henrique.dias/kraken_organelle_db/kraken2_build_%j.err
#SBATCH --time=2-00:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --partition=compute
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=henrique.dias@sdstate.edu

# Load Conda
source ~/miniconda3/etc/profile.d/conda.sh
conda activate preemu_env

# Set environment variables
export OMP_NUM_THREADS=8

# Define paths
DBDIR="/home/jacks.local/henrique.dias/scratch/emu_pipeline/kraken2_db/db"
FASTA_DIR="/home/jacks.local/henrique.dias/scratch/kraken_env/data/organelle_refseq"
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
