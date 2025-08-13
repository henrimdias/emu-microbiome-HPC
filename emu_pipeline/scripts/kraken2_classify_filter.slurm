#!/bin/bash
#SBATCH --job-name=kraken2_filter_reads
#SBATCH --output=/home/jacks.local/henrique.dias/scratch/emu_pipeline/no_organelle_filtered_data/kraken2_filter_%j.out
#SBATCH --error=/home/jacks.local/henrique.dias/scratch/emu_pipeline/no_organelle_filtered_data/kraken2_filter_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --partition=compute
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=henrique.dias@sdstate.edu

# Load environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate preemu_env

# Define input/output
READ_DIR="/home/jacks.local/henrique.dias/scratch/emu_pipeline/filtered_data"
DBDIR="/home/jacks.local/henrique.dias/scratch/emu_pipeline/kraken2_db/db"
OUTDIR="/home/jacks.local/henrique.dias/scratch/emu_pipeline/no_organelle_filtered_data"
SCRIPT="/home/jacks.local/henrique.dias/scratch/emu_pipeline/scripts/remove_kraken2_organelle_reads.py"

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
