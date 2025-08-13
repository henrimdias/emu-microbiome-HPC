#!/bin/bash
#SBATCH --job-name=merge_fastqs
#SBATCH --output=merge_fastqs_%j.log
#SBATCH --error=merge_fastqs_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --partition=compute
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=henrique.dias@sdstate.edu

# Load conda environment 
source activate nanotools_env

# Set directories
BASE_DIR="/home/jacks.local/henrique.dias/scratch/16S_ITS_metabarcoding_2025/gsf0260_data/Chris_16S_Native_Barcoding/Chris_16S/20250529_1144_P2S-01560-A_PBC46328_c7a003dc/fastq_pass"
MERGED_FASTQS="/home/jacks.local/henrique.dias/scratch/16S_merged_metabarcoding_2025"

# Create output directory
mkdir -p "$MERGED_FASTQS"

# Merge FASTQs
for barcode_dir in "$BASE_DIR"/barcode*/; do
    barcode_name=$(basename "$barcode_dir")
    echo "Merging $barcode_name..."

    merged_fastq="$MERGED_FASTQS/${barcode_name}_merged.fastq.gz"
    cat "$barcode_dir"/*.fastq.gz > "$merged_fastq"
done

echo "Done merging FASTQ files!"
