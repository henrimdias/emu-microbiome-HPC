#!/bin/bash
#SBATCH --job-name=nanostat_reports
#SBATCH --output=/home/jacks.local/henrique.dias/16S_ITS_2025/NanoStat_reports/nanostat_reports_%j.log
#SBATCH --error=/home/jacks.local/henrique.dias/16S_ITS_2025/NanoStat_reports/nanostat_reports_%j.err
#SBATCH --time=03:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --partition=compute
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=henrique.dias@sdstate.edu

# Load conda environment
source activate nanotools_env

# Set directories
INPUT_DIR="/home/jacks.local/henrique.dias/scratch/16S_ITS_metabarcoding_2025/gsf0260_data/gsf0260_chris_16s_demuxed sequences/chris_16s_plate1_WHWLdemux"
OUTPUT_DIR="/home/jacks.local/henrique.dias/16S_ITS_2025/NanoStat_reports"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Look for .fastq rather than .fastq.gz
shopt -s nullglob
fastq_files=("$INPUT_DIR"/*.fastq)
if [[ ${#fastq_files[@]} -eq 0 ]]; then
    echo "No .fastq files found in $INPUT_DIR. Exiting."
    exit 1
fi

for fastq_path in "${fastq_files[@]}"; do
    filename=$(basename "$fastq_path")
    # strip the .fastq extension
    barcode_name="${filename%.fastq}"

    echo "Running NanoStat for $barcode_name..."
    barcode_output="$OUTPUT_DIR/$barcode_name"
    mkdir -p "$barcode_output"

    NanoStat \
      --fastq "$fastq_path" \
      --outdir "$barcode_output" \
      --name "$barcode_name"
done

echo "All NanoStat reports generated successfully!"