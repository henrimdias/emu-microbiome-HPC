#!/bin/bash
#SBATCH --job-name=nanoplot_reports
#SBATCH --output=/home/jacks.local/henrique.dias/scratch/emu_pipeline/nanoplot_reports/nanoplot_log/nanoplot_reports_%j.log
#SBATCH --error=/home/jacks.local/henrique.dias/scratch/emu_pipeline/nanoplot_reports/nanoplot_log/nanoplot_reports_%j.err
#SBATCH --time=03:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --partition=compute
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=henrique.dias@sdstate.edu

# initialize conda for non-interactive shells
source /home/jacks.local/henrique.dias/miniconda3/etc/profile.d/conda.sh  
# now you can activate
conda activate nanotools_env


# Set directories
INPUT_DIR="/home/jacks.local/henrique.dias/scratch/emu_pipeline/filtered_data"
OUTPUT_DIR="/home/jacks.local/henrique.dias/scratch/emu_pipeline/nanoplot_after_nanofilt_data"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Enable nullglob so that unmatched globs disappear instead of staying literal
shopt -s nullglob

# Collect both .fastq and .fastq.gz files
fastq_files=("$INPUT_DIR"/*.{fastq,fastq.gz})

if [[ ${#fastq_files[@]} -eq 0 ]]; then
    echo "No .fastq or .fastq.gz files found in $INPUT_DIR. Exiting."
    exit 1
fi

for fastq_path in "${fastq_files[@]}"; do
    filename=$(basename "$fastq_path")
    # Strip the correct extension
    if [[ "$filename" == *.fastq.gz ]]; then
        sample="${filename%.fastq.gz}"
    else
        sample="${filename%.fastq}"
    fi

    echo "Running NanoPlot for $sample…"
    sample_out="$OUTPUT_DIR/$sample"
    mkdir -p "$sample_out"

    NanoPlot \
      --outdir "$sample_out" \
      --fastq "$fastq_path" \
      --threads ${SLURM_CPUS_PER_TASK} \
      --prefix "$sample" \
      --loglength \
      --plots kde hex dot

    echo "→ Finished NanoPlot report for $sample"
done

echo "All NanoPlot reports generated successfully!"



