#!/bin/bash
#SBATCH --job-name=nanofilt_filtering
#SBATCH --output=/home/jacks.local/henrique.dias/scratch/emu_pipeline/nanofilt_reports/nanofilt_filtering_%j.log
#SBATCH --error=/home/jacks.local/henrique.dias/scratch/emu_pipeline/nanofilt_reports/nanofilt_filtering_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --partition=compute
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=henrique.dias@sdstate.edu

### 1) Activate your conda env properly ###
source /home/jacks.local/henrique.dias/miniconda3/etc/profile.d/conda.sh
conda activate nanotools_env

### 2) Paths (note quotes around any path with spaces) ###
MERGED_FASTQS="/home/jacks.local/henrique.dias/scratch/emu_pipeline/raw_data"
FILTERED_FASTQS="/home/jacks.local/henrique.dias/scratch/emu_pipeline/filtered_data"

mkdir -p "$FILTERED_FASTQS"

### 3) Gather both .fastq.gz *and* .fastq ###
shopt -s nullglob
inputs=( "$MERGED_FASTQS"/*.fastq.gz "$MERGED_FASTQS"/*.fastq )

if (( ${#inputs[@]} == 0 )); then
    echo "No FASTQ or FASTQ.GZ files found in $MERGED_FASTQS; exiting."
    exit 1
fi

### 4) Loop and filter ###
for infile in "${inputs[@]}"; do
    # derive a clean sample name:
    if [[ $infile == *.fastq.gz ]]; then
        sample=$(basename "$infile" .fastq.gz)
        decompress="zcat"
    else
        sample=$(basename "$infile" .fastq)
        decompress="cat"
    fi

    echo "Filtering $sample …"

    $decompress "$infile" \
      | NanoFilt -q 10 -l 1000 --maxlength 1850 \
      | gzip > "$FILTERED_FASTQS/${sample}_filtered.fastq.gz"

    echo "→ wrote $FILTERED_FASTQS/${sample}_filtered.fastq.gz"
done

echo "NanoFilt filtering complete!"
