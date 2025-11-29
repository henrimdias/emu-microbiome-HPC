#!/bin/bash
# =============================================================================
# Article:
# Dias, H.M., et al. Reproducible Emu-based workflow for high-fidelity soil and
# plant microbiome profiling on HPC clusters. Bio-protocol. 2025.
#
# Script:
# End-to-end NanoTools QC workflow (NanoPlot → NanoFilt → NanoPlot → QC summary)
# for full-length 16S Nanopore reads.
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
# This SLURM pipeline:
#   - Runs NanoPlot on raw Nanopore FASTQ files.
#   - Filters reads by length and quality using NanoFilt.
#   - Runs NanoPlot again on filtered reads.
#   - Consolidates QC metrics into CSV/PDF summaries with a Python script.
#
# Inputs (user configuration section below):
#   RAW_DIR          : directory with raw FASTQ files
#   FILTERED_DIR     : directory for filtered FASTQ files
#   QC_DIR           : directory for NanoPlot QC reports
#   SUMMARY_DIR      : directory for summary outputs
#   SCRIPT_QC_SUMMARY: path to the Python QC summary script
#
# Outputs:
#   - NanoPlot reports for raw and filtered reads.
#   - Filtered FASTQ files.
#   - QC summary tables/plots in SUMMARY_DIR.
#
# Usage:
#   1) Edit paths and filtering parameters in the USER CONFIGURATION section.
#   2) Submit with:
#        sbatch nanotools_pipeline.slurm
#
# For full reproducibility, versions of NanoPlot, NanoFilt, NanoStat,
# Python, and dependencies are documented in the manuscript / docs.
# =============================================================================
#SBATCH --job-name=nanotools_pipeline
#SBATCH --output=/put/your/path/emu_pipeline/nanotools_logs/nanotools_pipeline_%j.out
#SBATCH --error=/put/your/path//emu_pipeline/nanotools_logs/nanotools_pipeline_%j.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --partition=compute
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=put-your-email@email.edu

##############################################################
#                NANOTOOLS QC WORKFLOW PIPELINE
#    NanoPlot → NanoFilt → NanoPlot → Summary consolidation
#
#  VERSION: 1.0
#  AUTHOR: Henrique M. Dias (South Dakota State University)
#  UPDATED: November 2025
#
#  DESCRIPTION:
#    This SLURM pipeline performs complete preprocessing and QC
#    of full-length 16S rRNA reads generated on Nanopore devices.
#
#  STEPS:
#    1. Generate QC reports for RAW reads (NanoPlot)
#    2. Filter reads by quality and length (NanoFilt)
#    3. Generate QC reports for FILTERED reads (NanoPlot)
#    4. Consolidate QC metrics into a master summary (Python)
#
#  REQUIREMENTS:
#    - Conda environment: nanotools_env
#      containing NanoPlot, NanoFilt, NanoStat, and Python libs
#    - Input: gzipped FASTQ files (*.fastq.gz)
#    - Output: QC reports (HTML, PNG) + filtered FASTQ + summary CSV
#
#  USAGE:
#    1. Edit ONLY the section below ("USER CONFIGURATION")
#    2. Submit job:
#         sbatch run_nanotools_pipeline.slurm
#    3. Review logs in: logs/nanotools_pipeline_<JOBID>.out
#
#  NOTE:
#    - All folders will be created automatically if missing.
#    - Modify MIN_LEN, MAX_LEN, and MIN_QUAL if needed.
#    - QC summary ensures expected read-length window (1450–1550 bp)
#      and read retention (60–90%) after filtering.
##############################################################

# ==============================
# USER CONFIGURATION
# ==============================
# ⚠️ EDIT ONLY THIS SECTION ⚠️

RAW_DIR="/put/your/path/emu_pipeline/raw_data"                 # Folder with raw FASTQ files
FILTERED_DIR="/put/your/path/emu_pipeline/filtered_data"       # Output folder for filtered FASTQ
QC_DIR="/put/your/path/emu_pipeline/nanoplot_reports"          # QC output folder (NanoPlot)
SUMMARY_DIR="/put/your/path/emu_pipeline/qc_summary_reports"   # Folder for summary outputs
SCRIPT_QC_SUMMARY="/put/your/path/emu_pipeline/scripts/run_QC_summary.py"  # Path to QC summary script

# Filtering parameters
MIN_LEN=1000   # Minimum read length to retain
MAX_LEN=1850   # Maximum read length to retain
MIN_QUAL=10    # Minimum mean read quality (Q-score)

# ==============================
# ENVIRONMENT SETUP
# ==============================

# initialize conda for non-interactive shells
source /put/your/path/miniconda3/etc/profile.d/conda.sh  
# now you can activate
conda activate nanotools_env

echo "=============================================================="
echo " NanoTools QC Workflow started on $(date)"
echo "--------------------------------------------------------------"
echo " RAW reads directory:       $RAW_DIR"
echo " FILTERED reads directory:  $FILTERED_DIR"
echo " QC reports directory:      $QC_DIR"
echo " SUMMARY directory:         $SUMMARY_DIR"
echo "=============================================================="

# ==============================
# STEP A: NanoPlot on RAW READS
# ==============================
echo ""
echo ">>> STEP A: Generating QC reports for RAW reads (NanoPlot)"
echo "--------------------------------------------------------------"

# Collect all FASTQ files (supports .fastq, .fq, and gzipped versions)
RAW_FILES=($(find "$RAW_DIR" -maxdepth 1 -type f \( -name "*.fastq" -o -name "*.fq" -o -name "*.fastq.gz" -o -name "*.fq.gz" \)))

# Check if any files were found
if [ ${#RAW_FILES[@]} -eq 0 ]; then
    echo "❌ ERROR: No FASTQ files found in $RAW_DIR"
    echo "   Expected extensions: .fastq, .fastq.gz, .fq, or .fq.gz"
    exit 1
fi

for FILE in "${RAW_FILES[@]}"; do
    SAMPLE=$(basename "$FILE")
    SAMPLE="${SAMPLE%%.*}"   # strip extension safely
    OUTPUT_DIR="${QC_DIR}/raw/${SAMPLE}"
    mkdir -p "$OUTPUT_DIR"

    echo "Processing sample: $SAMPLE"
    NanoPlot --fastq "$FILE" --threads 8 \
        --outdir "$OUTPUT_DIR" \
        --plots hex dot \
        --maxlength 2000
done
echo ">>> NanoPlot (raw) completed successfully."
echo ""

# ==============================
# STEP B: NanoFilt Filtering
# ==============================
echo ">>> STEP B: Filtering FASTQ files with NanoFilt"
echo "--------------------------------------------------------------"
echo "Filtering parameters:"
echo "  - Min length: $MIN_LEN bp"
echo "  - Max length: $MAX_LEN bp"
echo "  - Min quality: Q$MIN_QUAL"
echo ""

for FILE in "${RAW_FILES[@]}"; do
    SAMPLE=$(basename "$FILE")
    SAMPLE="${SAMPLE%%.*}"   # Remove extension
    OUTFILE="${FILTERED_DIR}/${SAMPLE}_filtered.fastq.gz"

    echo "Filtering sample: $SAMPLE"

    # Detect compression automatically
    if [[ "$FILE" == *.gz ]]; then
        zcat "$FILE" | NanoFilt -q $MIN_QUAL --length $MIN_LEN --maxlength $MAX_LEN | gzip > "$OUTFILE"
    else
        cat "$FILE" | NanoFilt -q $MIN_QUAL --length $MIN_LEN --maxlength $MAX_LEN | gzip > "$OUTFILE"
    fi
done
echo ">>> NanoFilt filtering completed successfully."
echo ""

# ==============================
# STEP C: NanoPlot on FILTERED READS
# ==============================
echo ">>> STEP C: Generating QC reports for FILTERED reads (NanoPlot)"
echo "--------------------------------------------------------------"

FILTERED_FILES=($(find "$FILTERED_DIR" -maxdepth 1 -type f \( -name "*_filtered.fastq" -o -name "*_filtered.fq" -o -name "*_filtered.fastq.gz" -o -name "*_filtered.fq.gz" \)))

if [ ${#FILTERED_FILES[@]} -eq 0 ]; then
    echo "❌ ERROR: No filtered FASTQ files found in $FILTERED_DIR"
    exit 1
fi

for FILE in "${FILTERED_FILES[@]}"; do
    SAMPLE=$(basename "$FILE" | sed 's/_filtered.*//')
    OUTPUT_DIR="${QC_DIR}/filtered/${SAMPLE}"
    mkdir -p "$OUTPUT_DIR"

    echo "Processing filtered sample: $SAMPLE"
    NanoPlot --fastq "$FILE" --threads 8 \
        --outdir "$OUTPUT_DIR" \
        --plots hex dot \
        --maxlength 2000
done
echo ">>> NanoPlot (filtered) completed successfully."


# ==============================
# STEP D: Consolidate QC Reports
# ==============================
echo ">>> STEP D: Consolidating QC reports with Python summary script"
echo "--------------------------------------------------------------"
python "$SCRIPT_QC_SUMMARY" --input "$QC_DIR" --output "$SUMMARY_DIR"

if [ $? -eq 0 ]; then
    echo ">>> QC summary completed successfully."
else
    echo "⚠️  WARNING: QC summary encountered an error. Check the log output."
fi

echo ""
echo "=============================================================="
echo " NanoTools QC Workflow finished on $(date)"
echo " All results saved under: $SUMMARY_DIR"
echo " Review logs and QC outputs for verification."
echo "=============================================================="
