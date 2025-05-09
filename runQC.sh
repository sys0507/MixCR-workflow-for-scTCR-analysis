#!/bin/bash -l
#SBATCH --job-name=fastqc
#SBATCH --ntasks=32
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --output=/scratch/kmxg762/mixcr/logs/fastqc_%j.out
#SBATCH --error=/scratch/kmxg762/mixcr/logs/fastqc_%j.err

# Purpose: Run FastQC and MultiQC on raw FASTQ files for quality control
# Input: Raw FASTQ files in specified directory
# Output: FastQC reports and MultiQC summary

# Exit on any error
set -e

# Change to working directory
cd /scratch/kmxg762/mixcr || { echo "Error: Cannot change to /scratch/kmxg762/mixcr"; exit 1; }

# Load required modules
module load GCC/10.3.0
module load FastQC/0.12.1-Java-11.0.2
module load MultiQC/1.25.1-foss-2021a

# Define directories
INPUT_DIR="/projects/be/us_techdev/TCR_discovery_kmxg762/Jurkat_mart1_alldata/Fastq"
FASTQC_DIR="/projects/be/us_techdev/TCR_discovery_kmxg762/Jurkat_mart1_alldata/FastQC"
MULTIQC_DIR="/projects/be/us_techdev/TCR_discovery_kmxg762/Jurkat_mart1_alldata/MultiQC"
LOG_DIR="/scratch/kmxg762/mixcr/logs"

# Create output and log directories
mkdir -p "$FASTQC_DIR" "$MULTIQC_DIR" "$LOG_DIR" || { echo "Error: Failed to create directories"; exit 1; }

# Check if directories are writable
for dir in "$FASTQC_DIR" "$MULTIQC_DIR" "$LOG_DIR"; do
    if [ ! -w "$dir" ]; then
        echo "Error: Directory $dir is not writable"
        exit 1
    fi
done

# Check for input FASTQ files
if ! ls "$INPUT_DIR"/*.fastq.gz >/dev/null 2>&1; then
    echo "Error: No FASTQ files found in $INPUT_DIR"
    exit 1
fi

# Run FastQC
echo "Running FastQC on FASTQ files..."
fastqc -o "$FASTQC_DIR" -t "$SLURM_NTASKS" "$INPUT_DIR"/Jurkat-P0321-Mart1*.fastq.gz > "$LOG_DIR/fastqc.stdout.log" 2> "$LOG_DIR/fastqc.stderr.log" || {
    echo "Error: FastQC failed. See $LOG_DIR/fastqc.stderr.log for details"
    exit 1
}

# Run MultiQC
echo "Running MultiQC on FastQC results..."
multiqc -o "$MULTIQC_DIR" "$FASTQC_DIR" -n multiqc_report --pdf > "$LOG_DIR/multiqc.stdout.log" 2> "$LOG_DIR/multiqc.stderr.log" || {
    echo "Error: MultiQC or PDF export failed. See $LOG_DIR/multiqc.stderr.log for details"
    exit 1
}

echo "FastQC and MultiQC completed successfully"
