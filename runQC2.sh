#!/bin/bash -l
#SBATCH --job-name=trim_qc
#SBATCH --ntasks=32
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --output=/scratch/kmxg762/mixcr/logs/trim_qc_%j.out
#SBATCH --error=/scratch/kmxg762/mixcr/logs/trim_qc_%j.err

# Purpose: Trim FASTQ files with fastp, then run FastQC and MultiQC
# Input: Raw FASTQ files in specified directory
# Output: Trimmed FASTQ files, FastQC reports, and MultiQC summary

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
OUTPUT_DIR="/projects/be/us_techdev/TCR_discovery_kmxg762/Jurkat_mart1_alldata"
TRIMMED_DIR="${OUTPUT_DIR}/Trimmed_Fastq"
FASTQC_DIR="${OUTPUT_DIR}/FastQC_Trimmed"
MULTIQC_DIR="${OUTPUT_DIR}/MultiQC_Trimmed"
LOG_DIR="/scratch/kmxg762/mixcr/logs"

# Create output and log directories
mkdir -p "$TRIMMED_DIR" "$FASTQC_DIR" "$MULTIQC_DIR" "$LOG_DIR" || { echo "Error: Failed to create directories"; exit 1; }

# Check if directories are writable
for dir in "$TRIMMED_DIR" "$FASTQC_DIR" "$MULTIQC_DIR" "$LOG_DIR"; do
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

# Ensure fastp is available and executable
if [ ! -f "./fastp" ]; then
    echo "Downloading fastp..."
    wget -O fastp http://opengene.org/fastp/fastp || { echo "Error: Failed to download fastp"; exit 1; }
    chmod +x ./fastp
elif [ ! -x "./fastp" ]; then
    echo "Error: fastp is not executable. Please check permissions."
    exit 1
fi

echo "Starting FASTQ processing..."

# Process paired-end and single-end FASTQ files
for R1 in "$INPUT_DIR"/*_R1_*.fastq.gz; do
    [ -f "$R1" ] || continue
    R2="${R1/_R1_/_R2_}"
    BASE=$(basename "$R1" _R1_001.fastq.gz)

    if [ -f "$R2" ]; then
        # Paired-end processing
        OUT_R1="${TRIMMED_DIR}/${BASE}_R1_trimmed.fastq.gz"
        OUT_R2="${TRIMMED_DIR}/${BASE}_R2_trimmed.fastq.gz"
        REPORT="${TRIMMED_DIR}/${BASE}_fastp.html"
        JSON="${TRIMMED_DIR}/${BASE}_fastp.json"

        echo "Processing paired-end: $BASE"
        ./fastp \
            --in1 "$R1" \
            --in2 "$R2" \
            --out1 "$OUT_R1" \
            --out2 "$OUT_R2" \
            --thread "$SLURM_NTASKS" \
            --cut_tail \
            --cut_tail_mean_quality 28 \
            --qualified_quality_phred 28 \
            --length_required 50 \
            --adapter_sequence auto \
            --html "$REPORT" \
            --json "$JSON" > "$LOG_DIR/fastp_${BASE}.stdout.log" 2> "$LOG_DIR/fastp_${BASE}.stderr.log" || {
            echo "Error: fastp failed for paired-end $BASE. See $LOG_DIR/fastp_${BASE}.stderr.log"
            exit 1
        }
    else
        # Single-end processing
        OUT="${TRIMMED_DIR}/${BASE}_trimmed.fastq.gz"
        REPORT="${TRIMMED_DIR}/${BASE}_fastp.html"
        JSON="${TRIMMED_DIR}/${BASE}_fastp.json"

        echo "Processing single-end: $BASE"
        ./fastp \
            --in1 "$R1" \
            --out1 "$OUT" \
            --thread "$SLURM_NTASKS" \
            --cut_tail \
            --cut_tail_mean_quality 28 \
            --qualified_quality_phred 28 \
            --length_required 50 \
            --adapter_sequence auto \
            --html "$REPORT" \
            --json "$JSON" > "$LOG_DIR/fastp_${BASE}.stdout.log" 2> "$LOG_DIR/fastp_${BASE}.stderr.log" || {
            echo "Error: fastp failed for single-end $BASE. See $LOG_DIR/fastp_${BASE}.stderr.log"
            exit 1
        }
    fi
    echo "Trimming completed for $BASE"
done

# Process remaining single-end files
for FASTQ in "$INPUT_DIR"/*.fastq.gz; do
    [ -f "$FASTQ" ] || continue
    if [[ "$FASTQ" == *_R1_* || "$FASTQ" == *_R2_* ]]; then
        continue
    fi
    BASE=$(basename "$FASTQ" .fastq.gz)
    OUT="${TRIMMED_DIR}/${BASE}_trimmed.fastq.gz"
    REPORT="${TRIMMED_DIR}/${BASE}_fastp.html"
    JSON="${TRIMMED_DIR}/${BASE}_fastp.json"

    echo "Processing single-end: $BASE"
    ./fastp \
        --in1 "$FASTQ" \
        --out1 "$OUT" \
        --thread "$SLURM_NTASKS" \
        --cut_tail \
        --cut_tail_mean_quality 28 \
        --qualified_quality_phred 28 \
        --length_required 50 \
        --adapter_sequence auto \
        --html "$REPORT" \
        --json "$JSON" > "$LOG_DIR/fastp_${BASE}.stdout.log" 2> "$LOG_DIR/fastp_${BASE}.stderr.log" || {
        echo "Error: fastp failed for single-end $BASE. See $LOG_DIR/fastp_${BASE}.stderr.log"
        exit 1
    }
    echo "Trimming completed for $BASE"
done

# Run FastQC on trimmed files
echo "Running FastQC on trimmed files..."
fastqc -o "$FASTQC_DIR" -t "$SLURM_NTASKS" "$TRIMMED_DIR"/*.fastq.gz > "$LOG_DIR/fastqc_trimmed.stdout.log" 2> "$LOG_DIR/fastqc_trimmed.stderr.log" || {
    echo "Error: FastQC failed on trimmed files. See $LOG_DIR/fastqc_trimmed.stderr.log"
    exit 1
}

# Run MultiQC on FastQC and fastp results
echo "Running MultiQC on trimmed files..."
multiqc "$FASTQC_DIR" "$TRIMMED_DIR"/*.json -o "$MULTIQC_DIR" -n multiqc_report --pdf > "$LOG_DIR/multiqc_trimmed.stdout.log" 2> "$LOG_DIR/multiqc_trimmed.stderr.log" || {
    echo "Error: MultiQC or PDF export failed. See $LOG_DIR/multiqc_trimmed.stderr.log"
    exit 1
}

echo "Trimming, FastQC, and MultiQC completed successfully"
