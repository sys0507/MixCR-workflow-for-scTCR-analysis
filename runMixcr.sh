#!/bin/bash -l
#SBATCH --job-name=mixcr_analysis
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --partition=medium
#SBATCH --output=/scratch/kmxg762/mixcr/logs/mixcr_%j.out
#SBATCH --error=/scratch/kmxg762/mixcr/logs/mixcr_%j.err

# Purpose: Run MiXCR analysis on trimmed FASTQ files, rename outputs, and generate QC reports
# Input: Trimmed FASTQ files in specified directory
# Output: MiXCR analysis results, renamed files, and QC reports (align, chainUsage, coverage, tags)

# Exit on any error
set -e

# Change to working directory
cd /scratch/kmxg762/mixcr || { echo "Error: Cannot change to /scratch/kmxg762/mixcr"; exit 1; }

# Load MiXCR module
module load MiXCR/4.7.0-Java-11.0.2 || { echo "Error: Failed to load MiXCR module"; exit 1; }
export PATH=/scratch/kmxg762/mixcr:$PATH

# Define directories
FASTQ_DIR="/projects/be/us_techdev/TCR_discovery_kmxg762/Jurkat_mart1_alldata/Trimmed_Fastq"
OUTPUT_DIR="/projects/be/us_techdev/TCR_discovery_kmxg762/Jurkat_mart1_alldata/data/mixcr_output_scTCR"
LOG_DIR="${OUTPUT_DIR}/logs"
QC_DIR="${OUTPUT_DIR}/qc_reports"

# Create output, log, and QC directories
mkdir -p "$OUTPUT_DIR" "$LOG_DIR" "$QC_DIR" || { echo "Error: Failed to create directories"; exit 1; }

# Check if directories are writable
for dir in "$OUTPUT_DIR" "$LOG_DIR" "$QC_DIR"; do
    if [ ! -w "$dir" ]; then
        echo "Error: Directory $dir is not writable"
        exit 1
    fi
done

# Step 1: Verify FASTQ files
echo "Verifying FASTQ files in $FASTQ_DIR..."
FASTQ_FILES=($(find "$FASTQ_DIR" -type f -name '*_L001_R[1-2]_trimmed.fastq.gz'))
if [ ${#FASTQ_FILES[@]} -eq 0 ]; then
    echo "Error: No FASTQ files found in $FASTQ_DIR"
    exit 1
fi

# Log FASTQ files and verify R1/R2 pairs
SAMPLE_IDS=()
declare -A SAMPLE_R1_FILES
declare -A SAMPLE_R2_FILES
{
    echo "Found ${#FASTQ_FILES[@]} FASTQ files:"
    for file in "${FASTQ_FILES[@]}"; do
        echo "  $file"
        if [[ $(basename "$file") =~ (.+_S[0-9]{1,4})_L001_R([1-2])_trimmed\.fastq\.gz$ ]]; then
            sample_id="${BASH_REMATCH[1]}"
            read_type="R${BASH_REMATCH[2]}"
            SAMPLE_IDS+=("$sample_id")
            if [ "$read_type" == "R1" ]; then
                SAMPLE_R1_FILES["$sample_id"]="$file"
            else
                SAMPLE_R2_FILES["$sample_id"]="$file"
            fi
        fi
    done
} > "$LOG_DIR/fastq_files.log"

# Remove duplicates and log unique sample IDs
UNIQUE_SAMPLE_IDS=($(printf '%s\n' "${SAMPLE_IDS[@]}" | sort -u))
echo "Found ${#UNIQUE_SAMPLE_IDS[@]} unique sample IDs:" >> "$LOG_DIR/fastq_files.log"
printf '  %s\n' "${UNIQUE_SAMPLE_IDS[@]}" >> "$LOG_DIR/fastq_files.log"

# Verify R1/R2 pairs
for sample_id in "${UNIQUE_SAMPLE_IDS[@]}"; do
    if [ -z "${SAMPLE_R1_FILES[$sample_id]}" ] || [ -z "${SAMPLE_R2_FILES[$sample_id]}" ]; then
        echo "Error: Missing R1 or R2 file for sample $sample_id"
        exit 1
    fi
done

# Step 2: Run MiXCR analysis
echo "Running MiXCR analysis..."
mixcr -Xmx64g analyze generic-lt-single-cell-amplicon \
    --species hs \
    --rna \
    --rigid-left-alignment-boundary FR3Begin \
    --floating-right-alignment-boundary C \
    --assemble-clonotypes-by CDR3 \
    --split-by-sample \
    "${FASTQ_DIR}/Jurkat-P0321-Mart1-{{n}}-{{CELL:a}}_{{SAMPLE:a}}_L001_{{R}}_trimmed.fastq.gz" \
    "$OUTPUT_DIR/results" \
    > "$LOG_DIR/mixcr_run.stdout.log" 2> "$LOG_DIR/mixcr_run.stderr.log" || {
    echo "Error: MiXCR analysis failed. See $LOG_DIR/mixcr_run.stderr.log"
    exit 1
}

# Step 2.1: Verify .vdjca files
VDJCA_FILES=($(find "$OUTPUT_DIR" -type f -name 'results.*.vdjca'))
if [ ${#VDJCA_FILES[@]} -eq 0 ]; then
    echo "Error: No .vdjca files generated in $OUTPUT_DIR"
    exit 1
fi

# Step 2.2: Extract sample metadata
echo "Extracting sample metadata..."
for file in "${VDJCA_FILES[@]}"; do
    sample=$(basename "$file" .vdjca)
    mixcr -Xmx64g exportReports \
        --json \
        "$file" \
        "${LOG_DIR}/${sample}_metadata.json" > "${LOG_DIR}/${sample}_metadata.stdout.log" 2> "${LOG_DIR}/${sample}_metadata.stderr.log" || {
        echo "Warning: Failed to extract metadata for $sample"
    }
done

# Step 3: Rename output files
echo "Renaming MiXCR output files..."
declare -A sample_mapping
while IFS= read -r fastq_file; do
    filename=$(basename "$fastq_file")
    if [[ $filename =~ (.+_S[0-9]{1,4})_L001_R1_trimmed\.fastq\.gz$ ]]; then
        full_sample_name="${BASH_REMATCH[1]}"
        short_sample_name=$(echo "$full_sample_name" | grep -o 'S[0-9]\{1,4\}')
        sample_mapping["$short_sample_name"]="$full_sample_name"
    fi
done < <(find "$FASTQ_DIR" -type f -name '*_S[0-9]*_L001_R1_trimmed.fastq.gz')

for ext in refined.vdjca alignments.vdjca assembledCells.clns clns clna txt tsv json; do
    while IFS= read -r output_file; do
        filename=$(basename "$output_file")
        if [[ $filename =~ results\.(S[0-9]{1,4})\.(.*)$ ]]; then
            short_sample_name="${BASH_REMATCH[1]}"
            suffix="${BASH_REMATCH[2]}"
            if [ -n "${sample_mapping[$short_sample_name]}" ]; then
                new_filename="results.${sample_mapping[$short_sample_name]}.${suffix}"
                mv "$output_file" "${OUTPUT_DIR}/${new_filename}" || {
                    echo "Error: Failed to rename $output_file"
                    exit 1
                }
            else
                echo "Warning: No full sample name for $short_sample_name in $filename"
            fi
        fi
    done < <(find "$OUTPUT_DIR" -type f -name "results.S[0-9]*.$ext")
done

# Step 4: Generate QC reports
REFINED_VDJCA_FILES=($(find "$OUTPUT_DIR" -type f -name '*.refined.vdjca'))
ALIGNMENTS_VDJCA_FILES=($(find "$OUTPUT_DIR" -type f -name '*.alignments.vdjca')))
ASSEMBLED_CELLS_CLNS_FILES=($(find "$OUTPUT_DIR" -type f -name '*.assembledCells.clns')))
CLNS_FILES=($(find "$OUTPUT_DIR" -type f -name '*.clns')))
CLNA_FILES=($(find "$OUTPUT_DIR" -type f -name '*.clna')))

# Generate align reports
if [ ${#REFINED_VDJCA_FILES[@]} -gt 0 ]; then
    mixcr -Xmx64g exportQc align \
        --absolute-values --width 1200 --height 800 --force-overwrite --verbose \
        "${REFINED_VDJCA_FILES[@]}" "${QC_DIR}/all_samples_align_refined.pdf" > "${LOG_DIR}/align_refined.stdout.log" 2> "${LOG_DIR}/align_refined.stderr.log" || {
        echo "Error: Failed to generate align QC for refined.vdjca"
        exit 1
    }
    for file in "${REFINED_VDJCA_FILES[@]}"; do
        sample=$(basename "$file" .refined.vdjca)
        mixcr -Xmx64g exportQc align \
            --absolute-values --width 1200 --height 800 --force-overwrite --verbose \
            "$file" "${QC_DIR}/${sample}_align_refined.pdf" > "${LOG_DIR}/${sample}_align_refined.stdout.log" 2> "${LOG_DIR}/${sample}_align_refined.stderr.log" || {
            echo "Error: Failed to generate align QC for $sample"
        }
    done
fi

if [ ${#ALIGNMENTS_VDJCA_FILES[@]} -gt 0 ]; then
    mixcr -Xmx64g exportQc align \
        --absolute-values --width 1200 --height 800 --force-overwrite --verbose \
        "${ALIGNMENTS_VDJCA_FILES[@]}" "${QC_DIR}/all_samples_align_alignments.pdf" > "${LOG_DIR}/align_alignments.stdout.log" 2> "${LOG_DIR}/align_alignments.stderr.log" || {
        echo "Error: Failed to generate align QC for alignments.vdjca"
        exit 1
    }
    for file in "${ALIGNMENTS_VDJCA_FILES[@]}"; do
        sample=$(basename "$file" .alignments.vdjca)
        mixcr -Xmx64g exportQc align \
            --absolute-values --width 1200 --height 800 --force-overwrite --verbose \
            "$file" "${QC_DIR}/${sample}_align_alignments.pdf" > "${LOG_DIR}/${sample}_align_alignments.stdout.log" 2> "${LOG_DIR}/${sample}_align_alignments.stderr.log" || {
            echo "Error: Failed to generate align QC for $sample"
        }
    done
fi

# Generate chainUsage reports
if [ ${#ASSEMBLED_CELLS_CLNS_FILES[@]} -gt 0 ]; then
    mixcr -Xmx64g exportQc chainUsage \
        --absolute-values --align-chain-usage --hide-non-functional --width 1200 --height 800 --force-overwrite --verbose \
        "${ASSEMBLED_CELLS_CLNS_FILES[@]}" "${QC_DIR}/all_samples_chainUsage_assembledCells.pdf" > "${LOG_DIR}/chainUsage_assembledCells.stdout.log" 2> "${LOG_DIR}/chainUsage_assembledCells.stderr.log" || {
        echo "Error: Failed to generate chainUsage QC for assembledCells.clns"
        exit 1
    }
    for file in "${ASSEMBLED_CELLS_CLNS_FILES[@]}"; do
        sample=$(basename "$file" .assembledCells.clns)
        mixcr -Xmx64g exportQc chainUsage \
            --absolute-values --align-chain-usage --hide-non-functional --width 1200 --height 800 --force-overwrite --verbose \
            "$file" "${QC_DIR}/${sample}_chainUsage_assembledCells.pdf" > "${LOG_DIR}/${sample}_chainUsage_assembledCells.stdout.log" 2> "${LOG_DIR}/${sample}_chainUsage_assembledCells.stderr.log" || {
            echo "Error: Failed to generate chainUsage QC for $sample"
        }
    done
fi

if [ ${#CLNS_FILES[@]} -gt 0 ] || [ ${#CLNA_FILES[@]} -gt 0 ]; then
    CLN_FILES=("${CLNS_FILES[@]}" "${CLNA_FILES[@]}")
    mixcr -Xmx64g exportQc chainUsage \
        --absolute-values --align-chain-usage --hide-non-functional --width 1200 --height 800 --force-overwrite --verbose \
        "${CLN_FILES[@]}" "${QC_DIR}/all_samples_chainUsage_clns.pdf" > "${LOG_DIR}/chainUsage_clns.stdout.log" 2> "${LOG_DIR}/chainUsage_clns.stderr.log" || {
        echo "Error: Failed to generate chainUsage QC for clns/clna"
        exit 1
    }
    for file in "${CLN_FILES[@]}"; do
        sample=$(basename "$file" | sed 's/\.\(clns\|clna\)$//')
        mixcr -Xmx64g exportQc chainUsage \
            --absolute-values --align-chain-usage --hide-non-functional --width 1200 --height 800 --force-overwrite --verbose \
            "$file" "${QC_DIR}/${sample}_chainUsage_clns.pdf" > "${LOG_DIR}/${sample}_chainUsage_clns.stdout.log" 2> "${LOG_DIR}/${sample}_chainUsage_clns.stderr.log" || {
            echo "Error: Failed to generate chainUsage QC for $sample"
        }
    done
fi

# Generate coverage and tags reports
FAILED_SAMPLES=()
for file in "${REFINED_VDJCA_FILES[@]}" "${ALIGNMENTS_VDJCA_FILES[@]}" "${CLNA_FILES[@]}"; do
    [ -f "$file" ] || continue
    sample=$(basename "$file" | sed 's/\.\(refined\.vdjca\|alignments\.vdjca\|clna\)$//')

    if [[ "$file" == *.refined.vdjca || "$file" == *.alignments.vdjca ]]; then
        mixcr -Xmx64g exportQc coverage \
            --show-boundaries --force-overwrite --verbose \
            "$file" "${QC_DIR}/${sample}_coverage.pdf" > "${LOG_DIR}/${sample}_coverage.stdout.log" 2> "${LOG_DIR}/${sample}_coverage.stderr.log" || {
            echo "Error: Failed to generate coverage QC for $sample"
            FAILED_SAMPLES+=("${sample}_coverage")
        }
    fi

    mixcr -Xmx64g exportQc tags \
        --force-overwrite --verbose \
        "$file" "${QC_DIR}/${sample}_tags.pdf" > "${LOG_DIR}/${sample}_tags.stdout.log" 2> "${LOG_DIR}/${sample}_tags.stderr.log" || {
        echo "Error: Failed to generate tags QC for $sample"
        FAILED_SAMPLES+=("${sample}_tags")
    }
done

# Step 5: Log read counts
echo "Logging read counts per sample..."
for file in "${REFINED_VDJCA_FILES[@]}" "${ALIGNMENTS_VDJCA_FILES[@]}" "${ASSEMBLED_CELLS_CLNS_FILES[@]}" "${CLNS_FILES[@]}" "${CLNA_FILES[@]}"; do
    sample=$(basename "$file" | sed 's/\.\(refined\.vdjca\|alignments\.vdjca\|assembledCells\.clns\|clns\|clna\)$//')
    mixcr -Xmx64g exportReports \
        --tables \
        "$file" \
        "${LOG_DIR}/${sample}_report.txt" > "${LOG_DIR}/${sample}_report.stdout.log" 2> "${LOG_DIR}/${sample}_report.stderr.log" || {
        echo "Warning: Failed to generate report for $sample"
    }
done

# Summary
echo "MiXCR analysis and QC completed."
if [ ${#FAILED_SAMPLES[@]} -eq 0 ]; then
    echo "All tasks completed successfully."
else
    echo "Errors occurred for reports: ${FAILED_SAMPLES[*]}"
    exit 1
fi
