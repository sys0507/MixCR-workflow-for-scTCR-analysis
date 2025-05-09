# MixcR Workflow for scTCR-seq Analysis

## Overview
This repository contains a comprehensive workflow for analyzing single-cell T-cell receptor sequencing (scTCR-seq) data using MiXCR, with preprocessing steps for quality control, trimming, and post-processing to convert MiXCR outputs to CSV. The pipeline is designed to run on a high-performance computing (HPC) cluster with SLURM job scheduling, processing raw FASTQ files to generate MiXCR analysis results, detailed quality control (QC) reports, and a summarized CSV of clonotype data.

The workflow includes four main scripts:
1. **runQC.sh**: Performs initial quality control on raw FASTQ files using FastQC and MultiQC.
2. **runQC2.sh**: Trims FASTQ files with fastp, followed by FastQC and MultiQC on trimmed files.
3. **runMixcr.sh**: Runs MiXCR analysis, renames output files, and generates extensive QC reports (align, chainUsage, coverage, tags).
4. **convert_mixcr_to_csv.py**: Converts MiXCR `clone.groups_TRAB.tsv` files into a summarized CSV containing TRA and TRB clonotype data for the top clone per sample.

## Features
- **Quality Control**: Uses FastQC and MultiQC to assess raw and trimmed FASTQ file quality.
- **Read Trimming**: Employs fastp for adapter removal and quality-based trimming of paired-end and single-end reads.
- **MiXCR Analysis**: Processes scTCR-seq data with MiXCR to identify and assemble TCR clonotypes, tailored for human samples.
- **Comprehensive QC Reports**: Generates detailed visualizations for alignment, chain usage, coverage, and tags using MiXCR’s exportQc functionality.
- **Post-Processing**: Converts MiXCR output to a CSV summarizing key TRA and TRB clonotype metrics (abundance, CDR3 sequences, V/D/J genes).
- **Robust Error Handling**: Includes extensive checks for file existence, directory permissions, and job failures.
- **SLURM Integration**: Optimized for HPC environments with SLURM job scheduling.

## Prerequisites
- **HPC Cluster**: Access to a SLURM-based HPC system (for shell scripts).
- **Modules** (for shell scripts):
  - GCC/10.3.0
  - FastQC/0.12.1-Java-11.0.2
  - MultiQC/1.25.1-foss-2021a
  - MiXCR/4.7.0-Java-11.0.2
- **Software**:
  - fastp (downloaded automatically by `runQC2.sh` if not present)
  - Python 3.6+ (for `convert_mixcr_to_csv.py`)
  - Python packages: `pandas`, `argparse`, `logging` (install via `pip install pandas`)
- **Input**: Paired-end or single-end FASTQ files in the specified input directory.
- **Storage**: Sufficient disk space for input, output, and log directories.

## Directory Structure
```
mixcr_workflow/
├── runQC.sh                  # Initial QC with FastQC and MultiQC
├── runQC2.sh                 # Trimming with fastp, followed by FastQC and MultiQC
├── runMixcr.sh               # MiXCR analysis and QC report generation
├── convert_mixcr_to_csv.py   # Convert MiXCR output to CSV
├── logs/                     # Log files for job outputs and errors
├── data/                     # Input FASTQ and output directories (user-defined)
└── README.md                 # This file
```

## Usage
1. **Clone the Repository**:
   ```bash
   git clone https://github.com/your-username/mixcr_workflow.git
   cd mixcr_workflow
   ```

2. **Set Up Directories**:
   - Update the directory paths in each shell script (`INPUT_DIR`, `OUTPUT_DIR`, etc.) to match your file system.
   - Ensure input FASTQ files follow the naming convention: `Jurkat-P0321-Mart1-{{n}}-{{CELL:a}}_{{SAMPLE:a}}_L001_{{R}}_trimmed.fastq.gz`.
   - For the Python script, specify the MiXCR output directory containing `clone.groups_TRAB.tsv` files.

3. **Run the Pipeline**:
   - **Shell Scripts**:
     ```bash
     chmod +x runQC.sh runQC2.sh runMixcr.sh
     sbatch runQC.sh
     sbatch runQC2.sh
     sbatch runMixcr.sh
     ```
   - **Python Script**:
     ```bash
     python convert_mixcr_to_csv.py --output-dir /path/to/mixcr_output_scTCR
     ```
     Replace `/path/to/mixcr_output_scTCR` with the path to the directory containing MiXCR `clone.groups_TRAB.tsv` files.

4. **Monitor Jobs** (for shell scripts):
   ```bash
   squeue -u $USER
   ```

5. **Check Outputs**:
   - FastQC and MultiQC reports: `$OUTPUT_DIR/FastQC`, `$OUTPUT_DIR/MultiQC`, `$OUTPUT_DIR/FastQC_Trimmed`, `$OUTPUT_DIR/MultiQC_Trimmed`
   - Trimmed FASTQ files: `$OUTPUT_DIR/Trimmed_Fastq`
   - MiXCR results: `$OUTPUT_DIR/data/mixcr_output_scTCR`
   - QC reports: `$OUTPUT_DIR/data/mixcr_output_scTCR/qc_reports`
   - Summary CSV: `$OUTPUT_DIR/data/mixcr_output_scTCR/mixcr_summary.csv`
   - Logs: `$LOG_DIR` (for shell scripts) and `convert_mixcr_to_csv.log` (for Python script)

## Notes
- **File Naming**: Ensure FASTQ files match the expected pattern for MiXCR processing, and `clone.groups_TRAB.tsv` files follow the `results.*.clone.groups_TRAB.tsv` naming convention.
- **Resource Allocation**: Adjust SLURM parameters (`--ntasks`, `--mem`, `--time`) based on your cluster’s configuration and dataset size.
- **Error Logs**: Check log files in `$LOG_DIR` (shell scripts) or `convert_mixcr_to_csv.log` (Python script) for troubleshooting.
- **Expected Samples**: The `runMixcr.sh` script includes a check for 5 samples; modify the `EXPECTED_SAMPLES` variable if needed.
- **Python Environment**: Ensure Python and required packages are installed. On an HPC, you may need to load a Python module or use a virtual environment.

## Contributing
Contributions are welcome! Please submit a pull request or open an issue for bug reports, feature requests, or suggestions.

## License
This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## References
- [MiXCR Documentation](https://mixcr.com/mixcr/)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [MultiQC](https://multiqc.info/)
- [fastp](https://github.com/OpenGene/fastp)
- [Pandas](https://pandas.pydata.org/)
