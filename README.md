# MixcR Workflow for scTCR-seq Analysis

## Overview
This repository contains a comprehensive workflow for analyzing single-cell T-cell receptor sequencing (scTCR-seq) data using MiXCR, with preprocessing steps for quality control and trimming. The pipeline is designed to run on a high-performance computing (HPC) cluster with SLURM job scheduling, processing raw FASTQ files to generate MiXCR analysis results and detailed quality control (QC) reports.

The workflow includes three main scripts:
1. **runQC.sh**: Performs initial quality control on raw FASTQ files using FastQC and MultiQC.
2. **runQC2.sh**: Trims FASTQ files with fastp, followed by FastQC and MultiQC on trimmed files.
3. **runMixcr.sh**: Runs MiXCR analysis, renames output files, and generates extensive QC reports (align, chainUsage, coverage, tags).

## Features
- **Quality Control**: Uses FastQC and MultiQC to assess raw and trimmed FASTQ file quality.
- **Read Trimming**: Employs fastp for adapter removal and quality-based trimming of paired-end and single-end reads.
- **MiXCR Analysis**: Processes scTCR-seq data with MiXCR to identify and assemble TCR clonotypes, tailored for human samples.
- **Comprehensive QC Reports**: Generates detailed visualizations for alignment, chain usage, coverage, and tags using MiXCR’s exportQc functionality.
- **Robust Error Handling**: Includes extensive checks for file existence, directory permissions, and job failures.
- **SLURM Integration**: Optimized for HPC environments with SLURM job scheduling.

## Prerequisites
- **HPC Cluster**: Access to a SLURM-based HPC system.
- **Modules**:
  - GCC/10.3.0
  - FastQC/0.12.1-Java-11.0.2
  - MultiQC/1.25.1-foss-2021a
  - MiXCR/4.7.0-Java-11.0.2
- **Software**:
  - fastp (downloaded automatically by `runQC2.sh` if not present)
- **Input**: Paired-end or single-end FASTQ files in the specified input directory.
- **Storage**: Sufficient disk space for input, output, and log directories.

## Directory Structure
```
mixcr_workflow/
├── runQC.sh          # Initial QC with FastQC and MultiQC
├── runQC2.sh         # Trimming with fastp, followed by FastQC and MultiQC
├── runMixcr.sh       # MiXCR analysis and QC report generation
├── logs/             # Log files for job outputs and errors
├── data/             # Input FASTQ and output directories (user-defined)
└── README.md         # This file
```

## Usage
1. **Clone the Repository**:
   ```bash
   git clone https://github.com/your-username/mixcr_workflow.git
   cd mixcr_workflow
   ```

2. **Set Up Directories**:
   - Update the directory paths in each script (`INPUT_DIR`, `OUTPUT_DIR`, etc.) to match your file system.
   - Ensure input FASTQ files follow the naming convention: `Jurkat-P0321-Mart1-{{n}}-{{CELL:a}}_{{SAMPLE:a}}_L001_{{R}}_trimmed.fastq.gz`.

3. **Run the Pipeline**:
   ```bash
   chmod +x runQC.sh runQC2.sh runMixcr.sh
   sbatch runQC.sh
   sbatch runQC2.sh
   sbatch runMixcr.sh
   ```

4. **Monitor Jobs**:
   ```bash
   squeue -u $USER
   ```

5. **Check Outputs**:
   - FastQC and MultiQC reports: `$OUTPUT_DIR/FastQC`, `$OUTPUT_DIR/MultiQC`, `$OUTPUT_DIR/FastQC_Trimmed`, `$OUTPUT_DIR/MultiQC_Trimmed`
   - Trimmed FASTQ files: `$OUTPUT_DIR/Trimmed_Fastq`
   - MiXCR results: `$OUTPUT_DIR/data/mixcr_output_scTCR`
   - QC reports: `$OUTPUT_DIR/data/mixcr_output_scTCR/qc_reports`
   - Logs: `$LOG_DIR`

## Notes
- **File Naming**: Ensure FASTQ files match the expected pattern for MiXCR processing. Mismatched names may cause errors.
- **Resource Allocation**: Adjust SLURM parameters (`--ntasks`, `--mem`, `--time`) based on your cluster’s configuration and dataset size.
- **Error Logs**: Check log files in `$LOG_DIR` for troubleshooting.
- **Expected Samples**: The `runMixcr.sh` script includes a check for 5 samples; modify the `EXPECTED_SAMPLES` variable if needed.

## Contributing
Contributions are welcome! Please submit a pull request or open an issue for bug reports, feature requests, or suggestions.

## License
This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## References
- [MiXCR Documentation](https://mixcr.com/mixcr/)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [MultiQC](https://multiqc.info/)
- [fastp](https://github.com/OpenGene/fastp)
