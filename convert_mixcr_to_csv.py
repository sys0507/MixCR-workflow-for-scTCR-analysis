import pandas as pd
import os
import glob
import argparse
import logging
from datetime import datetime

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler('convert_mixcr_to_csv.log')
    ]
)
logger = logging.getLogger(__name__)

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description='Convert MiXCR clone.groups_TRAB.tsv files to a summary CSV.')
    parser.add_argument(
        '--output-dir',
        type=str,
        required=True,
        help='Path to the MiXCR output directory containing clone.groups_TRAB.tsv files'
    )
    return parser.parse_args()

def main():
    """Main function to process MiXCR output and generate CSV."""
    args = parse_arguments()
    output_dir = args.output_dir

    # Validate output directory
    if not os.path.isdir(output_dir):
        logger.error(f"Output directory does not exist: {output_dir}")
        exit(1)

    # Find all clone.groups_TRAB.tsv files
    clone_group_files = glob.glob(os.path.join(output_dir, "results.*.clone.groups_TRAB.tsv"))
    logger.info(f"Found {len(clone_group_files)} clone.groups_TRAB.tsv files")

    if not clone_group_files:
        logger.error(f"No clone.groups_TRAB.tsv files found in {output_dir}")
        exit(1)

    # Initialize list to store data
    data = []

    for clone_group_file in clone_group_files:
        # Extract sample ID from filename
        sample_id = os.path.basename(clone_group_file).replace("results.", "").replace(".clone.groups_TRAB.tsv", "")
        logger.info(f"Processing {clone_group_file} (Sample ID: {sample_id})")

        # Read clone.groups_TRAB.tsv
        try:
            df_clones = pd.read_csv(clone_group_file, sep='\t')
            logger.info(f"  Rows: {len(df_clones)}, Columns: {list(df_clones.columns)}")
        except Exception as e:
            logger.error(f"Failed to read {clone_group_file}: {e}")
            continue

        # Calculate total reads
        total_reads = df_clones['groupReadCount'].sum() if 'groupReadCount' in df_clones.columns else pd.NA

        # Select top clone based on combined TRA/TRB read counts
        if df_clones.empty:
            logger.warning(f"No clonotypes found in {clone_group_file}")
            continue

        required_columns = ['TRA.primary.readCount', 'TRB.primary.readCount']
        if not all(col in df_clones.columns for col in required_columns):
            logger.warning(f"Missing required columns in {clone_group_file}: {required_columns}")
            continue

        df_clones['total_read_count'] = df_clones['TRA.primary.readCount'] + df_clones['TRB.primary.readCount']
        top_clone = df_clones.sort_values(by='total_read_count', ascending=False).head(1)

        # Initialize row dictionary
        row = {'Well': sample_id, 'Num_Total_Reads': total_reads}

        # Extract TRA data
        row['Abundance'] = top_clone['TRA.primary.readCount'].iloc[0] if pd.notna(top_clone['TRA.primary.readCount'].iloc[0]) else pd.NA
        row['CDR3_nt_sequence'] = top_clone['TRA.primary.nSeqCDR3'].iloc[0] if pd.notna(top_clone['TRA.primary.nSeqCDR3'].iloc[0]) else pd.NA
        row['CDR3_aa_sequence'] = top_clone['TRA.primary.aaSeqCDR3'].iloc[0] if pd.notna(top_clone['TRA.primary.aaSeqCDR3'].iloc[0]) else pd.NA
        row['V_Gene'] = top_clone['TRA.primary.bestVHit'].iloc[0] if pd.notna(top_clone['TRA.primary.bestVHit'].iloc[0]) else pd.NA
        row['D_Gene'] = top_clone['TRA.primary.bestDHit'].iloc[0] if pd.notna(top_clone['TRA.primary.bestDHit'].iloc[0]) else pd.NA
        row['J_Gene'] = top_clone['TRA.primary.bestJHit'].iloc[0] if pd.notna(top_clone['TRA.primary.bestJHit'].iloc[0]) else pd.NA

        # Extract TRB data
        row['Abundance.1'] = top_clone['TRB.primary.readCount'].iloc[0] if pd.notna(top_clone['TRB.primary.readCount'].iloc[0]) else pd.NA
        row['CDR3_nt_sequence.1'] = top_clone['TRB.primary.nSeqCDR3'].iloc[0] if pd.notna(top_clone['TRB.primary.nSeqCDR3'].iloc[0]) else pd.NA
        row['CDR3_aa_sequence.1'] = top_clone['TRB.primary.aaSeqCDR3'].iloc[0] if pd.notna(top_clone['TRB.primary.aaSeqCDR3'].iloc[0]) else pd.NA
        row['V_Gene.1'] = top_clone['TRB.primary.bestVHit'].iloc[0] if pd.notna(top_clone['TRB.primary.bestVHit'].iloc[0]) else pd.NA
        row['D_Gene.1'] = top_clone['TRB.primary.bestDHit'].iloc[0] if pd.notna(top_clone['TRB.primary.bestDHit'].iloc[0]) else pd.NA
        row['J_Gene.1'] = top_clone['TRB.primary.bestJHit'].iloc[0] if pd.notna(top_clone['TRB.primary.bestJHit'].iloc[0]) else pd.NA

        data.append(row)
        logger.info(f"Added data for {sample_id}")

    # Check if data is empty
    if not data:
        logger.error("No data processed. Check if files are empty or lack required columns.")
        exit(1)

    # Create dataframe
    df = pd.DataFrame(data)

    # Define column order
    columns = [
        'Well', 'Num_Total_Reads',
        'Abundance', 'CDR3_nt_sequence', 'CDR3_aa_sequence', 'V_Gene', 'D_Gene', 'J_Gene',
        'Abundance.1', 'CDR3_nt_sequence.1', 'CDR3_aa_sequence.1', 'V_Gene.1', 'D_Gene.1', 'J_Gene.1'
    ]
    df = df[columns]

    # Sort by Well
    df = df.sort_values(by='Well')

    # Save to CSV
    output_csv = os.path.join(output_dir, 'mixcr_summary.csv')
    df.to_csv(output_csv, index=False)
    logger.info(f"Saved summary to {output_csv}")

if __name__ == "__main__":
    main()
