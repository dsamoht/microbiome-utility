import argparse
import logging
from pathlib import Path
import sys
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("collapse_asv.log", mode="w"),
        logging.StreamHandler()
    ]
)

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("asv_in", help="Input ASV table")
    parser.add_argument("asv_collapsed", help="Collapsed ASV table")
    return parser.parse_args()

def main():
    args = parse_arguments()
    INPUT = Path(args.asv_in)
    OUTPUT = Path(args.asv_collapsed)

    logging.info("Opening ASV table...")
    # rows = samples, columns = sequence strings
    asv_table = pd.read_csv(INPUT, index_col=0, header=0)
    logging.info(f"Loaded table has {asv_table.shape[1]} ASVs and {asv_table.shape[0]} samples.")

    # 1. First, instantly merge 100% exact duplicate sequence columns if any exist
    logging.info("Collapsing exact duplicate sequence strings...")
    asv_table = asv_table.groupby(asv_table.columns, axis=1).sum()
    
    # 2. Sort ASVs by length (descending) so longer sequences can act as targets for shorter sub-sequences
    logging.info("Sorting sequences by length to map fragments...")
    sorted_sequences = sorted(asv_table.columns, key=len, reverse=True)
    
    # Map shorter sequences to longer matching parents
    mapping = {}
    resolved_seqs = [] # Keeps track of validated unique/parent sequences

    for seq in sorted_sequences:
        matched = False
        for parent in resolved_seqs:
            # If the shorter seq is perfectly contained within the longer parent sequence
            if seq in parent: 
                mapping[seq] = parent
                matched = True
                logging.info(f"Mapping fragment sequence (len {len(seq)}) to parent (len {len(parent)})")
                break
        if not matched:
            mapping[seq] = seq
            resolved_seqs.append(seq)

    # 3. Fast vector collapse using pandas groupby mapping
    logging.info("Merging abundances...")
    asv_table = asv_table.groupby(mapping, axis=1).sum()

    logging.info(f"Saving {OUTPUT} ({asv_table.shape[1]} ASVs, {asv_table.shape[0]} samples) to disk.")
    asv_table.to_csv(OUTPUT)
    logging.info("Program finished.")

if __name__ == "__main__":
    main()
