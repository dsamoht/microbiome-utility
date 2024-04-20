"""
usage:
python collapse_asv.py [SEQTAB] [COLLAPSED_SEQTAB]
"""
import argparse
from collections import defaultdict
from datetime import datetime
from glob import glob
import logging
from pathlib import Path
import shutil
import subprocess
import sys

import pandas as pd


__author__ = "deschenes.thomas@gmail.com"


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("collapse_asv.log".format(datetime.now()), mode="w"),
        logging.StreamHandler()
    ]
)

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("asv_in", help="Input ASV table")
    parser.add_argument("asv_collapsed", help="Collapsed ASV table")
    return parser.parse_args()


def main():

    INPUT = Path(parse_arguments().asv_in)
    OUTPUT = Path(parse_arguments().asv_collapsed)
    SESSION_ID = "tmp-collapse-asv"

    logging.info("Checking if `cd-hit` is installed and accessible...")
    if not shutil.which('cd-hit'):
        logging.error("`cd-hit` (https://github.com/weizhongli/cdhit) does not seem to be installed. Please install it and make it accessible in your PATH before running the script.")
        sys.exit(1)
    logging.info("Done.")

    # Open ASV table. A seqtab matrix from DADA2 pipeline. rows are samples; columns are ASVs.
    logging.info("Opening ASV table...")
    asv_table = pd.read_csv(INPUT, index_col=0, header=0)
    asv_id_to_seq = {f"asv_{i}": asv for i, asv in enumerate(asv_table.columns)}
    logging.info("Done.")
    logging.info(f"Loaded table has {asv_table.shape[1]} ASVs and {asv_table.shape[0]} samples.")

    with open(f"{SESSION_ID}.fna", "w") as fna_out:
        for i, asv in enumerate(asv_table.columns):
            fna_out.write(f">asv_{i}\n{asv}\n")


    logging.info("Running `cd-hit-est` to cluster ASVs at 100% sequence identity...")
    try:
        subprocess.run(["cd-hit-est",
                        "-i", f"{SESSION_ID}.fna",
                        "-o", SESSION_ID,
                        "-T", "0",
                        "-M", "0",
                        "-c", "1",
                        "-d", "0",
                        "-g", "1"],
                        check=True)
    except subprocess.CalledProcessError:
        logging.error("Clustering failed. Is your input ASV table in the good format ?")
        sys.exit()

    logging.info("Done.")

    # Create 1 list of ASVs per cluster
    asv_groups = defaultdict(list)
    with open(f"{SESSION_ID}.clstr", "r") as cdhit_out:
        for line in cdhit_out:
            if line.strip().startswith(">"):
                current_cluster = "cluster_" + line.strip().split()[1]
            else:
                asv = line.strip().split()[2].split(">")[1].split(r"...")[0]
                asv_groups[current_cluster].append(asv)
    
    # Collapse identical ASVs
    # The representative ASV is the most abundant accross samples
    for asv_list in asv_groups.values():
        if len(asv_list) > 1:
            most_abundant_val = 0
            most_abundant_seq = ""
            for asv in asv_list:
                asv_abundance = asv_table[asv_id_to_seq[asv]].sum()
                if asv_abundance > most_abundant_val:
                    most_abundant_val = asv_abundance
                    most_abundant_seq = asv_id_to_seq[asv]

            asv_to_remove = [asv_id_to_seq[asv] for asv in asv_list if most_abundant_seq != asv_id_to_seq[asv]]
            for asv in asv_to_remove:
                asv_table[most_abundant_seq] = asv_table[most_abundant_seq] + asv_table[asv]
                logging.info(f"ASV `{asv}` was collapsed into ASV `{most_abundant_seq}`")
            asv_table = asv_table.drop(asv_to_remove, axis=1)

    asv_table.to_csv(OUTPUT)
    for tmp_file in glob(f"{SESSION_ID}*"):
        Path(tmp_file).unlink()

    logging.info(f"Saving {OUTPUT} ({asv_table.shape[1]} ASVs, {asv_table.shape[0]} samples) to disk.")
    logging.info("Program finished.")


if __name__ == "__main__":
    main()
