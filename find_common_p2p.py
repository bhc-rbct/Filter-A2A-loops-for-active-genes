import argparse
import json
from pathlib import Path
from pybedtools import BedTool

from time import time
import pandas as pd

P2P_SCORE = "p2p2_score"
COL_NAMES_BEDPE_GB = ["chrom1", "start1", "end1",
                      "chrom2", "start2", "end2", P2P_SCORE, "genes"]



def get_interaction_per_enhancer_df(p2p_files, min_score):
    """
    Summarizes the number of interactions per gene, with optional score filtering.

    Args:
        p2p_files (pd.DataFrame): DataFrame of P2P interactions.
        min_score (float): Minimum score threshold for filtering interactions.

    Returns:
        pd.DataFrame: DataFrame summarizing interactions per gene.
    """
    # Apply minimum score filter if required
    if min_score > 1:
        p2p_files = p2p_files[p2p_files[P2P_SCORE] >= min_score]

    # Group by gene, summing scores, and sort by interaction score
    interactions_summary = (
        p2p_files
        .groupby(["TF_chrom", "TF_start", "TF_end"])[["p2p_1", "p2p_2"]]
        .sum()
        .reset_index()
        #.sort_values(by=[P2P_SCORE], ascending=False)
    )

    return interactions_summary

def get_P2P_for_with_TSS(p2p_bedpe, tF_bed):
    """
    Filters P2P interactions for TSS and transcription factors (TFs).

    Args:
        p2p_bedpe (BedTool): P2P BEDPE file.
        tSS_iag_bed (BedTool): BED file containing TSS of genes of interest.
        tF_bed (BedTool): BED file of transcription factor regions.

    Returns:
        BedTool: Filtered BEDPE file for P2P interactions with TSS and TF.
    """
    # Further filter interactions with TF regions that do not overlap TSS
    p2p_for_TSS_iag_tf = p2p_bedpe.pair_to_bed(tF_bed, type="xor")

    return p2p_for_TSS_iag_tf

# Utility Functions
def read_bed_as_df(file, custom_col_name=[], set_numcols=None):
    """
    Reads a BED file into a Pandas DataFrame and sets column names.
    """
    df = pd.read_table(file).dropna(how="all", axis=1)
    num_used_cols = min(set_numcols, len(df.columns)) if set_numcols else len(df.columns)
    col_names = custom_col_name[:num_used_cols] + [
        f"unnamed_{i}" for i in range(num_used_cols - len(custom_col_name))
    ]
    df_cleaned = df.iloc[:, :num_used_cols]
    df_cleaned.columns = col_names
    return df_cleaned


def main(args):
    start_time = time()
    
    p2p_1 = read_bed_as_df(args.p2p_files[0], COL_NAMES_BEDPE_GB).drop('genes', axis=1)
    p2p_2 = read_bed_as_df(args.p2p_files[1], COL_NAMES_BEDPE_GB).drop('genes', axis=1)
    print(f"Number P2Ps of first file: {len(p2p_1)}")
    print(f"Number P2Ps of second file: {len(p2p_2)}")


    # Merge on the first 6 columns (inner join keeps only matching rows)
    p2p_3 = pd.merge(
    p2p_1, 
    p2p_2,
    on=["chrom1", "start1", "end1",
                      "chrom2", "start2", "end2"],
    how='inner',
    suffixes=('_p2p_1', '_p2p_2')  # Renames duplicate columns
    )

    p2p_3.drop(P2P_SCORE + '_p2p_2', axis=1).to_csv(args.output_files[0], sep='\t', header=None, index=None)
    p2p_3.drop(P2P_SCORE + '_p2p_1', axis=1).to_csv(args.output_files[1], sep='\t', header=None, index=None)

    print(f"Number P2Ps of result P2Ps: {len(p2p_3)}")
    
    print(f"--- Total execution time: {time() - start_time:.2f} seconds ---")


if __name__ == "__main__":
    # Argument parsing
    parser = argparse.ArgumentParser(
        description="Performs filtering of a P2P file using a TSS/TES file. Optionally, filters TSS/TES for interesting and active genes based on a gene list, H3K4me3, and TF files.",
        epilog="NOTE: BedTools must be installed!"
    )
    
    parser.add_argument('-p2ps', '--p2p_files', type=str, nargs=2, required=True,
                    help='Paths to the two P2P files (BEDPE format) to find intersection. Provide two file paths separated by a space.')
    parser.add_argument('-o', '--output_files', type=str, nargs=2, required=True,
                        help='Provide two file names separated by a space. The first one is used for the interactions from p2p files one and the second from file two')

    
    args = parser.parse_args()
    
    # Display parameters
    print(f"Parameters:\n- Paths to the two P2P files: {args.p2p_files}\n"
          f"Output Names: {args.output_files}")
    
    # Run the main function
    main(args)
