#!/usr/bin/env python3

"""
Author: Nadine Schacherer
Last Modified: 10.04.2025 - 11:46
Description: Script to filter P2P files using TSS/TES data with optional gene filtering.
"""

import argparse
import json
from typing import List
import pandas as pd
from pathlib import Path
import os
from pybedtools import BedTool
from time import time
# import openpyxl

# Intermediate file names
INT_P2P_NAME = '_for_TSS_iag.bedpe'

# Colnames for bedpe file if associated TS --> Change if bedpe files consist of other names
TES_TES_OFFSET = 9
COL_NAMES_GENE = "TSS_gene_name"
P2P_SCORE = "p2p2_score"
COL_NAMES_BEDPE = ["chrom1", "start1", "end1", "chrom2", "start2", "end2", P2P_SCORE, "TSS_chrom",
                   "TSS_start", "TSS_end", COL_NAMES_GENE, "TF_chrom", "TF_start", "TF_end"]
COL_NAMES_BEDPE_OLD = ["chrom1", "start1", "end1", "chrom2", "start2", "end2", P2P_SCORE, "TSS_chrom",
                   "TSS_start", "TSS_end", COL_NAMES_GENE, "TSS_score", "TSS_strand", "TF_chrom", "TF_start", "TF_end"]
COL_NAMES_BED_INTERSECT = ["TSS_chrom", "TSS_start", "TSS_end", COL_NAMES_GENE,
                           "TSS_score", "TSS_strand", "H3K4me3_chrom", "H3K4me3_start", "H3K4me3_end"]
COL_NAMES_BED_INTERSECT_KEEP = ["H3K4me3_chrom",
                                "H3K4me3_start", "H3K4me3_end", COL_NAMES_GENE]
COL_NAMES_BEDPE_GB = ["chrom1", "start1", "end1",
                      "chrom2", "start2", "end2", P2P_SCORE]
COL_SIMP_BED = ["chrom", "start", "end"]
COL_NAME_STRAND = "TSS_strand"

BEDPE_COL_NAMES_G = ["chrom1", "start1", "end1", "chrom2",
                     "start2", "end2", P2P_SCORE, "strand1", "strand2"]
BED_COL_NAMES_G = ["chrom", "chromStart", "chromEnd", "name", "score", "strand",
                   "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStart"]

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


def provide_TSS_or_TES_per_gene(row, ts_type):
    """
    Adjusts start and end positions for TSS or TES based on strand and offset.
    """
    if ts_type == "TSS":
        if row["strand"] == "+":
            row["chromStart"] = max(row["chromStart"] - TES_TES_OFFSET, 0)
            row["chromEnd"] = row["chromStart"] + TES_TES_OFFSET
        else:
            row["chromEnd"] = max(row["chromEnd"] - TES_TES_OFFSET, 0)
            row["chromStart"] = row["chromEnd"] - TES_TES_OFFSET
    elif ts_type == "TES":
        if row["strand"] == "+":
            row["chromEnd"] = row["chromEnd"] + TES_TES_OFFSET
            row["chromStart"] = max(row["chromEnd"] - TES_TES_OFFSET, 0)
        else:
            row["chromStart"] = row["chromStart"] + TES_TES_OFFSET
            row["chromEnd"] = max(row["chromStart"] - TES_TES_OFFSET, 0)
    return row


def get_TSS_for_interesting_genes(gene_list, tSS_TSE_df, ts_kind):
    """
    Filters TSS/TES for a list of genes and recalculates positions.
    """
    tSS_TSE_ip_df = (
        tSS_TSE_df.loc[tSS_TSE_df["name"].isin(gene_list)]
        if gene_list
        else tSS_TSE_df
    )
    if ts_kind != "both":
        tSS_TSE_ip_df = tSS_TSE_ip_df.apply(
            lambda row: provide_TSS_or_TES_per_gene(row, ts_kind), axis=1
        )
    return tSS_TSE_ip_df


def get_TSS_for_active_genes(tSS_TSE_ig_bed, h3k4me3_bed):
    """
    Identifies active genes by intersecting TSS with H3K4me3 peaks.

    Args:
        tSS_TSE_ig_bed (BedTool): TSS/TES BED file for genes of interest.
        h3k4me3_bed (BedTool): H3K4me3 BED file.

    Returns:
        BedTool: Filtered BED file containing TSS of active genes.
    """
    return tSS_TSE_ig_bed.intersect(h3k4me3_bed, wb=True)


def get_P2P_for_with_TSS(p2p_bedpe, tSS_iag_bed, tF_bed):
    """
    Filters P2P interactions for TSS and transcription factors (TFs).

    Args:
        p2p_bedpe (BedTool): P2P BEDPE file.
        tSS_iag_bed (BedTool): BED file containing TSS of genes of interest.
        tF_bed (BedTool): BED file of transcription factor regions.

    Returns:
        BedTool: Filtered BEDPE file for P2P interactions with TSS and TF.
    """
    # Identify TF regions that do not intersect with TSS
    tf_wo_tss = tF_bed.intersect(tSS_iag_bed, v=True, wa=True)

    # Filter P2P interactions with TSS of genes of interest
    p2p_for_TSS_iag = p2p_bedpe.pair_to_bed(tSS_iag_bed, type="xor")

    # Further filter interactions with TF regions that do not overlap TSS
    p2p_for_TSS_iag_tf = p2p_for_TSS_iag.pair_to_bed(tf_wo_tss, type="xor")

    return p2p_for_TSS_iag_tf


def agg_unique_p2p_plus_filter(p2p_for_TSS_iag_df, min_score):
    """
    Aggregates unique P2P interactions, associates genes, and filters by minimum score.

    Args:
        p2p_for_TSS_iag_df (pd.DataFrame): DataFrame of P2P interactions.
        min_score (float): Minimum score threshold for filtering interactions.

    Returns:
        pd.DataFrame: Aggregated and filtered P2P interactions.
    """
    # Group by P2P columns, aggregate unique genes
    grouped_df = (
        p2p_for_TSS_iag_df
        .groupby(COL_NAMES_BEDPE_GB)[COL_NAMES_GENE]
        .agg(lambda x: ', '.join(x.unique()))
        .reset_index()
    )

    return grouped_df


def get_is_per_gene_df(p2p_for_TSS_iag_df, min_score):
    """
    Summarizes the number of interactions per gene, with optional score filtering.

    Args:
        p2p_for_TSS_iag_df (pd.DataFrame): DataFrame of P2P interactions.
        min_score (float): Minimum score threshold for filtering interactions.

    Returns:
        pd.DataFrame: DataFrame summarizing interactions per gene.
    """
    # Apply minimum score filter if required
    if min_score > 1:
        p2p_for_TSS_iag_df = p2p_for_TSS_iag_df[p2p_for_TSS_iag_df[P2P_SCORE] >= min_score]

    # Group by gene, summing scores, and sort by interaction score
    interactions_summary = (
        p2p_for_TSS_iag_df
        .groupby(COL_NAMES_GENE)[P2P_SCORE]
        .sum()
        .reset_index()
        .sort_values(by=[P2P_SCORE], ascending=False)
    )

    return interactions_summary

def get_interaction_per_enhancer_df(p2p_for_TSS_iag_df, min_score):
    """
    Summarizes the number of interactions per gene, with optional score filtering.

    Args:
        p2p_for_TSS_iag_df (pd.DataFrame): DataFrame of P2P interactions.
        min_score (float): Minimum score threshold for filtering interactions.

    Returns:
        pd.DataFrame: DataFrame summarizing interactions per gene.
    """
    # Apply minimum score filter if required
    if min_score > 1:
        p2p_for_TSS_iag_df = p2p_for_TSS_iag_df[p2p_for_TSS_iag_df[P2P_SCORE] >= min_score]

    # First remove duplicate TF-TSS pairs keeping the first P2P_SCORE (or sum them if needed)
    unique_scores = p2p_for_TSS_iag_df.drop_duplicates(["chrom1", "start1", "end1", "chrom2", "start2", "end2", "TSS_chrom", "TSS_start", "TSS_end"])

    # Then group by TF coordinates
    sum_scores = unique_scores.groupby(['TF_chrom', 'TF_start', 'TF_end'])[P2P_SCORE].sum().reset_index()
    all_genes = p2p_for_TSS_iag_df.groupby(['TF_chrom', 'TF_start', 'TF_end'])[COL_NAMES_GENE].agg(lambda x: ', '.join(x.unique())).reset_index()

    # Merge the results
    interactions_summary = pd.merge(sum_scores, all_genes, on=['TF_chrom', 'TF_start', 'TF_end'])

    return interactions_summary


def runAllSteps(
    genes_list: List[str],
    tSS_TSE_df: pd.DataFrame,
    h3k4me3_bed: BedTool,
    tF_bed: BedTool,
    p2p_bedpe: BedTool,
    ts_kind: str,
    min_P2P_SCORE: int,
    resFPrefix: str,
    TS_fPrefix: str,
    output_dir: str
):
    """
    Runs the full pipeline for filtering P2P interactions based on TSS/TES, histone marks, and transcription factors.

    Args:
        genes_list (List[str]): List of genes of interest.
        tSS_TSE_df (pd.DataFrame): DataFrame of TSS/TES data.
        h3k4me3_bed (BedTool): BedTool object for H3K4me3 regions.
        tF_bed (BedTool): BedTool object for transcription factor regions.
        p2p_bedpe (BedTool): BedTool object for P2P interactions.
        ts_kind (str): Type of transcription site to consider ('TSS', 'TES', or 'both').
        min_P2P_SCORE (int): Minimum score threshold for filtering P2P interactions.
        resFPrefix (str): Prefix for result files.
        TS_fPrefix (str): Prefix for transcription site-related files.
        output_dir (str): Directory for saving output files.

    Returns:
        None
    """
    # Step 1: Process TSS/TES filtering if no pre-filtered file is provided

    print('-------- STEP 1: Get TSS/TES of interesting genes ------------')
    tss_ip_df = get_TSS_for_interesting_genes(genes_list, tSS_TSE_df, ts_kind)
    print(f"Reduced from {len(tSS_TSE_df)} to {len(tss_ip_df)} TSS/TES")

    # Save filtered TSS/TES file for BedTool processing
    tSS_TSE_df_ig_path = output_dir / f"{TS_fPrefix}_ig.bed"
    tss_ip_df.to_csv(tSS_TSE_df_ig_path, sep='\t', header=None, index=None)
    tSS_TSE_ig_bed = BedTool(tSS_TSE_df_ig_path)

    # Step 2: Filter TSS/TES for active genes using histone marks
    print('-------- STEP 2: Get TSS/TES of active genes ------------')
    tSS_iag_bed = get_TSS_for_active_genes(tSS_TSE_ig_bed, h3k4me3_bed)
    print(f"Reduced from {len(tSS_TSE_ig_bed)} to {len(tSS_iag_bed)} TSS/TES (with duplicates)")

    # Save filtered file, remove duplicates, and update BedTool object
    tSS_iag_bed.saveas(output_dir / f"{TS_fPrefix}_iag.bed")
    tss_iag_df_right_format = (
        read_bed_as_df(output_dir / f"{TS_fPrefix}_iag.bed", custom_col_name=COL_NAMES_BED_INTERSECT)
        [COL_NAMES_BED_INTERSECT_KEEP]
        .drop_duplicates()
    )
    tss_iag_df_right_format.to_csv(
        output_dir / f"{TS_fPrefix}_iag.bed", header=None, index=None, sep='\t'
    )

    print(f"Reduced from {len(tSS_TSE_ig_bed)} to {len(tss_iag_df_right_format)} TSS/TES (without duplicates)")
    tSS_iag_bed = BedTool(output_dir / f"{TS_fPrefix}_iag.bed")


    # Step 3: Filter P2P interactions based on filtered TSS/TES
    print('-------- STEP 3: P2P filtering ------------')
    p2p_for_TSS_iag = get_P2P_for_with_TSS(p2p_bedpe, tSS_iag_bed, tF_bed)
    p2p_for_TSS_iag.saveas(output_dir / f"{resFPrefix}{INT_P2P_NAME}")

    # Load P2P interactions into a DataFrame for further processing
    p2p_for_TSS_iag_df = read_bed_as_df(
        output_dir / f"{resFPrefix}{INT_P2P_NAME}", custom_col_name=COL_NAMES_BEDPE
    )


    # Step 4: Aggregate and filter interactions
    is_per_gene_df = get_is_per_gene_df(p2p_for_TSS_iag_df, min_P2P_SCORE)
    is_per_gene_df.to_csv(
        output_dir / f"interaction_per_iag.txt",
        sep='\t', header=None, index=None
    )

    # Step 5: Aggregate and filter interactions
    interaction_per_enhencer_region = get_interaction_per_enhancer_df(p2p_for_TSS_iag_df, min_P2P_SCORE)
    interaction_per_enhencer_region.to_csv(
        output_dir / f"interaction_per_enhancer_region_iag.bed",
        sep='\t', header=None, index=None
    )

    p2p_for_TSS_iag_df_g = agg_unique_p2p_plus_filter(p2p_for_TSS_iag_df, min_P2P_SCORE)
    p2p_for_TSS_iag_df_g.to_csv(
        output_dir / f"p2p_plus_numInteractions_plus_associated_iaGeneList.bedpe",
        sep='\t', header=False, index=False
    )
    print(f"Reduced from {len(p2p_bedpe)} to {len(p2p_for_TSS_iag_df_g)} P2P interactions")

    # Save active genes and clean up intermediate files
    active_genes = is_per_gene_df.iloc[:, 0].tolist()
    TSS_active_genes = tSS_TSE_df[tSS_TSE_df.iloc[:, 3].isin(active_genes)].drop_duplicates(subset=[tSS_TSE_df.columns[3]], keep='first')
    TSS_active_genes.to_csv(output_dir / 'TSS_TES_for_iag_only_first_appearance_per_gene.bed', sep='\t', header=False, index=False)

    # os.remove(output_dir / f"{resFPrefix}{INT_P2P_NAME}")
    os.remove(tSS_TSE_df_ig_path)
    os.remove(output_dir / f"{TS_fPrefix}_iag.bed")



def p2p_filter_main(args):
    start_time = time()
    
    # Prepare output directory
    output_directory = Path(args.output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)
    
    # Save parameters to a file for reproducibility
    param_file = output_directory / "parameters_used.json"
    with open(param_file, "w") as f:
        json.dump(vars(args), f, indent=4)
    
    # Read input gene list if provided
    i_genes_list = None
    if args.genes_filename:
        with open(args.genes_filename, 'r') as f:
            i_genes_list = f.read().splitlines()

    # Load input data
    tSS_TSE_df = read_bed_as_df(args.tSS_TSE_uf, BED_COL_NAMES_G) if args.tSS_TSE_uf else None
    h3k4me3_bed = BedTool(args.histone_bed_file) if args.histone_bed_file else None

    # Handle transcription factor (TF) file
    tf_basename = os.path.basename(args.tF_bed_file)
    tf_temp_path = output_directory / f"{os.path.splitext(tf_basename)[0]}_temp.bed"
    tf_df = read_bed_as_df(args.tF_bed_file, custom_col_name=COL_SIMP_BED, set_numcols=3)
    tf_df.to_csv(tf_temp_path, sep='\t', header=None, index=None)
    tF_bed = BedTool(tf_temp_path)

    # Load P2P file
    p2p_temp_path = output_directory / f"{os.path.splitext(tf_basename)[0]}_temp.bedpe"
    tf_df = read_bed_as_df(args.p2p_file, custom_col_name=COL_NAMES_BEDPE_GB, set_numcols=7)
    tf_df.to_csv(p2p_temp_path, sep='\t', header=None, index=None)
    p2p_bedpe = BedTool(p2p_temp_path)

    # Generate file prefixes for output
    resFPrefix = f"{os.path.basename(args.p2p_file).split('.bedpe')[0].replace('A2A', '')}x_{os.path.basename(args.tF_bed_file).split('.bed')[0]}"
    TS_fPrefix = (f"{os.path.basename(args.tSS_TSE_uf).split('.bed')[0]}_c{args.kind_TS}_"
                  f"{os.path.basename(args.histone_bed_file).split('.bed')[0]}" if args.tSS_TSE_uf else None)

    # Run main steps
    runAllSteps(i_genes_list, tSS_TSE_df, h3k4me3_bed, tF_bed, p2p_bedpe, args.kind_TS,
                args.min_P2P_SCORE, resFPrefix, TS_fPrefix, output_directory)
    
    # Clean up temporary files
    os.remove(tf_temp_path)
    os.remove(p2p_temp_path)

    
    print(f"--- Total execution time: {time() - start_time:.2f} seconds ---")


if __name__ == "__main__":
    # Argument parsing
    parser = argparse.ArgumentParser(
        description="Performs filtering of a P2P file using a TSS/TES file. Optionally, filters TSS/TES for interesting and active genes based on a gene list, H3K4me3, and TF files.",
        epilog="NOTE: BedTools must be installed!"
    )
    
    parser.add_argument('-g', '--genes_filename', type=str, required=False,
                        help='Path to a text file containing gene names of interest')
    parser.add_argument('-uts', '--tSS_TSE_uf', type=str, required=False,
                        help='Path to the unfiltered TSS/TES file (BED format)')
    parser.add_argument('-hf', '--histone_bed_file', type=str, required=False,
                        help='Path to the Histone BED file for filtering')
    parser.add_argument('-p2p', '--p2p_file', type=str, required=True,
                        help='Path to the P2P file (BEDPE format)')
    parser.add_argument('-tff', '--tF_bed_file', type=str, required=True,
                        help='Path to the Transcription Factor BED file')
    # parser.add_argument('-fts', '--tSS_TSE_f', type=str, required=False,
    #                     help='Path to the filtered TSS/TES file (BED format)')
    parser.add_argument('-o', '--output_directory', type=str, default='./', required=False,
                        help='Path to the output directory (default: current directory)')
    parser.add_argument('-mS', '--min_P2P_SCORE', type=int, default=1, required=False,
                        help='Minimum score to retain P2P interactions (default: 1)')
    parser.add_argument('-kindTS', '--kind_TS', type=str, choices=['both', 'TSS', 'TES'], default='TSS', required=False,
                        help="Type of Transcription Sites to use (default: 'TSS')")
    
    args = parser.parse_args()
    
    # Display parameters
    print(f"Parameters:\n- Lowest score to keep P2P: {args.min_P2P_SCORE}\n"
          f"- TS type: {args.kind_TS}\n- Output Directory: {args.output_directory}")
    
    # Run the main function
    p2p_filter_main(args)

