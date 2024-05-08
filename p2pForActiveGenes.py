#!/usr/bin/env python3

"""
author: Nadine Schacherer
last modified: 23.05.2023
"""

import argparse
from typing import List
import pandas as pd
from pathlib import Path
import os
from pybedtools import BedTool
from time import time
# import openpyxl

# Intermediate file names
int_p2P_name = '_for_TSS_iag.bedpe'

# Colnames for bedpe file if associated TS --> Change if bedpe files consist of other names
TES_TES_offset = 9
col_names_gene = "TSS_gene_name"
p2p_score = "p2p2_score"
col_names_bedpe = ["chrom1", "start1", "end1", "chrom2", "start2", "end2", p2p_score, "TSS_chrom",
                   "TSS_start", "TSS_end", "TSS_gene_name", "TSS_score", "TSS_strand", "TF_chrom", "TF_start", "TF_end"]
col_names_bed_intersect = ["TSS_chrom", "TSS_start", "TSS_end", col_names_gene,
                           "TSS_score", "TSS_strand", "H3K4me3_chrom", "H3K4me3_start", "H3K4me3_end"]
col_names_bed_intersect_keep = ["H3K4me3_chrom",
                                "H3K4me3_start", "H3K4me3_end", col_names_gene]
col_names_bedpe_gb = ["chrom1", "start1", "end1",
                      "chrom2", "start2", "end2", p2p_score]
col_simp_bed = ["chrom", "start", "end"]
col_name_strand = "TSS_strand"

bedpe_col_names_g = ["chrom1", "start1", "end1", "chrom2",
                     "start2", "end2", p2p_score, "strand1", "strand2"]
bed_col_names_g = ["chrom", "chromStart", "chromEnd", "name", "score", "strand",
                   "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStart"]


def read_bed_as_df(file, custom_col_name=[], set_numcols=None):
    """ read in bed file as panda df; set column names and data types for columns
    Args:
        file_path (_type_): 
    """

    df = pd.read_table(file).dropna(how='all', axis=1, inplace=False)
    num_used_cols = min(set_numcols, len(df.columns)) if set_numcols is not None else len(df.columns)
    df_cleaned = df.iloc[:, :num_used_cols]
    df_cleaned.columns = custom_col_name[:num_used_cols] + ["unnamed_"+str(i) for i in range(num_used_cols - len(custom_col_name))]
    return df_cleaned


def provide_TSS_or_TES_per_gene(row, ts_type):
    # calculated new start and end positions depending on the strand and the required TS (also using a offset to account for experimental inaccuracies)
    if ts_type == 'TSS':
        row["chromStart"] = max(row["chromStart"] -
                                TES_TES_offset, 0) if row['strand'] == '+' else max(row["chromEnd"] -
                                                                                    TES_TES_offset, 0)
        row["chromEnd"] = row["chromStart"] + \
            TES_TES_offset if row['strand'] == '+' else row["chromEnd"] + \
            TES_TES_offset
    elif ts_type == 'TES':
        row["chromStart"] = max(row["chromEnd"] -
                                TES_TES_offset, 0) if row['strand'] == '+' else max(row["chromStart"] -
                                                                                    TES_TES_offset, 0)
        row["chromEnd"] = row["chromEnd"] + \
            TES_TES_offset if row['strand'] == '+' else row["chromStart"] + \
            TES_TES_offset
    return row


def get_TSS_for_interesting_genes(gene_list, tSS_TSE_df, ts_kind):
    # filter for TSS/TES of genes supplied in list and recalculate start & end position of required TSS and or TES

    tSS_TSE_ip_df = tSS_TSE_df.loc[tSS_TSE_df["name"].isin(
        gene_list)] if gene_list is not None else tSS_TSE_df

    if ts_kind != "both":
        tSS_TSE_ip_df = tSS_TSE_ip_df.apply(
            lambda row: provide_TSS_or_TES_per_gene(row, ts_kind), axis=1)

    return tSS_TSE_ip_df


def get_TSS_for_active_genes(tSS_TSE_ig_bed, h3k4me3_bed):

    tSS_iag_bed = tSS_TSE_ig_bed.intersect(h3k4me3_bed, wb=True)

    return tSS_iag_bed


def get_P2P_for_with_TSS(p2p_bedpe, tSS_iag_bed, tF_bed):

    tf_wo_tss = tF_bed.intersect(tSS_iag_bed, v=True, wa=True)
    p2p_for_TSS_iag = p2p_bedpe.pair_to_bed(tSS_iag_bed, type="xor")
    p2p_for_TSS_iag_tf = p2p_for_TSS_iag.pair_to_bed(tf_wo_tss, type="xor")

    return p2p_for_TSS_iag_tf


def agg_unique_p2p_plus_filter(p2p_for_TSS_iag_df, min_score):
    # Group by P2P, save the list of genes associated + if wanted filter for min_score
    p2p_for_TSS_iag_df_g = p2p_for_TSS_iag_df.groupby(
        col_names_bedpe_gb)[col_names_gene].unique().reset_index()
    if min_score > 1:
        p2p_for_TSS_iag_df_g = p2p_for_TSS_iag_df_g[p2p_for_TSS_iag_df_g[p2p_score] >= min_score]
    return p2p_for_TSS_iag_df_g


def get_is_per_gene_df(p2p_for_TSS_iag_df, min_score):
    # summerize number of interactions per gene
    if min_score > 1:
        p2p_for_TSS_iag_df = p2p_for_TSS_iag_df[p2p_for_TSS_iag_df[p2p_score] >= min_score]

    p2p_for_TSS_iag_df_is = p2p_for_TSS_iag_df.groupby(col_names_gene)[p2p_score].sum(
    ).reset_index().sort_values(by=[p2p_score], ascending=False)
    return p2p_for_TSS_iag_df_is


def runAllSteps(genes_list: List[str], tSS_TSE_df: pd.DataFrame, h3k4me3_bed: BedTool, tF_bed: BedTool, p2p_bedpe: BedTool, ts_kind: str, min_p2p_score: int, ts_iag_file: str, resFPrefix: str, TS_fPrefix: str, output_dir: str):

    # Perform filtering of Transcription sides (TS) only if wanted --> ts_iag_file needes to be supplied by user

    if ts_iag_file is None:
        # Perform first TS filtering Step --> filter for TS of genes one is interested in (user can also choose if he wants only TSS, TES or both)
        print('-------- STEP 1: Get TSS/TES of interesting genes ------------')
        tss_ip_df = get_TSS_for_interesting_genes(
            genes_list, tSS_TSE_df, ts_kind)
        print("Reduced from %i to %i TSS/TES" %
              (len(tSS_TSE_df.index), len(tss_ip_df.index)))

        # Save filtered TS file as bed to be able to load it with BedTools
        tSS_TSE_df_ig_path = output_dir / (TS_fPrefix + '_ig.bed')
        print(tSS_TSE_df_ig_path, TS_fPrefix)
        tss_ip_df.to_csv(tSS_TSE_df_ig_path, sep='\t', header=None, index=None)
        tSS_TSE_ig_bed = BedTool(tSS_TSE_df_ig_path)

        # Perform second TS filtering Step --> filter for active genes using the supplied h3k4me3 and tf file using intersect from BedTools
        print('-------- STEP 2: Get TSS/TES of active genes ------------')
        tSS_iag_bed = get_TSS_for_active_genes(tSS_TSE_ig_bed, h3k4me3_bed)
        print("Reduced from %i to %i TSS/TES (with duplicates)" %
              (len(tSS_TSE_ig_bed), len(tSS_iag_bed)))

        # Save filtered TS file for later usage (use filename from original TS file + suffix)
        tSS_iag_bed.saveas(output_dir / (TS_fPrefix + '_iag.bed'))
        tss_iag_df_right_format = read_bed_as_df(output_dir / (TS_fPrefix + '_iag.bed'), custom_col_name=col_names_bed_intersect)[
            col_names_bed_intersect_keep].drop_duplicates(keep='first')

        tss_iag_df_right_format.to_csv(
            output_dir / (TS_fPrefix + '_iag.bed'), header=None, index=None, sep='\t')
        print("Reduced from %i to %i TSS/TES (without duplicates)" %
              (len(tSS_TSE_ig_bed), len(tss_iag_df_right_format.index)))
        tSS_iag_bed = BedTool(output_dir / (TS_fPrefix + '_iag.bed'))
    else:
        # Load already filtered TS file
        tSS_iag_bed = BedTool(ts_iag_file)

    # Filter P2P file for genes which TS is supplied in tSS_iag_bed (using 'xor' intersection type)
    print('-------- STEP 3: P2P filtering ------------')
    p2p_for_TSS_iag = get_P2P_for_with_TSS(p2p_bedpe, tSS_iag_bed, tF_bed)

    # Save p2p_for_TSS_iag in order to load with pandas
    p2p_for_TSS_iag.saveas(output_dir / (resFPrefix + int_p2P_name))
    p2p_for_TSS_iag_df = read_bed_as_df(
        output_dir / (resFPrefix + int_p2P_name), custom_col_name=col_names_bedpe)

    # Group same P2P and summarize all gene names of associated TS
    is_per_gene_df = get_is_per_gene_df(p2p_for_TSS_iag_df, min_p2p_score)
    is_per_gene_df.to_csv(output_dir / (resFPrefix + "_ms%i_TS_%s_perGene.txt" %
                          (min_p2p_score, ts_kind)), sep='\t', header=None, index=None)
    p2p_for_TSS_iag_df_g = agg_unique_p2p_plus_filter(
        p2p_for_TSS_iag_df, min_p2p_score)

    p2p_for_TSS_iag_df_g.to_csv(output_dir / (resFPrefix + "_ms%i_TS_%s_P2P.bedpe" %
                                (min_p2p_score, ts_kind)), sep='\t', header=None, index=None)
    print("Reduced from %i to %i P2P" %
          (len(p2p_bedpe), len(p2p_for_TSS_iag_df_g.index)))
    # remove intermediate files
    os.remove(output_dir / (resFPrefix + int_p2P_name))
    # os.remove(output_dir / (TS_fPrefix + '_iag.bed'))
    os.remove(tSS_TSE_df_ig_path)


def main(args):
    start_time = time()
    output_directory = Path(args.output_directory)
    # read in all given files
    # if args.genes_filename is not None and args.genes_filename.split('.')[-1]  == "xlsx":
    #     i_genes_df = read_xlxs(path_folder / args.genes_filename, sheet_num=args.sheet_num-1, read_hidden = False)
    #     i_genes_list = list(i_genes_df.loc[:,args.gene_colname])
    if args.genes_filename is not None:
        with open(args.genes_filename, 'r') as f:
            i_genes_list = f.read().split("\n")
    else:
        i_genes_list = None

    tSS_TSE_df = read_bed_as_df(
        args.tSS_TSE_uf, bed_col_names_g) if args.tSS_TSE_uf is not None else None
    h3k4me3_bed = BedTool(
        args.histone_bed_file) if args.histone_bed_file is not None else None

    tf_path = (args.tF_bed_file.split('.bed')[0] + '_temp.bed')
    tf_df = read_bed_as_df(output_directory / args.tF_bed_file, custom_col_name=col_simp_bed, set_numcols=3)
    tf_df.to_csv(output_directory / tf_path,
                 sep='\t', header=None, index=None)

    tF_bed = BedTool(tf_path) if args.tF_bed_file is not None else None

    p2p_bedpe = BedTool(args.p2p_file)

    # get prefix of TS and bedpe files
    resFPrefix = os.path.basename(args.p2p_file).split('.bedpe')[0].replace(
        'A2A', '') + 'x' + os.path.basename(args.tF_bed_file).split('.bed')[0]
    TS_fPrefix = os.path.basename(args.tSS_TSE_uf).split(
        '.bed')[0] + '_c' + args.kind_TS + '_' + os.path.basename(args.histone_bed_file).split('.bed')[0] if args.tSS_TSE_uf is not None else None

    runAllSteps(i_genes_list, tSS_TSE_df, h3k4me3_bed, tF_bed, p2p_bedpe, args.kind_TS,
                args.min_p2p_score, args.tSS_TSE_f, resFPrefix, TS_fPrefix, output_directory)
    os.remove(output_directory / tf_path)
    print("--- Total execution time is %s seconds ---" % (time() - start_time))


# Argument Parsing
description = "Performs filtering of a P2P file using a TSS/TES file, optional one can filter the TSS/TES file beforehand for interessting and active genes (by supplying a gene list, a H3k4me3 and transcription Factor bed file)."
parser = argparse.ArgumentParser(
    description=description, epilog="NOTE: BedTools has to be installed!")
parser.add_argument('-g', nargs='?', dest='genes_filename',  type=str,
                    help='Filename of the txt file with the important genes', required=False)
# parser.add_argument('-sn', nargs='?', default=1, dest='sheet_num', metavar='--xlxsSheetNum',  type=int,
#                     help="Sheet number in the xlsx file, which contains gene names stored in a collumn called 'GeneID' (default: 1)", required=False)
# parser.add_argument('-g_cn', nargs='?', default="GeneID", dest='gene_colname', metavar='--xlxsGeneColname',  type=str,
#                     help="Name of the column containing the genen names in the xlsx file, which contains gene names stored in a collumn called 'GeneID' (default: 'GeneID')", required=False)

parser.add_argument('-uts', nargs='?', dest='tSS_TSE_uf',
                    type=str, help='Path to unfiltered TSS/TES file (bed file)', required=False)
parser.add_argument('-hf', nargs='?', dest='histone_bed_file',
                    type=str, help='Path to Histone bed file for filtering (bed file)', required=False)
parser.add_argument('-p2p', nargs='?', dest='p2p_file',
                    type=str, help='Path to P2P file (bedpe file)', required=True)
parser.add_argument('-tff', nargs='?', dest='tF_bed_file', type=str,
                    help='Path to Transcription Factor bed file for filtering (bed file)', required=True)
parser.add_argument('-fts', nargs='?', dest='tSS_TSE_f',
                    type=str, help='Path to filtered TSS/TES file (bed file)', required=False)
parser.add_argument('-o', nargs='?', dest='output_directory', default='./',
                    type=str, help='Name of the output directory (default: current directory)', required=False)

parser.add_argument('-mS', nargs='?', default=1, dest='min_p2p_score',
                    type=int, help='Lowest score to keep P2P (default: 1)', required=False)
parser.add_argument('-kindTS', nargs='?', dest='kind_TS', type=str, default='TSS',
                    help="Kind of Transcription site, which should be used (default: 'TSS')", choices=['both', 'TSS', 'TES'], required=False)

args = parser.parse_args()
print("Parameter: \n Lowest score to keep P2P: %i \n TS type: %s \n Output Directory: %s" %
      (args.min_p2p_score, args.kind_TS, args.output_directory))
main(args)


# def read_xlxs(filename, sheet_num, read_hidden= False):

#     # Read Excel file as Pandas DataFrame
#     df = pd.read_excel(filename)

#     # Open an Excel workbook
#     workbook = openpyxl.load_workbook(filename)

#     # Create a `Worksheet` object
#     worksheet = workbook[workbook.sheetnames[sheet_num]]

#     # List of indices corresponding to all hidden rows
#     if not read_hidden:
#         hidden_rows_idx = [
#             row - 2
#             for row, dimension in worksheet.row_dimensions.items()
#             if dimension.hidden
#         ]

#         df.drop(hidden_rows_idx, axis=0, inplace=True)

#         # Reset the index
#         df.reset_index(drop=True, inplace=True)
#     return df
