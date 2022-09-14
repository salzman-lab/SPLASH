#!/usr/bin/env python3

import argparse
import pandas as pd
from Bio.Seq import Seq
import numpy as np
from Bio.Blast.Applications import NcbiblastnCommandline
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
pd.options.mode.chained_assignment = None  # default='warn'


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--anchors_pvals",
        type=str
    )
    # parser.add_argument(
    #     "--genome_annotations",
    #     type=str
    # )
    # parser.add_argument(
    #     "--element_annotations",
    #     type=str
    # )
    parser.add_argument(
        "--outfile",
        type=str
    )
    args = parser.parse_args()
    return args


def get_ann(row, seq_type):
    priority_list = [
        f'{seq_type}_hits_UniVec',
        f'{seq_type}_hits_illumina_adapters',
        f'{seq_type}_hits_escherichia_phage_phiX174',
        f'{seq_type}_hits_RF_all',
        f'{seq_type}_hits_dfam_te_eukaryota',
        f'{seq_type}_hits_tncentral_te_prokaryotes_final',
        f'{seq_type}_hits_mge_aclame_genes_all_0',
        f'{seq_type}_hits_ice_iceberg',
        f'{seq_type}_hits_direct_repeats',
        f'{seq_type}_hits_spacers',
        f'{seq_type}_hits_eukaryota_its1_itstonedb',
        f'{seq_type}_hits_eukaryota_its2_biozentrum',
        f'{seq_type}_hits_WBcel235',
        f'{seq_type}_hits_TAIR10',
        f'{seq_type}_hits_grch38_1kgmaj'
    ]

    for ref in priority_list:
        if row.loc[ref] != '*':
            row.loc[f"{seq_type}_top_ann"] = ref
            row.loc[f"{seq_type}_top_ann_hit"] = row.loc[ref]
            row.loc[f"{seq_type}_top_ann_hit_pos"] = row.loc[ref.replace("hits", "hits_pos")]
            break
    return row


def add_summary(df, ann_table, seq_type):
    c_list = [c for c in ann_table.columns if 'hits_pos' not in c]

    ann_table['all_unannotated'] = (
        ann_table[c_list]
        .drop(seq_type, axis=1)
        .replace('*', np.NaN)
        .isnull()
        .all(axis=1)
    )

    # df with bowtie2 ann
    ann_df = ann_table[ann_table['all_unannotated'] == False]

    ann_df[f"{seq_type}_top_ann"] = None
    ann_df[f"{seq_type}_top_ann_hit"] = None
    ann_df[f"{seq_type}_top_ann_hit_pos"] = None

    ann_df = ann_df.apply(
        lambda row:
        get_ann(row, seq_type), axis=1
    )

    ann_df[f"{seq_type}_num_ann"] = (
            ann_df[c_list]
            .drop(seq_type, axis=1)
            .replace('*', np.nan)
            .notnull()
            .astype(int)
            .sum(axis=1)
        )

    ann_df[f"{seq_type}_annotation_source"] = "bowtie2"
    bowtie_df = ann_df[[seq_type, f"{seq_type}_top_ann", f"{seq_type}_top_ann_hit", f"{seq_type}_top_ann_hit_pos", f"{seq_type}_num_ann", f"{seq_type}_annotation_source"]]

    return bowtie_df


def main():
    args = get_args()

    ## read in anchors and their anchors_pvals
    anchors_pvals = (
        pd.read_csv("anchors_pvals.tsv", sep='\t')
        .drop_duplicates()
    )

    ## account for decoy anchors_pvals
    if len(anchors_pvals.columns) == 2:
        anchors_pvals.columns = ['anchor', 'control_pvalue']

    anchors_pvals_cols = [c for c in anchors_pvals.columns if "cj_" not in c]
    anchors_pvals = anchors_pvals[anchors_pvals_cols]

    genome_annotations = pd.read_csv("genome_annotations_anchors.tsv", sep='\t')

    df = pd.merge(anchors_pvals, genome_annotations, on='anchor')

    ann_anchors = pd.read_csv("element_annotations_anchors.tsv", sep='\t')

    summarized_anchors = add_summary(df, ann_anchors, 'anchor')

    df = pd.merge(df, summarized_anchors, on='anchor', how='left')

    df.to_csv(args.outfile, sep='\t', index=False)


main()
