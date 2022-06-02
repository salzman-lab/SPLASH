#!/usr/bin/env python3

import argparse
import pandas as pd
import glob

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--outfile_anchor_fasta",
        type=str
    )
    parser.add_argument(
        "--outfile_target_fasta",
        type=str
    )
    parser.add_argument(
        "--outfile_anchor_target_counts",
        type=str
    )
    args = parser.parse_args()
    return args


def main():
    args = get_args()

    # read in target_counts paths
    df_paths = glob.glob("*_target_counts.tsv")

    # read in all target_counts files
    dfs = []
    for df_path in df_paths:
        try:
            dfs.append(pd.read_csv(df_path.strip(), sep='\t'))
        except pd.errors.EmptyDataError:
            pass

    # merge all target_counts files, pivot on samples, and replace all NA with 0
    df = (
        pd.concat(dfs)
        .pivot(index=['anchor', 'target'], columns='sample', values='count')
        .reset_index()
        .fillna(0
    ))

    # output anchor targets counts file
    df.to_csv(args.outfile_anchor_target_counts, sep='\t', index=False)

    # output anchor and target fastas
    anchors = (
        df['anchor']
        .drop_duplicates()
        .tolist()
    )

    targets = (
        df['target']
        .drop_duplicates()
        .tolist()
    )

    with open(args.outfile_anchor_fasta, "w") as anchor_fasta:
        for seq in anchors:
            anchor_fasta.write(f'>{seq}\n{seq}\n')
    anchor_fasta.close()

    with open(args.outfile_target_fasta, 'w') as target_fasta:
        for seq in targets:
            target_fasta.write(f'>{seq}\n{seq}\n')
    target_fasta.close()


main()
