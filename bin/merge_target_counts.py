#!/usr/bin/env python3

import gzip
import argparse
import pandas as pd
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import utils

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--targets_samplesheet",
        type=str,
        help='input targets_samplesheet'
    )
    parser.add_argument(
        "--outfile_counts_distances",
        type=str
    )
    args = parser.parse_args()
    return args

def main():
    args = get_args()

    # read in target_counts paths
    with open(args.targets_samplesheet) as file:
        df_paths = file.readlines()

    # append dfs to list, setting indices to merge on later
    dfs = []
    for df_path in df_paths:
        dfs.append(
            pd.read_csv(df_path.strip(), sep='\t')
            .set_index(['anchor', 'target'])
        )

    # perform outer merges, filling NA with 0
    counts = (
        dfs[0]
        .join(dfs[1:])
        .fillna(0)
        .reset_index()
    )

    # output anchor targets counts file
    counts.to_csv(args.outfile_counts_distances, sep='\t', index=False)


main()
