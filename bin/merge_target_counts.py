#!/usr/bin/env python3

import gzip
import argparse
import pandas as pd
import warnings
import nltk
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
    # merge all target_counts files
    counts = pd.read_csv(df_paths[0].strip(), sep='\t')
    for df_path in df_paths[1:]:
        try:
            df = pd.read_csv(df_path.strip(), sep='\t')
            counts = counts.merge(
                df,
                on=['anchor', 'target'],
                how='outer'
            )
            del(df)
        except pd.errors.EmptyDataError:
            pass

    # fill NA with 0
    counts = counts.fillna(0)

    # output anchor targets counts file
    counts.to_csv(args.outfile_counts_distances, sep='\t', index=False)


main()
