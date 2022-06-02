#!/usr/bin/env python3

import pandas as pd
import statsmodels.api as sm
import sys
import argparse
import glob

### writes out file:
### args.base_dir + '/pvals_all_ann.csv'

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--fdr_threshold",
        type=float
    )
    parser.add_argument(
        "--base_dir",
        type=str
    )
    # parser.add_argument(
    #     "--annFile",
    #     type=str
    # )
    parser.add_argument(
        "--outfile_scores",
        type=str
    )
    args = parser.parse_args()
    return args


def main():
    args = get_args()

    dfs = []
    for df_path in glob.glob("scores*tsv"):
        try:
            dfs.append(
                pd.read_csv(df_path.strip(), sep='\t')
            )

        except pd.errors.EmptyDataError:
            pass

    df = pd.concat(dfs)

    # ### read in annotations and merge
    # ann_genome = pd.read_csv(args.annFile,sep='\t')
    # df = df.merge(ann_genome)

    outdf = df.copy()

    if not df.empty:
        _, pv_hash_corrected,_, _ = sm.stats.multipletests(df.pv_hash, alpha=.05, method='fdr_by')
        outdf['pv_hash_corrected'] = pv_hash_corrected

        outdf = outdf[outdf.pv_hash_corrected < args.fdr_threshold]

    outdf.to_csv(args.outfile_scores, sep='\t', index=False)


main()
