#!/usr/bin/env python3

import pandas as pd
import statsmodels.api as sm
import sys
import argparse
import glob


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--fdr_threshold",
        type=float,
        default = .05
    )
    parser.add_argument(
        "--outfile_scores",
        type=str
    )
    parser.add_argument(
        "--all_anchors_pvals_file",
        action='store_true'
    )

    args = parser.parse_args()
    return args


def main():
    args = get_args()

    print('aggregating')
    dfs = []
    for df_path in glob.glob("scores_*tsv"):

        try:
            dfs.append(pd.read_csv(df_path.strip(), sep='\t'))
        except pd.errors.EmptyDataError:
            pass

    df = pd.concat(dfs)


    outdf = df.copy()

    if not df.empty:
        _, pv_corrected,_, _ = sm.stats.multipletests(df.pv_hash_both, alpha=.05, method='fdr_by')
        outdf['pv_hash_both_corrected'] = pv_corrected

        if args.all_anchors_pvals_file:
            outdf.to_csv(f"all_{args.outfile_scores}", sep='\t', index=False)

        outdf = outdf[outdf.pv_hash_both_corrected < args.fdr_threshold]

    outdf.to_csv(args.outfile_scores, sep='\t', index=False)




main()
