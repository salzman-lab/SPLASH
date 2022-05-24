#!/usr/bin/env python3

import pandas as pd
import statsmodels.api as sm
import sys
import argparse

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--infile",
        type=str
    )
    parser.add_argument(
        "--pval_threshold",
        type=float
    )
    parser.add_argument(
        "--outfile_anchors",
        type=str
    )
    parser.add_argument(
        "--outfile_scores",
        type=str
    )
    args = parser.parse_args()
    return args


def main():
    args = get_args()

    with open(args.infile) as file:
        df_paths = file.readlines()

    dfs = []
    for df_path in df_paths:

        try:
            dfs.append(
                pd.read_csv(df_path.strip(), sep='\t')
            )

        except pd.errors.EmptyDataError:
            pass

    df = pd.concat(dfs)

    outdf = pd.DataFrame(columns=['anchor', 'pv_Rand'])

    if not df.empty:
        reject, pvals_corrected,_, _ = sm.stats.multipletests(df.pv_Rand, alpha=.05, method='fdr_by')

        outdf = df[reject]

        outdf = (
            outdf[outdf['pv_Rand'] < args.pval_threshold]
            .sort_values('pv_Rand')
            .head(5000)
        )

    outdf.to_csv(args.outfile_scores, sep='\t', index=False, header=False)
    outdf[['anchor']].to_csv(args.outfile_anchors, sep='\t', header=False, index=False)


main()

