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
        "--samplesheet",
        type=str
    )
    parser.add_argument(
        "--outfile_scores",
        type=str
    )
    parser.add_argument(
        "--outfile_all_anchors_pvals",
        type=str
    )
    parser.add_argument(
        "--outfile_Cjs",
        type=str
    )

    args = parser.parse_args()
    return args


def main():
    args = get_args()

    # Determine if we are using sheet Cjs
    with open(args.samplesheet, 'r') as f:
        cols = f.readline().split(',')

    if len(cols) == 1:
        useSheetCjs = False
    elif len(cols) == 2:
        useSheetCjs = True

    print('aggregating')
    dfs = []
    for df_path in glob.glob("scores_*tsv"):
        try:
            dfs.append(pd.read_csv(df_path.strip(), sep='\t'))
        except pd.errors.EmptyDataError:
            pass

    df = pd.concat(dfs)

    ## if it's an empty df, just output the df
    out_pvals = df.copy()

    if not df.empty:

        # Correct random Cjs
        _, pval_random_corrected, _, _ = sm.stats.multipletests(df['pval_random'], alpha=.05, method='fdr_by')
        out_pvals['pval_random_corrected'] = pval_random_corrected

        # If using sheet Cjs, also correct samplesheet and aggregated Cjs
        if useSheetCjs:
            _, pval_samplesheet_corrected, _, _ = sm.stats.multipletests(df['pval_samplesheet'], alpha=.05, method='fdr_by')
            out_pvals['pval_samplesheet_corrected'] = pval_samplesheet_corrected

            _, pval_aggregated_corrected, _, _ = sm.stats.multipletests(df['pval_aggregated'], alpha=.05, method='fdr_by')
            out_pvals['pval_aggregated_corrected'] = pval_aggregated_corrected

        # Output all anchors and pvals
        out_pvals.to_csv(args.outfile_all_anchors_pvals, sep='\t', index=False)

        # Reject
        if useSheetCjs:
            out_pvals = out_pvals[out_pvals['pval_aggregated_corrected'] < args.fdr_threshold]
        else:
            out_pvals = out_pvals[out_pvals['pval_random_corrected'] < args.fdr_threshold]

        # Drop optHash
        out_pvals.drop(columns=['optHash'], axis=1, inplace=True)

        # Extract columns for cjs and stats files
        cj_cols = [c for c in out_pvals.columns if c == "anchor" or "cj_rand_opt_" in c]
        stats_cols = [c for c in out_pvals.columns if c == "anchor" or "cj_rand_opt_" not in c]

        # Subset and output
        out_cjs = out_pvals[cj_cols]
        cj_cols_mod = [c.replace("cj_rand_opt_", "") for c in out_cjs.columns]
        out_cjs.columns = cj_cols_mod
        out_cjs.to_csv(args.outfile_Cjs, sep='\t', index=False)

        # Subset
        out_pvals = out_pvals[stats_cols]

    # Output
    out_pvals.to_csv(args.outfile_scores, sep='\t', index=False)


main()
