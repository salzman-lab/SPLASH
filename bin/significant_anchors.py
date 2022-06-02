#!/usr/bin/env python3

import pandas as pd
import statsmodels.api as sm
import sys
import argparse
import glob

### inputs:
##### outputs of compute_pvals.py: args.base_dir+"/pvals_stratified/pvals_{}.csv"
##### genome annotation file: args.annFile (e.g. /oak/stanford/groups/horence/kaitlin/results_nomad/bulk_RNAseq/paper/AA_antibody_secreting_cells/genome_annotations/genome_annotations_anchor.tsv)
### writes out file:
##### args.base_dir + '/pvals_all_ann.csv'

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
        _, pv_hash_corrected,_, _ = sm.stats.multipletests(df.pv_hash, alpha=.05, method='fdr_by')
        outdf['pv_hash_corrected'] = pv_hash_corrected

        outdf = outdf[outdf.pv_hash_corrected < args.fdr_threshold]


    print('writing')
    outdf.to_csv(args.outfile_scores, sep='\t', index=False)


main()
