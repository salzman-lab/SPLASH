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
        "--base_dir",
        type=str
    )
    parser.add_argument(
        "--ann_file",
        type=str,
        default =""
    )
        
    args = parser.parse_args()
    return args


def main():
    args = get_args()

    print('aggregating')
    dfs = []
    for df_path in glob.glob(args.base_dir+"/pvals_*.csv"):

        try:
            tmpdf = pd.read_csv(df_path.strip(), sep='\t')
            if 'pvals_all_ann.csv' in df_path.strip():#### skip over pvals_all_ann.csv
                print('found pvals_all_ann.csv file,overwriting')
                continue 
            dfs.append(tmpdf)
        except pd.errors.EmptyDataError:
            pass

    df = pd.concat(dfs)
#     print(df.columns)
    ### read in annotations and merge
    if args.ann_file != "":
        print("using annotation file")
        ann_genome = pd.read_csv(args.ann_file,sep='\t')
        df = df.merge(ann_genome)

    outdf = df.copy()

    if not df.empty:
        _, pv_hash_corrected,_, _ = sm.stats.multipletests(df.pv_hash, alpha=.05, method='fdr_by')
        outdf['pv_hash_corrected'] = pv_hash_corrected

        outdf = outdf[outdf.pv_hash_corrected < args.fdr_threshold]
        
        
    print('writing')
    outdf.to_csv(args.base_dir + '/pvals_all_ann.csv', sep='\t', index=False)


main()