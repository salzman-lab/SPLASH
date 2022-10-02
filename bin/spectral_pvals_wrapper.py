import numpy as np
import pandas as pd
from tqdm import tqdm,tqdm_notebook, tqdm_pandas
import os
import glob
import pickle
import argparse
import scipy
import scipy.stats
import statsmodels.api as sm
from pathlib import Path 
import itertools

### Simple wrapper script to run compute_spectral_pvals.py on all
###   abundant stratified files, then correct the p-values.


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument( ### folder with abundant_stratified_{}.txt.gz files
        "--in_folder",
        type=str
    )
    parser.add_argument( ### folder containing intermediate p-value files
        "--intermediate_fldr",
        type=str
    )
    parser.add_argument( ### output file path, required
        "--outfile_scores",
        type=str
    )
    parser.add_argument( ### samplesheet file, if samplesheet-based cj (metadata) are to be used
        "--samplesheet",
        type=str,
        default=""
    )
    parser.add_argument( ### columns to output. Options are default, metadata, full
        "--output_verbosity",
        type=str,
        default="default"
    )
    args = parser.parse_args()
    return args


def main():
    args = get_args()

    print('running main')
    for f_seq in itertools.product(['A','T','C','G'],repeat=3):
        f_seq = "".join(f_seq)
        command = """python /oak/stanford/groups/horence/tavorb/stringstats_tavor/stringstats/bin/compute_spectral_pvals.py --outfile_scores {}/int_pvals_{}.tsv --samplesheet {} --infile {}/abundant_stratified_{}.txt.gz --output_verbosity {}""".format(args.intermediate_fldr,f_seq
        ,args.samplesheet,
        args.in_folder,f_seq,
        args.output_verbosity)
        print("running: ", command)
        os.system(command)


    print('aggregating')
    dfs = []
    for df_path in glob.glob(args.intermediate_fldr+"/*.tsv"):
        try:
            dfs.append(pd.read_csv(df_path.strip(), sep='\t'))
        except pd.errors.EmptyDataError:
            pass


    df = pd.concat(dfs)

    ## if it's an empty df, just output the df
    outdf = df.copy()

    if not df.empty:

        ### create corrected pvalues
        newDF = outdf[outdf.columns[~outdf.columns.str.startswith('pval')]]

        for c in outdf.columns[outdf.columns.str.startswith('pval')]:
            _, c_corrected, _, _ = sm.stats.multipletests(outdf[c], method='fdr_by')
            newDF.loc[:,c+'_corrected'] = c_corrected
            
        newDF.sort_values('pval_SVD_corrAnalysis_corrected').to_csv(args.outfile_scores, sep='\t', index=False)





main()