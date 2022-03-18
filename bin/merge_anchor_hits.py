#!/usr/bin/env python

import gzip
import argparse
import pandas as pd
import numpy as np
import os


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--samplesheet",
        type=str
    )
    parser.add_argument(
        "--summary_table",
        type=str
    )
    parser.add_argument(
        "--outfile",
        type=str
    )    
    args = parser.parse_args()
    return args


def main():
    args = get_args()

    # read in all the dfs
    with open(args.samplesheet) as file:
        hits_paths = file.readlines()

    # read in summary table
    df = pd.read_csv(args.summary_table, sep='\t')

    # iteratively marge in hits on anchor column
    for hits_path in hits_paths:
        hits = pd.read_csv(hits_path.strip(), sep='\t')
        df = df.merge(hits, on='anchor')
    
    df.to_csv(
        args.outfile,
        sep='\t',
        index=False
    )
    

main()