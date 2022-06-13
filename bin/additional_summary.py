#!/usr/bin/env python3

import argparse
import pandas as pd

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--summary",
        type=str
    )
    parser.add_argument(
        "--consensus_called_exons",
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

    try:
        info = pd.read_csv(args.consensus_called_exons, sep='\t').drop_duplicates()
        info = info[['anchor', 'anchor_local_gene', 'anchor_end_to_end_gene', 'consensus_gene']]
        info = (
            info
            .groupby(['anchor','anchor_local_gene', 'anchor_end_to_end_gene'])
            .agg(pd.Series.mode)
            .reset_index()
        )
        info.columns = ['anchor', 'anchor_local_gene', 'anchor_end_to_end_gene', 'consensus_gene_mode']
    except:
        info = pd.DataFrame(columns = ['anchor', 'anchor_local_gene', 'anchor_end_to_end_gene', 'consensus_gene_mode'])

    summary = pd.read_csv(args.summary, sep='\t')

    outdf = pd.merge(summary, info, on='anchor', how='left')

    outdf.to_csv(args.outfile, sep='\t', index=False)


main()
