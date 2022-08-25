#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--samplesheet",
        type=str
    )
    parser.add_argument(
        "--is_10X",
        action='store_true'
    )
    parser.add_argument(
        "--run_unsupervised_pvals",
        action='store_true'
    )
    parser.add_argument(
        "--outfile",
        type=str,
        default=''
    )
    args = parser.parse_args()
    return args


def main():
    args = get_args()

    if args.is_10X:

        # read in fastq samplesheet
        samplesheet = pd.read_csv(
            args.samplesheet,
            names=['id', 'metadata']
        )

        if args.run_unsupervised_pvals:
            # create one samplesheet without Cjs, to be assigned during computation
            samplesheet[['id']].to_csv(
                args.outfile,
                index=False,
                header=False
            )

        else:
            for m in samplesheet['metadata'].unique():
                p = samplesheet.copy()
                p['cj'] = np.where(
                    samplesheet['metadata']==m,
                    1,
                    -1
                )
                # output each pairwise samplesheet
                p[['id', 'cj']].to_csv(
                    f'pairwise_samplesheet_{m}.csv',
                    index=False,
                    header=False
                )

    else:

        # read in fastq samplesheet, which contains per-sample metadata
        samplesheet = pd.read_csv(
            args.samplesheet,
            names=['metadata', 'id', 'fastq']
        )

        if args.run_unsupervised_pvals:
            # create one samplesheet with all Cj=1
            samplesheet[['id']].to_csv(
                args.outfile,
                index=False,
                header=False
            )

        else:
            # for each unique metadata value, create a new samplesheet where:
            #   1. metadata value = 1
            #   2. all other values = -1
            for m in samplesheet['metadata'].unique():
                p = samplesheet.copy()
                p['cj'] = np.where(
                    samplesheet['metadata']==m,
                    1,
                    -1
                )
                # output each pairwise samplesheet
                p[['id', 'cj']].to_csv(
                    f'pairwise_samplesheet_{m}.csv',
                    index=False,
                    header=False
                )


main()
