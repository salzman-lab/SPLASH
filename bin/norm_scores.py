#!/usr/bin/env python3

import os
import argparse
import pandas as pd
import numpy as np

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--anchor_scores",
        type=str
    )
    parser.add_argument(
        "--anchor_target_counts",
        type=str
    )
    parser.add_argument(
        "--samplesheet",
        type=str
    )
    parser.add_argument(
        "--kmer_size",
        type=int
    )
    parser.add_argument(
        "--outfile_norm_scores",
        type=str
    )
    args = parser.parse_args()
    return args


def main():
    args = get_args()

    # read in files
    counts = pd.read_csv(args.anchor_target_counts, sep='\t')
    scores = pd.read_csv(args.anchor_scores, sep='\t')

    # read in samples from samplesheet
    # [fastq_file, optional group_id]
    sample_list = pd.read_csv(
        args.samplesheet,
        header=None
    )

    # get list of samples from fastq_files
    # if fastq_file = "file1.fastq.gz", sample = "file1"
    samples = (
        sample_list
        .iloc[:,0]
        .apply(
            lambda x:
            os.path.basename(x).split('.')[0]
        )
        .tolist()
    )

    # make group_ids dict, default to 1 if no group_ids provided
    group_ids_dict = {}
    if sample_list.shape[1] == 1:
        for i in range(0, len(samples)):
            group_ids_dict[samples[i]] = 1
    else:
        group_ids = sample_list.iloc[:,1].tolist()
        for i in range(0, len(samples)):
            group_ids_dict[samples[i]] = group_ids[i]

    # define terms
    n_j = (
        counts                          # file of anchor-target-sample counts
        .drop('min_distance', axis=1)   # df of anchor-targets x sample counts
        .groupby('anchor')              # over all anchors
        .sum()                          # get the per-sample sums
    )

    sqrt_n_j = (
        n_j                             # per-sample sums
        .apply(np.sqrt)                 # square-root the per anchor-sample sums
    )

    scores = (
        scores                          # file of scores per anchor and per anchor-sample
        .set_index('anchor')            # df of anchor x S_j
        .drop('anchor_score', axis=1)   # drop anchor_score column
    )

    norm_summary_score = (
        scores                          # df of anchor x S_j
        .assign(**group_ids_dict)       # map group_ids to samples
        .mul(scores)                    # multiply S_j by group_ids
        .mul(sqrt_n_j, axis=1)          # multiply by n_j (square-root of per anchor-sample sums)
        .sum(axis=1)                    # sum over samples, to get a per-anchor score
    )

    sum_sqrt_n_j_c_j = (
        sqrt_n_j                        # square-root of per anchor-sample sums
        .assign(**group_ids_dict)       # map sums to group_ids
        .mul(sqrt_n_j)                  # multiply n_j by group_ids
        .sum(axis=1)                    # sum over samples, to get a per-anchor score
    )

    # intialise the scores table
    scores_table = (
        pd.merge(
            n_j.copy(),
            pd.DataFrame(norm_summary_score, columns=['norm_summary_score']),
            on='anchor'
        )
    )

    # initialise
    scores_table['expectation'] = None
    scores_table['variance_distance'] = None

    # for each anchor, get expectation and variance_distance
    for anchor, df in counts.groupby('anchor'):

        # subset and rename min_distance as i
        df = (
            df[['min_distance']]
            .rename(columns={'min_distance' : 'i'})
        )

        # get fraction of times each distance occurs as p_hat
        dist_dict = (
            df['i']
            .value_counts(normalize=True)
            .to_dict()
        )

        # map p_hat to i
        df['p_hat'] = df['i'].map(dist_dict)

        # get i (targets are already sorted by abundance)
        df['i'] = df.index

        # define X
        expectation = sum((df['i'] * df['p_hat'] / args.kmer_size) * sum_sqrt_n_j_c_j[anchor])

        # define variance_d
        variance_distance_first_sum = sum(df['p_hat'] * df['i']**2 / args.kmer_size)
        variance_distance_second_sum = (sum(df['p_hat'] * df['i'] / args.kmer_size)) ** 2
        variance_distance = variance_distance_first_sum - variance_distance_second_sum

        # add these values to the score table
        scores_table.loc[anchor, 'expectation'] = expectation
        scores_table.loc[anchor, 'variance_distance'] = variance_distance

    # output scores table
    (
        scores_table
        .reset_index()
        .to_csv(args.outfile_norm_scores, sep='\t', index=False)
    )


main()
