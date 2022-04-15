#!/usr/bin/env python3

import os
import argparse
import pandas as pd
import numpy as np

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--use_std",
        type=str
    )
    parser.add_argument(
        "--anchor_sample_scores",
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
    anchor_sample_scores = pd.read_csv(args.anchor_sample_scores, sep='\t')

    # read in samples from samplesheet
    # [fastq_file, optional group_id]
    sample_list = pd.read_csv(
        args.samplesheet,
        header=None,
        sep='\t'
    )

    # redefine bool
    if args.use_std == "true":
        use_std = True
    elif args.use_std == "false":
        use_std = False
    if sample_list.shape[1] == 1:
        use_std = True

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
        .drop('distance', axis=1)       # df of anchor-targets x sample counts
        .groupby('anchor')              # over all anchors
        .sum()                          # get the per-sample sums
    )

    sqrt_n_j = (
        n_j                             # per-sample sums
        .apply(np.sqrt)                 # square-root the per anchor-sample sums
    )

    anchor_sample_scores = (
        anchor_sample_scores                    # file of anchor_sample_scores per anchor and per anchor-sample
        .set_index('anchor')                    # df of anchor x S_j
        .drop('anchor_sample_std', axis=1)      # drop anchor_score column
    )

    if use_std:
        norm_summary_score = (
            anchor_sample_scores                # df of anchor x S_j
            .mul(pd.Series(group_ids_dict))     # multiply by group_ids
            .mul(sqrt_n_j, axis=1)              # multiply by n_j (square-root of per anchor-sample sums)
            .std(axis=1)                        # std over samples, to get a per-anchor score
        )
    else:
        norm_summary_score = (
            anchor_sample_scores                # df of anchor x S_j
            .mul(pd.Series(group_ids_dict))     # multiply by group_ids
            .mul(sqrt_n_j, axis=1)              # multiply by n_j (square-root of per anchor-sample sums)
            .sum(axis=1)                        # sum over samples, to get a per-anchor score
        )

    sum_sqrt_n_j_c_j = (
        sqrt_n_j                                # square-root of per anchor-sample sums
        .mul(pd.Series(group_ids_dict))         # multiply by group_ids
        .sum(axis=1)                            # sum over samples, to get a per-anchor score
    )

    # intialise the scores table
    scores_table = n_j.copy()
    scores_table['total_anchor_counts'] = scores_table.sum(axis=1)
    scores_table = (
        pd.merge(
            scores_table,
            pd.DataFrame(norm_summary_score, columns=['norm_summary_score']),
            on='anchor'
        )
    )

    # initialise columns
    scores_table['first_moment'] = None
    scores_table['second_moment'] = None
    scores_table['third_moment'] = None
    scores_table['fourth_moment'] = None
    scores_table['fourth_central_moment'] = None
    scores_table['expectation'] = None
    scores_table['variance_distance'] = None

    # for each anchor, get expectation and variance_distance
    for anchor, df in counts.groupby('anchor'):

        # subset and rename distance as i
        distances = (
            df
            .drop(['anchor', 'target'], axis=1)
            .rename(columns={'distance' : 'i'})
        )

        # initialise
        p = pd.DataFrame()
        p['i'] = range(0, args.kmer_size+1)
        p['counts'] = 0

        for i, df in distances.groupby('i'):
            # get total counts for targets with that distance
            count = (
                df
                .drop('i', axis=1)
                .values
                .sum()
            )
            # set the count
            p.loc[p['i']==i, 'counts'] = count

        # get fraction of times each distance occurs as p_hat
        p['p_hat'] = p['counts']/p['counts'].sum()

        # define moments
        first_moment = (p['p_hat'] * p['i']).sum()
        second_moment = (p['p_hat'] * p['i']**2).sum()
        third_moment = (p['p_hat'] * p['i']**3).sum()
        fourth_moment = (p['p_hat'] * p['i']**4).sum()

        mu = (p['i'] * p['p_hat']).sum()
        fourth_central_moment = (p['i'] * (p['i'] - mu)**4).sum()

        # define variance_distance
        variance_distance = (p['p_hat'] * p['i']**2).sum() - ((p['p_hat'] * p['i']).sum()) ** 2

        # add these values to the score table
        scores_table.loc[anchor, 'variance_distance'] = variance_distance
        scores_table.loc[anchor, 'first_moment'] = first_moment
        scores_table.loc[anchor, 'second_moment'] = second_moment
        scores_table.loc[anchor, 'third_moment'] = third_moment
        scores_table.loc[anchor, 'fourth_moment'] = fourth_moment
        scores_table.loc[anchor, 'fourth_central_moment'] = fourth_central_moment

    # output scores table
    (
        scores_table
        .reset_index()
        .to_csv(args.outfile_norm_scores, sep='\t', index=False)
    )


main()
