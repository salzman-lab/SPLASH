#!/usr/bin/env python3

import gzip
import argparse
import pandas as pd
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--targets_samplesheet",
        type=str,
        help='input targets_samplesheet'
    )
    parser.add_argument(
        "--bound_distance",
        type=str
    )
    parser.add_argument(
        "--max_distance",
        type=int
    )
    parser.add_argument(
        "--kmer_size",
        type=int
    )
    parser.add_argument(
        "--outfile_counts_distances",
        type=str
    )
    parser.add_argument(
        "--outfile_anchor_sample_scores",
        type=str
    )
    args = parser.parse_args()
    return args


def distance(seq1, seq2):
    """
    Returns the Hamming distance between 2 strings
    """
    distance = 0
    for x,y in zip(seq1, seq2):
        distance = distance + (x != y)

    return distance


def get_distance_df(anchor, df, bound_distance):
    """
    Takes a df of anchor counts and returns a df of hamming distances per target, relative to target abundance
    """

    # initalise
    min_dists = []

    # in this dataframe of targets sorted by decreasing abundance,
    # compare hamming distances between each target and all preceding targets
    for index, row in df.iterrows():

        # current target
        sequence = row['target']

        # get all preceding targets
        candidates = (
            df['target']
            .loc[0:index-1]
            .to_numpy()
        )

        # if it's the first and most abundant target, it is initialized with a distance of 0
        if index == 0:
            min_dist = 0

        # else, assign min_dist to be the min dist between itself and all preceding targets
        else:
            min_dist = min([distance(x, sequence) for x in candidates])

        # if we are bounding the distance, use the min of the max_distance and the current min_dist
        if bound_distance:
            min_dist = min(bound_distance, min_dist)

        # record target and its distance
        min_dists.append((sequence, min_dist))

    # convert to df, to be used as scalar multiplier of the numerator
    distance_df = (
        pd.DataFrame(
            min_dists,
            columns=['target', 'distance']
        )
        .set_index('target')
    )
    distance_df.index.name = None

    return distance_df


def get_distance_scores(anchor, counts, bound_distance, kmer_size):
    """
    Get distance scores for one anchor
    """
    # create df of targets x sample_counts
    counts = (
        counts
        .drop(
            ['anchor'],
            axis=1
        )
        .set_index('target')
    )

    # make count column to sort anchors by abundance
    counts['count'] = counts.sum(axis=1)

    # take the top 50 most abundant targets
    counts = (
        counts
        .sort_values(       # sort by counts
            'count',
            ascending=False
        )
        .drop(              # delete it because you don't need it anymore
            'count',
            axis=1
        )
        .reset_index()      # reset the index
        .head(50)           # take top 50 most abundant targets
    )

    # now that the targets are sorted by abundance, get the distances per target
    distance_df = get_distance_df(anchor, counts, bound_distance)

    # make df of counts and distances
    counts_distances = counts.copy()

    # set anchor column
    counts_distances['anchor'] = anchor
    # merge anchor-target counts with target distnaces
    counts_distances = pd.merge(
        counts_distances,
        distance_df,
        left_on='target',
        right_index=True
    )
    # move anchor column to the first column
    first_column = counts_distances.pop('anchor')
    counts_distances.insert(0, 'anchor', first_column)

    # make anchor counts df for getting mu
    anchor_counts = counts_distances.drop('anchor', axis=1)

    # intialise
    p = pd.DataFrame()
    p['i'] = range(0, kmer_size+1)
    p['counts'] = 0

    # for each distance, get the the total number of occurrences of that distance
    for i, df in counts_distances.groupby('distance'):
        count = (
            df
            .drop(['anchor', 'target', 'distance'], axis=1)
            .values
            .sum()
        )
        p.loc[p['i']==i, 'counts'] = count

    # get proportion of counts with that distance
    p['p_hat'] = p['counts'] / p['counts'].sum()

    # get mu
    mu = (p['i'] * p['p_hat']).sum() / kmer_size

    # get sum(n_i * d_i) term
    numerator = (
        counts
        .set_index('target')                # make df of targets x sample_counts
        .multiply(                          # multiply each target count by its respective min_dist
            distance_df['distance'] - mu,
            axis=0
        )
        .sum(axis=0)                        # get sum of (counts * min_dist) across samples
    )

    # get sum(n_i) term
    denominator = (
        counts
        .set_index('target')                # make df of targets x sample_counts
        .sum(axis=0)                        # get sum of counts across samples
    )

    # get S_i = sum(n_i * d_i) / sum(n_i)
    row = (
        numerator
        .divide(
            denominator,
            fill_value=0,
            axis=0
        )
        .rename(anchor)
    )
    return row, counts_distances


def main():
    args = get_args()

    # define
    if args.bound_distance == 'true':
        bound_distance = args.max_distance
    elif args.bound_distance == 'false':
        bound_distance = None

    # read in target_counts paths
    with open(args.targets_samplesheet) as file:
        df_paths = file.readlines()
    # merge all target_counts files
    counts = pd.read_csv(df_paths[0].strip(), sep='\t')
    for df_path in df_paths[1:]:
        try:
            df = pd.read_csv(df_path.strip(), sep='\t')
            counts = counts.merge(
                df,
                on=['anchor', 'target'],
                how='outer'
            )
            del(df)
        except pd.errors.EmptyDataError:
            pass

    # fill NA with 0
    counts = counts.fillna(0)

    counts.to_csv('counts.tsv', sep='\t', index=False)

    # intialise
    anchor_sample_scores_df = pd.DataFrame()
    counts_distances_df = pd.DataFrame()

    # get scores for each anchor
    for anchor, df in counts.groupby('anchor'):
        anchor_sample_scores, counts_distances = get_distance_scores(
            anchor,
            df,
            bound_distance,
            args.kmer_size
        )
        # append anchor score row
        anchor_sample_scores_df = anchor_sample_scores_df.append(anchor_sample_scores)
        # append coutns_distances df
        counts_distances_df = counts_distances_df.append(counts_distances)

    # get anchor_score value as std of all anchor sample scores
    anchor_sample_scores_df['anchor_sample_std'] = anchor_sample_scores_df.std(axis=1)
    anchor_sample_scores_df.index.name = 'anchor'
    # output anchor scores file
    anchor_sample_scores_df.reset_index().to_csv(args.outfile_anchor_sample_scores, sep='\t', index=False)

    # output anchor targets counts file
    counts_distances_df.to_csv(args.outfile_counts_distances, sep='\t', index=False)


main()
