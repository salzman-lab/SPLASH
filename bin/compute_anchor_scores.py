#!/usr/bin/env python

import gzip
import argparse
import pandas as pd
import numpy as np


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--samplesheet",
        type=str,
        help='input samplesheet'
    )
    parser.add_argument(
        "--outfile_counts_distances",
        type=str
    )
    parser.add_argument(
        "--outfile_anchor_scores",
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


def get_distance_df(anchor, df):
    """
    Takes a df of anchor counts and returns a df of hamming distances per target, relative to target abundance
    """
    # exhaustive method if anchor we are at 100K reads or if we have not seen this anchor before
    min_dists = []

    # for each targets, compare hamming distances between itself and all previous targets
    for index, row in df.iterrows():
        
        # current target
        sequence = row['target']

        # get all preceeding targets
        candidates = (
            df['target']
            .loc[0:index-1]
            .to_numpy()
        )            

        # if it's the first and most abundant target, it gets a distance of 0
        if index == 0:
            min_dist = 0
        # else, assign the minimum distance to all preceeding targets
        else:
            min_dist = min([distance(x, sequence) for x in candidates])

        # record target and it's distance
        min_dists.append((sequence, min_dist))
    
    # convert to df, to be used as scalar multiplier of the numerator
    distance_df = (
        pd.DataFrame(
            min_dists, 
            columns=['target', 'min_distance']
        )
        .set_index("target")
    )
    distance_df.index.name = None

    return distance_df


def get_distance_scores(anchor, counts):
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

    # take the top 50 most abundant anchors
    counts = (
        counts
        .sort_values( # sort by counts
            'count', 
            ascending=False
        )
        .drop( # delete it because you don't need it anymore
            'count',
            axis=1
        )
        .reset_index()
        .head(50)
    )

    # now that the targets are sorted by abundance, get the distances per target
    distance_df = get_distance_df(anchor, counts)
    
    # make version of counts with distance column
    counts_distances = counts.copy()
    counts_distances['anchor'] = anchor
    counts_distances = pd.merge(
        counts_distances, 
        distance_df,
        left_on='target',
        right_index=True
    )
    first_column = counts_distances.pop('anchor')
    counts_distances.insert(0, 'anchor', first_column)
    
    # get sum(n_i * d_i) term
    numerator = (
        counts
        .set_index('target')
        .multiply(                  # multiply each count by the minimum distance
            distance_df['min_distance'], 
            axis=0
        )
        .sum(axis=0)                # get sum of (counts * min_dist) across samples
    )

    # get sum(n_i) term
    denominator = (
        counts
        .set_index('target')
        .sum(axis=0)                # get sum of counts across samples
    )

    # return new row
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

    # read in all the dfs
    with open(args.samplesheet) as file:
        df_paths = file.readlines()

    counts = pd.read_csv(df_paths[0].strip(), sep='\t')
    counts=counts.rename(columns = {'adj_kmer':'target'})
    for df_path in df_paths[1:]:
        try:
            df = pd.read_csv(df_path.strip(), sep='\t')
            df=df.rename(columns = {'adj_kmer':'target'})

            counts = counts.merge(
                df,
                on=['anchor', 'target'],
                how='outer'
            )
            del(df)
        except pd.errors.EmptyDataError:
            pass

    counts = counts.fillna(0)

    # get distance scores of anchors that are not in the whitelist (if an anchor has ever been called significant, no need to compute its scores)
    anchor_scores_df = pd.DataFrame()
    counts_distances_df = pd.DataFrame()
    
    # get score for each anchor and append score to output df    
    for anchor, df in counts.groupby('anchor'):
        anchor_scores, counts_distances = get_distance_scores(anchor, df)
        anchor_scores_df = anchor_scores_df.append(anchor_scores)
        counts_distances_df = counts_distances_df.append(counts_distances)

    anchor_scores_df['anchor_score'] = anchor_scores_df.std(axis=1)
    anchor_scores_df.index.name = 'anchor'
    anchor_scores_df.reset_index().to_csv(args.outfile_anchor_scores, sep='\t', index=False)

    counts_distances_df.to_csv(args.outfile_counts_distances, sep='\t', index=False)


main()
