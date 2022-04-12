#!/usr/bin/env python3

import pandas as pd
import gzip
import os
import re
import sys
import argparse
import logging
import numpy as np
from Bio import SeqIO
import math
import nltk
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
from datetime import datetime
import Reads


def get_distance(seq1, seq2, distance_type):
    """
    Returns the distance between 2 strings
    """
    if distance_type.lower() == 'hamming':
        distance = 0
        for x,y in zip(seq1, seq2):
            distance = distance + (x != y)
    elif distance_type.lower() == 'jaccard':
        distance = 1 - (nltk.jaccard_distance(set(seq1), set(seq2)))

    return distance


def is_ignorelisted(self, anchor):
    """
    Returns true if an anchor has previously ignorelisted
    """
    if anchor in self.ignorelist:
        return True
    else:
        return False


def compute_phase_1_scores(anchor_counts, group_ids_dict, kmer_size, distance_type, score_type):
    # intialise a df for targets x distances
    distance_df = pd.DataFrame(index=anchor_counts['target'], columns=['distance'])

    # iterate over targets sorted by decreasing abundance
    # get min hamming distance of each target to its preceeding targets
    for index, row in anchor_counts.iterrows():
        target = row['target']
        candidates = (
            anchor_counts['target']
            .loc[0:index-1]
            .to_numpy()
        )

        # if it's the first aka most abundant target, initialise the distance to be 0
        if index == 0:
            min_dist = 0

        # else, get the min distance to all preceeding targets
        else:
            min_dist = min([get_distance(target, x, distance_type) for x in candidates])

        # update target distances
        distance_df.loc[distance_df.index==target, 'distance'] = min_dist
        anchor_counts.loc[anchor_counts['target']==target, 'distance'] = min_dist

    # intialise
    p = pd.DataFrame()
    p['i'] = range(0, kmer_size+1)
    p['counts'] = 0

    # for each distance, get the the total number of occurrences of that distance
    for i, df in anchor_counts.groupby('distance'):
        count = (
            df
            .drop(['target', 'distance'], axis=1)
            .values
            .sum()
        )
        p.loc[p['i']==i, 'counts'] = count

    # get proportion of counts with that distance
    p['p_hat'] = p['counts'] / p['counts'].sum()

    # get u
    mu = (p['i'] * p['p_hat']).sum()

    if score_type == 'slow':
        # center the distances with u
        anchor_counts['distance_centered'] = anchor_counts['distance'] - mu

        # prepare df to get terms
        anchor_counts = (
            anchor_counts
            .drop(['distance'], axis=1)
            .set_index('target')
            )

        m = anchor_counts.drop('distance_centered', axis=1)
        distance_centered = anchor_counts['distance_centered']

        numerator = (
            m
            .mul(pd.Series(group_ids_dict))
            .mul(distance_centered, axis=0)
            .sum(axis=0)
        )

        denominator = m.sum(axis=0)

        phase_1_scores = (
            numerator
            .divide(denominator, fill_value=0, axis=0)
        )

        distances = distance_df['distance']

    elif score_type == "fast":
        phase_1_scores = pd.DataFrame()
        distances = pd.DataFrame()

    return phase_1_scores, distances, mu


def compute_phase_2_scores(previous_score, n, distance, mu, c):
    """
    Returns phase_2 scores, given a distance of a new target
    """
    new_mu = mu * (n/(n+1)) + (distance / (n+1))
    new_score = previous_score * (n/(n+1)) + ((distance - new_mu) * c / (n + 1))

    return new_score


def get_read_chunk(iteration, samples, n_iterations):
    """
    Returns reads from all files for this iteration
    """
    read_chunk = []
    for sample in samples:

        try:
            n_digits = len(str(n_iterations))
            fastq_file = f'{sample}_{str(iteration).zfill(n_digits)}.fastq.gz'

            with gzip.open(fastq_file, 'rt') as fastq:
                for seqreacord in SeqIO.parse(fastq, 'fastq'):
                    read = Reads.Read(str(seqreacord.seq))
                    read_chunk.append((read, sample))
        except:
            print(f'{sample}, {iteration}')

    return read_chunk


def is_valid_anchor_target(anchor, read_counter_freeze, anchor_counts, anchor_freeze_threshold, status_checker):
    """
    Return True if an anchor is valid for further computation
    """
    # set anchor_counter_freeze if we have reached the anchor_freeze_threshold
    if len(anchor_counts.total_counts) >= anchor_freeze_threshold:
        anchor_counter_freeze = True
    else:
        anchor_counter_freeze = False

    # break condition if not at read_threshold
    lt_read_threshold_breaks = [
        status_checker.is_ignorelisted(anchor),                         # if anchor is ignorelisted
        not anchor_counts.contains(anchor) and read_counter_freeze,     # if we are no longer accumulating new anchors and this is a new anchor
        not anchor_counts.contains(anchor) and anchor_counter_freeze    # if we are no longer accumulating new anchors and this is a new anchor
    ]

    # break condition if at read threshold
    gt_read_threshold_breaks = [
        not anchor_counts.contains(anchor)                              # if this is a new anchor
    ]

    # assign correct break condition
    anchor_break_conditions = gt_read_threshold_breaks if read_counter_freeze else lt_read_threshold_breaks

    if any(anchor_break_conditions):
        return False
    else:
        return True


def is_phase_0(anchor, target_counts_threshold, anchor_targets):
    """
    Return True if we are not yet at the target_counts_threshold threshold, and we are still
    accumulating targets for this anchor
    """
    if anchor not in anchor_targets or anchor_targets.num_targets(anchor) < target_counts_threshold:
        return True
    else:
        return False


def is_phase_1(anchor, anchor_status, anchor_targets, target_counts_threshold, anchor_counts, anchor_counts_threshold):
    """
    Return True if the anchor qualifies for phase 1
    """
    # define conditions to enter phase_1
    phase_1_conditions = [
        anchor_targets.num_targets(anchor) == target_counts_threshold,
        anchor_counts.get_total_counts(anchor) > anchor_counts_threshold
    ]

    if not anchor_status.is_phase_2(anchor) and any(phase_1_conditions):
        return True
    else:
        return False


def get_iteration_summary_scores(
    score_type,
    iteration,
    read_chunk,
    kmer_size,
    lookahead,
    num_reads,
    anchor_counts,
    anchor_targets_samples,
    anchor_targets,
    anchor_topTargets_scores,
    anchor_topTargets_distances,
    anchor_status,
    status_checker,
    read_counter_freeze,
    target_counts_threshold,
    anchor_counts_threshold,
    anchor_freeze_threshold,
    anchor_mode,
    window_slide,
    compute_target_distance,
    group_ids_dict,
    distance_type
    ):
    """
    Return the summary scores for the reads of one iteration
    """
    phase_0 = 0
    phase_1 = 0
    phase_2 = 0
    phase_1_ignore_score = 0
    phase_1_compute_score = 0
    phase_2_fetch_distance = 0
    phase_2_compute_distance = 0
    valid_anchor = 0
    invalid_anchor = 0
    ignorelisted_anchor = 0
    new_anchor = 0

    # get reads for this iteration
    for read_tuple in read_chunk:

        # define read_tuple
        read, sample = read_tuple

        # get list of all anchors from each read
        anchor_list = read.get_anchors(anchor_mode, window_slide, lookahead, kmer_size)

        # loop over each anchor in the list of all anchors from each read
        for anchor in anchor_list:

            # if this anchor-target pair is eligible for computation, proceed
            if is_valid_anchor_target(anchor, read_counter_freeze, anchor_counts, anchor_freeze_threshold, status_checker):

                if is_phase_0(anchor, target_counts_threshold, anchor_targets):
                    phase_0 += 1

                # get target
                target = read.get_target(anchor, lookahead, kmer_size)

                # updates, always
                anchor_counts.update_total_counts(anchor, iteration)                # update anchor total counts
                anchor_counts.update_all_target_counts(anchor, sample, iteration)   # update anchor-sample counts for all targets

                # if we have not yet logged this target as a top target and we are not at target_counts_threshold
                if not anchor_targets.is_topTarget(anchor, target):
                    if is_phase_0(anchor, target_counts_threshold, anchor_targets):
                        anchor_targets.update_target(anchor, target)             # accumulate as a top target
                        anchor_targets_samples.update(anchor, target, sample)    # acummulate counts of top target

                # else, this is a top target that we have seen before, accumulate counts
                else:
                    anchor_targets_samples.update(anchor, target, sample)        # acummulate counts of top target

                # if we are at threshold, only accumulate anchor_samples count and anchor total counts
                if is_phase_1(anchor, anchor_status, anchor_targets, target_counts_threshold, anchor_counts, anchor_counts_threshold):

                    phase_1 += 1

                    # compute phase_1 score
                    scores, topTargets_distances, mu = compute_phase_1_scores(
                        anchor_targets_samples.get_anchor_counts_df(anchor),
                        group_ids_dict,
                        kmer_size,
                        distance_type,
                        score_type
                    )

                    # set mu_threshold
                    if distance_type == 'hamming':
                        mu_threshold = 2
                    elif distance_type == 'jaccard':
                        mu_threshold = 0.2

                    # if mu < mu_threshold and we have not entered read_counter_freeze, ignorelist this anchor
                    if mu < mu_threshold:

                        status_checker.update_ignorelist(anchor, read_counter_freeze)

                        phase_1_ignore_score += 1

                    # if mu < mu_threshold, proceed with updates and transition to phase_2
                    else:
                        # add this anchor to dict
                        anchor_topTargets_scores.initialise(anchor, scores)                        # update the topTargets for anchor

                        if score_type == 'slow':

                            # updates
                            anchor_topTargets_distances.update_distances(anchor, topTargets_distances) # update distances for topTargets for anchor
                            anchor_topTargets_distances.update_mu(anchor, mu)

                            # after phase_1/transition score computation, assign this anchor to phase_2 and only compute phase_2 score for this anchor
                            anchor_status.assign_phase_2(anchor)

                            phase_1_compute_score += 1


                # compute phase_2 score
                elif anchor_status.is_phase_2(anchor):

                    if score_type == 'slow':

                        phase_2 += 1

                        # if this is a target that we have seen before
                        if anchor_topTargets_distances.has_distance(anchor, target):

                            # get its previously-recorded distance
                            distance = anchor_topTargets_distances.get_distance(anchor, target)
                            phase_2_fetch_distance += 1

                        # if we have never seen this target before
                        else:

                            # compute target distance from topTargets
                            if compute_target_distance:
                                distance = anchor_topTargets_scores.compute_target_distance(anchor, target, distance_type)
                            else:
                                distance = 1

                            phase_2_compute_distance += 1

                        # compute phase_2 score
                        new_score = compute_phase_2_scores(
                            anchor_topTargets_scores.get_score(anchor, sample),
                            anchor_counts.get_all_target_counts(anchor, sample),
                            distance,
                            anchor_topTargets_distances.get_mu(anchor),
                            group_ids_dict[sample]
                        )

                        # update the score for this anchor
                        anchor_topTargets_scores.update(anchor, sample, new_score)

                valid_anchor += 1
            else:
                invalid_anchor += 1

    # return values for logging
    return phase_0, phase_1, phase_2, phase_1_compute_score, phase_1_ignore_score, phase_2_fetch_distance, phase_2_compute_distance, valid_anchor, invalid_anchor, ignorelisted_anchor, new_anchor



