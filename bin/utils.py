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
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
from datetime import datetime
import Reads


def get_distance(seq1, seq2):
    """
    Returns the Hamming distance between 2 strings
    """
    distance = 0
    for x,y in zip(seq1, seq2):
        distance = distance + (x != y)

    return distance


def is_ignorelisted(self, anchor):
    """
    Returns true if an anchor has previously ignorelisted
    """
    if anchor in self.ignorelist:
        return True
    else:
        return False


def compute_phase_1_score(anchor, df):
    """
    Computes the transition score for an anchor
    """
    # get list of targets for this anchor
    top_targets = list(df['target'])

    min_dists = []

    # iterate over targets sorted by decreasing abundance
    # get min hamming distance of each target to its preceeding targets, such that
    #   for target in n_1...n_5, get d_n_i = min_(1=j)(distance(n_i, n_j))
    for index, row in df.iterrows():
        sequence = row['target']
        candidates = (
            df['target']
            .loc[0:index-1]
            .to_numpy()
        )

        # if it's the first aka most abundant target, initialise the distance to be 0
        if index == 0:
            min_dist = 0

        # else, get the min distance to all preceeding targets
        else:
            min_dist = min([get_distance(sequence, x) for x in candidates])

        min_dists.append((sequence, min_dist))

    # dataframe of [target, target_distance]
    distance_df = (
        pd.DataFrame(
            min_dists,
            columns=['target', 'distance']
        )
        .set_index("target")
    )
    distance_df.index.name = None

    # get phase_1 score of sum(n_i * d_i) / sum(n_i)

    # sum(n_i * d_i)
    # scale by distances, and get sum of (counts * min_dist) across samples
    numerator = (
        df
        .set_index('target')
        .multiply(
            distance_df['distance'],
            axis=0
        )
        .sum(axis=0)
    )

    # sum(n_i)
    # get sum of counts across samples
    denominator = (
        df
        .set_index('target')
        .sum(axis=0)
    )

    scores = (
        numerator
        .divide(
            denominator,
            fill_value=0,
            axis=0
        )
    )

    return scores, top_targets, min_dists


def compute_phase_2_score(previous_scores, counts_df, distance):
    """
    Returns phase_2 scores, given a distance of a new target
    """
    n = (
        counts_df
        .set_index('target')
        .sum(axis=0)
    )

    new_scores = previous_scores * (n/(n+1)) + (distance/(n+1))

    return new_scores


def get_read_chunk(iteration, samples):
    """
    Returns reads from all files for this iteration
    """
    read_chunk = []
    for sample in samples:

        try:
            fastq_file = f'{sample}_{iteration:02d}.fastq.gz'

            with gzip.open(fastq_file, 'rt') as fastq:
                for seqreacord in SeqIO.parse(fastq, 'fastq'):
                    read = Reads.Read(str(seqreacord.seq))
                    read_chunk.append((read, sample))
        except:
            print(f'{sample}, {iteration}')

    return read_chunk


def get_iteration_summary_scores(
    iteration,
    read_chunk,
    kmer_size,
    looklength,
    num_reads,
    anchor_counts,
    anchor_targets_samples,
    anchor_targets,
    anchor_scores_topTargets,
    anchor_target_distances,
    anchor_status,
    status_checker,
    read_counter_freeze,
    target_counts_threshold,
    anchor_counts_threshold,
    anchor_freeze_threshold,
    anchor_mode,
    window_slide,
    compute_target_distance):
    """
    Return the summary scores for the reads of one iteration
    """
    phase_1 = 0
    phase_2 = 0
    phase_1_ignore_score = 0
    phase_1_compute_score = 0
    phase_2_fetch_distance = 0
    phase_2_compute_distance = 0

    # get reads for this iteration
    for read_tuple in read_chunk:

        # define read_tuple
        read, sample = read_tuple

        # get list of all anchors from each read
        anchor_list = read.get_anchors(kmer_size, anchor_mode, window_slide, looklength, kmer_size)

        # loop over each anchor in the list of all anchors from each read
        for anchor in anchor_list:

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
            break_conditions = gt_read_threshold_breaks if read_counter_freeze else lt_read_threshold_breaks

            # only proceed with calculations if we haven't met any break conditions
            if any(break_conditions):
                break

            else:
                # get the target for this anchor
                target = read.get_target(anchor, looklength, kmer_size)

                # if this target exists, proceed
                if target:

                    # updates
                    anchor_counts.update_total_counts(anchor, iteration)            # update anchor total counts
                    anchor_counts.update_sample_counts(anchor, sample, iteration)   # update anchor-sample counts

                    # if we are not at threshold, accumulate this target as a topTarget and accumulate anchor_samples count
                    if (anchor not in anchor_targets) or (anchor_targets.num_targets(anchor) < target_counts_threshold):

                        # accumulate targets
                        anchor_targets.update_target(anchor, target)                    # update anchor-targets
                        anchor_targets_samples.update(anchor, target, sample)           # update anchor-sample-target counts

                    # if we are at threshold, only accumulate anchor_samples count and anchor total counts
                    else:

                        # define conditions to enter phase_1
                        phase_1_conditions = [
                            anchor_targets.num_targets(anchor) == target_counts_threshold,
                            anchor_counts.get_total_counts(anchor) > anchor_counts_threshold
                        ]

                        # compute phase_1 score if not already in phase_2 and meets any phase_1_condition
                        if not anchor_status.is_phase_2(anchor) and any(phase_1_conditions):

                            phase_1 += 1

                            # compute phase_1 score
                            scores, top_targets, topTargets_distances = compute_phase_1_score(
                                anchor,
                                anchor_targets_samples.get_anchor_counts_df(anchor)
                            )

                            # if mean(S_i) < 3 and we have not entered read_counter_freeze, ignorelist this anchor
                            if scores.mean() < 3:

                                status_checker.update_ignorelist(anchor, read_counter_freeze)

                                phase_1_ignore_score += 1

                            # if mean(S_i) >= 3, proceed with updates and transition to phase_2
                            else:

                                # updates
                                anchor_scores_topTargets.update(anchor, scores, top_targets)                      # update the topTargets for anchor
                                anchor_target_distances.update_topTargets_distances(anchor, topTargets_distances) # update distances for topTargets for anchor

                                # after phase_1/transition score computation, assign this anchor to phase_2 and only compute phase_2 score for this anchor
                                anchor_status.assign_phase_2(anchor)

                                phase_1_compute_score += 1

                        # compute phase_2 score
                        else:

                            phase_2 += 1

                            # if this is a target that we have seen before
                            if anchor_target_distances.has_distance(anchor, target):

                                # get its previously-recorded distance
                                distance = anchor_target_distances.get_distance(anchor, target)
                                phase_2_fetch_distance += 1

                            # if we have never seen this target before
                            else:

                                # compute target distance from topTargets
                                if compute_target_distance:
                                    distance = anchor_scores_topTargets.compute_target_distance(anchor, target)
                                else:
                                    distance = 1

                                phase_2_compute_distance += 1

                            # compute phase_2 score
                            new_score = compute_phase_2_score(
                                anchor_scores_topTargets.get_score(anchor),
                                anchor_targets_samples.get_anchor_counts_df(anchor),
                                distance
                            )

                            # updates
                            anchor_scores_topTargets.update(anchor, new_score, anchor_scores_topTargets.get_topTargets(anchor)) # update the score for this anchor

    # return values for logging
    return phase_1, phase_2, phase_1_compute_score, phase_1_ignore_score, phase_2_fetch_distance, phase_2_compute_distance

