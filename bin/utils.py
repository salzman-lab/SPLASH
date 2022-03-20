#!/usr/bin/env python

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
import logging
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
    top_targets = list(df['target'].head(5))

    min_dists = []
    for index, row in df.iterrows():
        sequence = row['target']
        candidates = (
            df['target']
            .loc[0:index-1]
            .to_numpy()
        )

        if index == 0:
            min_dist = 0
        else:
            min_dist = min([get_distance(sequence, x) for x in candidates])

        min_dists.append((sequence, min_dist))

    distance_df = (
        pd.DataFrame(
            min_dists,
            columns=['target', 'distance']
        )
        .set_index("target")
    )
    distance_df.index.name = None

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
    Returns phase_2 scores, given a distnace of a new target
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
            pass

    return read_chunk


def get_iteration_summary_scores(
    iteration,
    read_chunk,
    kmer_size,
    num_reads,
    group_ids,
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
    window_slide):
    """
    Return the summary scores for the reads of one iteration
    """
    phase_1 = 0
    phase_2 = 0
    ignore_diversity = 0
    spacers_enc = 0


    logging.info(f'\t\tread_counter_freeze = {read_counter_freeze}')

    # get reads for this iteration
    for read_tuple in read_chunk:

        # define read_tuple
        read, sample = read_tuple

        # get list of all anchors from each read
        # ([anchor1, anchor2, anchor3], fastq_id)
        anchor_list = read.get_anchors(kmer_size, anchor_mode, window_slide)

        # loop over each anchor in the list of all anchors from each read
        for anchor in anchor_list:

            # set anchor_counter_freeze if we have reached the anchor_freeze_threshold
            if len(anchor_counts.counter) >= anchor_freeze_threshold:
                anchor_counter_freeze = True
            else:
                anchor_counter_freeze = False

            lt_read_threshold_breaks = [
                status_checker.is_ignorelisted(anchor), # if anchor is ignorelisted
                not anchor_counts.contains(anchor) and read_counter_freeze, # if we are no longer accumulating new anchors and this is a new anchor
                not anchor_counts.contains(anchor) and anchor_counter_freeze # if we are no longer accumulating new anchors and this is a new anchor
            ]

            gt_read_threshold_breaks = [
                not anchor_counts.contains(anchor) # if this is a new anchor
            ]

            break_conditions = gt_read_threshold_breaks if read_counter_freeze else lt_read_threshold_breaks

            # if not any(break_conditions), continue with calculations
            if any(break_conditions):
                break

            else:
                # get the target for this anchor
                target = read.get_target(anchor, target_dist=0, target_len=27)

                # if this target exists, proceed
                if target:

                    # updates
                    anchor_counts.update(anchor, iteration)               # increment counts for this anchor
                    anchor_targets_samples.update(anchor, target, sample) # increment counts for this anchor:sample:target
                    anchor_targets.update_target(anchor, target)          # update targets for this anchor

                    # for this anchor, get the number of unique targets for thresholding
                    num_targets = len(anchor_targets.get_targets(anchor))

                    # if this anchor only has 5 unique targets, continue and don't do anything else
                    if num_targets < target_counts_threshold:
                        continue

                    # compute phase_1/transition score only if anchor is not in phase_2 AND if either of these are true:
                        # if this anchor has exactly 5 unique targets and has NOT entered phase_2
                        # if this anchor has at least 10 total counts
                    elif not anchor_status.is_phase_2(anchor):
                        if (num_targets == 5) or (anchor_counts.get_count(anchor) > anchor_counts_threshold):

                            # compute phase_1/transition score
                            scores, top_targets, topTargets_distances = compute_phase_1_score(
                                anchor,
                                anchor_targets_samples.get_anchor_counts_df(anchor)
                            )

                            # if mean(S_i) < 3 and we have not entered read_counter_freeze, ignorelist this anchor
                            if scores.mean() < 3:

                                if not read_counter_freeze:
                                    status_checker.update_ignorelist(anchor)
                                    ignore_diversity += 1

                            # if mean(S_i) >= 3, proceed with updates and transition to phase_2
                            else:

                                # updates
                                anchor_scores_topTargets.update(anchor, scores, top_targets)                      # update the topTargets for anchor
                                anchor_target_distances.update_topTargets_distances(anchor, topTargets_distances) # update distances for topTargets for anchor

                                # after phase_1/transition score computation, assign this anchor to phase_2 and only compute phase_2 score for this anchor
                                anchor_status.assign_phase_2(anchor)

                                phase_1 += 1

                    # compute phase_2 score
                    else:

                        # if this is a target that we have seen before
                        if anchor_target_distances.has_distance(anchor, target):

                            # get its previously-recorded distance
                            distance = anchor_target_distances.get_distance(anchor, target)

                        # if we have never seen this target before
                        else:

                            # compute target distance from topTargets
                            distance = anchor_scores_topTargets.compute_target_distance(anchor, target)

                            # updates
                            anchor_target_distances.update_target_distance(anchor, target, distance) # update dict with this new (target:distance) pair

                        # # only do this step if this anchor:target:sample has less than 20 phase_2 score coputations
                        # if not anchor_targets_samples.reached_phase_2_threshold(anchor, target, sample):

                        # compute fast score
                        new_score = compute_phase_2_score(
                            anchor_scores_topTargets.get_score(anchor),
                            anchor_targets_samples.get_anchor_counts_df(anchor),
                            distance
                        )

                        # updates
                        anchor_scores_topTargets.update(anchor, new_score, anchor_scores_topTargets.get_topTargets(anchor)) # update the score for this anchor

                        phase_2 += 1

    # output summary scores
    if anchor_scores_topTargets:
        summary_scores = anchor_scores_topTargets.get_summary_scores(group_ids)

    else:
        summary_scores = pd.DataFrame(columns=['score'])

    return summary_scores, phase_1, phase_2, ignore_diversity

