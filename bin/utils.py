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
import time
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

        # hardcode kmer size for now
        k = 7
        set1 = set([seq1[i:i+k] for i in range(len(seq1)-k+1)])
        set2 = set([seq2[i:i+k] for i in range(len(seq2)-k+1)])

        distance = len(set1.intersection(set2)) / len(set1.union(set2))

    return distance


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
    p['i'] = anchor_counts['distance'].unique()
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

    if score_type == 'slow' or score_type == "force":
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


def is_valid_anchor_target(anchor, target, read_counter_freeze, anchor_counts, anchor_freeze_threshold, status_checker):
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

    # only return True if target is valid and none of the anchor break conditions are true
    if target and not any(anchor_break_conditions):
        return True
    else:
        return False


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


def get_sample_dict(sample_list, samples, use_std):
    """
    Return group_ids_dict and sample_index_dict
    """
    # if not using standard deviation score, make dict of {sample : group_id}
    group_ids_dict = {}
    if not use_std or sample_list.shape[1] != 1:
        group_ids = sample_list.iloc[:,1].tolist()
        for i in range(0, len(samples)):
            group_ids_dict[samples[i]] = group_ids[i]
    else:
        for i in range(0, len(samples)):
            group_ids_dict[samples[i]] = 1

    # create sample index dict, of {sample : sample_id}
    sample_index_dict = {}
    for i in range(len(samples)):
        sample_index_dict[samples[i]] = i

    return group_ids_dict, sample_index_dict


def ignorelist_scores(anchor_topTargets_scores, status_checker, use_std, read_counter_freeze):
    """
    Update ignorelist with an anchor if its scores are in the bottom 20% of either
    std(data has no metadata) or sum(data has metadata)
    """
    # ignorelist condition : ignorelist anchors that are in the bottom 20% of abs(sum(scores))
    ignorelist_anchors_score = anchor_topTargets_scores.get_blacklist_anchors(use_std)

    # update ignorelist
    for anchor in ignorelist_anchors_score:
        status_checker.update_ignorelist(anchor, read_counter_freeze)

    logging.info(f'\t\tanchors that have failed the score requirement = {len(ignorelist_anchors_score)}')

    return


def ignorelist_abundance(anchor_counts, status_checker, num_reads, c_type, samples, read_counter_freeze):
    """
    Update ignorelist with an anchor if its coutns do not meet an abundance requirement
    """
    # intialise
    ignore_abundance = 1

    # define
    k = math.floor(num_reads / 100000)
    if c_type == "num_samples":
        c = len(samples)
    elif c_type == "constant":
        c = 1

    # define min number of counts required per anchor to prevent ignorelisting
    anchor_min_count = k * c

    # get the anchors that do not meet the abundance requirement of anchor_min_count
    ignorelist_anchors_abundance = anchor_counts.get_ignorelist_anchors(anchor_min_count)

    # add those anchors to the ignorelist
    for anchor in ignorelist_anchors_abundance:
        status_checker.update_ignorelist(anchor, read_counter_freeze)
        ignore_abundance += 1

    # logging: update sizes for abundance and score ignorelisting
    logging.info(f'\t\tanchors that have failed the abundance requirement = {len(ignorelist_anchors_abundance)}')
    logging.info(f'\t\t\tabundance requirement = {anchor_min_count} minimum total anchors')

    return



def log_params(args, use_std):
    """
    Output input parameters
    """
    logging.info("--------------------------------------Parameters--------------------------------------")
    logging.info("--------------------------------------------------------------------------------------")
    logging.info(f'samplesheet = {args.samplesheet}')
    logging.info(f'Number of targets per anchor required to calculate phase_1 score = {args.target_counts_threshold}')
    logging.info(f'Number of total anchors required to calculate phase_1 score = {args.anchor_counts_threshold}')
    logging.info(f'Keeplist condition: top {args.num_keep_anchors} anchor significance scores if more than {args.anchor_score_threshold} anchors with scores')
    logging.info(f'Ignorelist condition: phase 1 scores with mean(S_i) < 3')
    logging.info(f'Ignorelist condition: if more than {args.anchor_score_threshold} anchors with scores, scores in 40-60% percentile')
    logging.info(f'Ignorelist condition: computing min_anchor_counts with c_type = {args.c_type}')
    logging.info(f'Anchor freeze threshold = {args.anchor_freeze_threshold} anchors')
    logging.info(f'Read freeze threshold = {args.read_freeze_threshold} anchors')
    logging.info('')

    logging.info(f'n_iterations             = {args.n_iterations}')
    logging.info(f'chunk_size               = {args.max_reads}')
    logging.info(f'kmer_size                = {args.kmer_size}')
    logging.info(f'target_counts_threshold  = {args.target_counts_threshold}')
    logging.info(f'anchor_counts_threshold  = {args.anchor_counts_threshold}')
    logging.info(f'anchor_freeze_threshold  = {args.anchor_freeze_threshold}')
    logging.info(f'read_freeze_threshold    = {args.read_freeze_threshold}')
    logging.info(f'anchor_mode              = {args.anchor_mode}')
    if args.anchor_mode == 'tile':
        logging.info(f'window_slide             = {args.window_slide}')
    logging.info(f'c_type                   = {args.c_type}')
    logging.info(f'lookahead                = {args.lookahead}')
    logging.info(f'num_keep_anchors         = {args.num_keep_anchors}')
    logging.info(f'use_std                  = {use_std}')
    logging.info(f'compute_target_distance  = {args.compute_target_distance}')
    logging.info(f'distance_type            = {args.distance_type}')
    logging.info(f'score_type               = {args.score_type}')
    logging.info("--------------------------------------------------------------------------------------")
    logging.info("--------------------------------------------------------------------------------------")
    logging.info("")
    logging.info("")
    return None


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
    # logging: intialise
    start_time = time.time()
    phase_0 = 0
    phase_1 = 0
    phase_2 = 0
    phase_1_ignore_score = 0
    phase_1_compute_score = 0
    phase_2_fetch_distance = 0
    phase_2_compute_distance = 0
    valid_anchor = 0
    invalid_anchor = 0

    # get reads for this iteration
    for read_tuple in read_chunk:

        # define read_tuple
        read, sample = read_tuple

        # get list of all anchors from each read
        anchor_list = read.get_anchors(anchor_mode, window_slide, lookahead, kmer_size)

        # loop over each anchor in the list of all anchors from each read
        for anchor in anchor_list:

            # get target
            target = read.get_target(anchor, lookahead, kmer_size)

            # if this anchor-target pair is eligible for computation, proceed
            if is_valid_anchor_target(anchor, target, read_counter_freeze, anchor_counts, anchor_freeze_threshold, status_checker):

                # logging
                if is_phase_0(anchor, target_counts_threshold, anchor_targets):
                    phase_0 += 1

                # updates, always
                anchor_counts.update_total_counts(anchor, iteration)
                anchor_counts.update_all_target_counts(anchor, sample, iteration)

                # if we are in phase 0 and we have not yet logged this target as a top target,
                # first accumulate this target as a top target and then accumulate its count (the order is important)
                if not anchor_targets.is_topTarget(anchor, target):
                    if is_phase_0(anchor, target_counts_threshold, anchor_targets):
                        anchor_targets.update_target(anchor, target)
                        anchor_targets_samples.update(anchor, target, sample)

                # else, this is a top target that we have seen before, simply accumulate counts
                else:
                    anchor_targets_samples.update(anchor, target, sample)

                # if we are at threshold, only accumulate anchor_samples count and anchor total counts
                if is_phase_1(anchor, anchor_status, anchor_targets, target_counts_threshold, anchor_counts, anchor_counts_threshold):

                    # logging
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

                        # update ignorelist
                        status_checker.update_ignorelist(anchor, read_counter_freeze)

                        # logging
                        phase_1_ignore_score += 1

                    # if mu < mu_threshold, proceed with updates and transition to phase_2
                    else:
                        # add this anchor to dict
                        anchor_topTargets_scores.initialise(anchor, scores)

                        if score_type == 'slow':

                            # updates
                            anchor_topTargets_distances.update_distances(anchor, topTargets_distances)
                            anchor_topTargets_distances.update_mu(anchor, mu)

                            # after phase_1/transition score computation, assign this anchor to phase_2 and only compute phase_2 score for this anchor
                            anchor_status.assign_phase_2(anchor)

                            # logging
                            phase_1_compute_score += 1


                # compute phase_2 score
                elif anchor_status.is_phase_2(anchor):

                    if score_type == 'slow':

                        # logging
                        phase_2 += 1

                        # if this is a target that we have seen before
                        if anchor_topTargets_distances.has_distance(anchor, target):

                            # get its previously-recorded distance
                            distance = anchor_topTargets_distances.get_distance(anchor, target)

                            # logging
                            phase_2_fetch_distance += 1

                        # if we have never seen this target before
                        else:

                            # compute target distance from topTargets
                            if compute_target_distance:
                                distance = anchor_topTargets_scores.compute_target_distance(anchor, target, distance_type)
                            else:
                                distance = 1

                            # logging
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

                # logging
                valid_anchor += 1

            else:
                # logging
                invalid_anchor += 1

    # logging: get total run time
    run_time = time.time() - start_time

    # logging: get anchor percentages
    try:
        valid_anchor_percent = round((valid_anchor / (valid_anchor + invalid_anchor)) * 100, 2)
        invalid_anchor_percent = round((invalid_anchor / (valid_anchor + invalid_anchor)) * 100, 2)
        phase_0_perecent = round((phase_0 / (phase_0+phase_1+phase_2)) * 100, 2)
        phase_1_percent = round((phase_1 / (phase_0+phase_1+phase_2)) * 100, 2)
        phase_2_percent = round((phase_2 / (phase_0+phase_1+phase_2)) * 100, 2)

    except:
        valid_anchor_percent = 0
        invalid_anchor_percent = 0
        phase_0_perecent = 0
        phase_1_percent = 0
        phase_2_percent = 0

    # logging: output for this iteration
    logging.info("")
    logging.info("-----------------------------------------------------------------------------------------------------------------")
    logging.info("")
    logging.info(f'i = {iteration}')
    logging.info(f'\tReads procesed per file = {num_reads}')
    logging.info(f'\tReads processed total = {len(read_chunk)}')
    logging.info(f'\tRead_counter_freeze = {read_counter_freeze}')
    logging.info("")
    logging.info(f'\tIteration run time = {round(run_time, 2 )} seconds')
    logging.info("")
    logging.info(f'\tinvalid anchors = {invalid_anchor} ({invalid_anchor_percent}%)')
    logging.info("")
    logging.info(f'\tvalid anchors = {valid_anchor} ({valid_anchor_percent}%)')
    logging.info(f'\t\tphase_0 anchors = {phase_0} ({phase_0_perecent}%)')
    logging.info(f'\t\tphase_1 anchors = {phase_1} ({phase_1_percent}%)')
    logging.info(f'\t\t\tscore passed diversity condition = {phase_1_compute_score}')
    logging.info(f'\t\t\tscore failed diversity condition = {phase_1_ignore_score}')
    logging.info(f'\t\tphase_2 anchors = {phase_2} ({phase_2_percent}%)')
    logging.info(f'\t\t\ttarget distance score is fetched = {phase_2_fetch_distance}')
    logging.info(f'\t\t\ttarget distance score is computed = {phase_2_compute_distance}')
    logging.info("")
    logging.info(f'\tnumber of anchors with candidate scores = {len(anchor_topTargets_scores)}')
    logging.info(f'\t\tsize of anchor_counts dict = {len(anchor_counts.total_counts)}')
    logging.info(f'\t\tsize of all_target_counts dict = {len(anchor_counts.all_target_counts)}')
    logging.info(f'\t\tsize of anchor_targets_samples dict = {len(anchor_targets_samples.counter)}')
    logging.info(f'\t\tsize of anchor_topTargets_scores dict = {len(anchor_topTargets_scores)}')
    logging.info(f'\t\tsize of anchor_topTargets_distances dict = {len(anchor_topTargets_distances)}')
    logging.info("")
    logging.info(f'\tignorelist size = {len(status_checker.ignorelist)}')
    logging.info(f'\t\tanchors that have failed the diversity requirement = {phase_1_ignore_score}')

