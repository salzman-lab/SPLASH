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
import time
import utils
import Anchors
import Reads


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--n_iterations",
        type=int
    )
    parser.add_argument(
        "--max_reads",
        type=int
    )
    parser.add_argument(
        "--kmer_size",
        type=int
    )
    parser.add_argument(
        "--looklength",
        type=int
    )
    parser.add_argument(
        "--samplesheet",
        type=str
    )
    parser.add_argument(
        "--target_counts_threshold",
        type=int
    )
    parser.add_argument(
        "--anchor_counts_threshold",
        type=int
    )
    parser.add_argument(
        "--anchor_freeze_threshold",
        type=int
    )
    parser.add_argument(
        "--anchor_mode",
        type=str
    )
    parser.add_argument(
        "--window_slide",
        type=int
    )
    parser.add_argument(
        "--num_keep_anchors",
        help='up or down',
        type=int
    )
    parser.add_argument(
        "--anchor_score_threshold",
        help='up or down',
        type=int
    )
    parser.add_argument(
        "--c_type",
        type=str
    )
    parser.add_argument(
        "--outfile",
        help='up or down',
        type=str
    )

    args = parser.parse_args()
    return args


def main():
    args = get_args()

    # set up logging
    logging.basicConfig(
        filename = f'get_anchors.log',
        format='%(asctime)s %(levelname)-8s %(message)s',
        level=logging.INFO,
        filemode='w',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

    logging.info('Keeplist if in top or bottom 100 scores')
    logging.info(f'Number of targets required to calculate phase_1 score = {args.target_counts_threshold}')
    logging.info(f'Number of total anchor counts required to calculate phase_1 score = {args.anchor_counts_threshold}')
    logging.info('min(Mean scores) to ignorelist = 3')
    logging.info(f'Percentile of scores to ignorelist = 40-60')
    logging.info(f'Read freeze threshold = 300k reads')
    logging.info(f'Anchor freeze threshold = {args.anchor_freeze_threshold} anchors')
    logging.info('')


    # read in samples
    sample_list = pd.read_csv(args.samplesheet)

    # get list of samples and group_ids
    samples = sample_list.iloc[:,0].tolist()
    group_ids = sample_list.iloc[:,2].tolist()

    # initialise objects
    anchor_counts = Anchors.AnchorCounts(len(samples))          # {anchor : counts}
    anchor_targets_samples = Anchors.AnchorTargetsSamples()     # {anchor : {samples : {targets : counts}}}
    anchor_targets = Anchors.AnchorTargets()                    # {anchor : targets}
    anchor_scores_topTargets = Anchors.AnchorScoresTopTargets() # {anchor : [scores, top_targets]}
    anchor_target_distances = Anchors.AnchorTargetDistances()   # {anchor : {targets : distances}}
    anchor_status = Anchors.AnchorStatus()                      # [anchor1, anchor2]
    status_checker = Anchors.StatusChecker(
        anchor_counts,
        anchor_targets_samples,
        anchor_targets,
        anchor_scores_topTargets,
        anchor_target_distances,
        anchor_status
    )

    # make dict of {sample : group_id}
    group_ids_dict = {}
    for i in range(0, len(samples)):
        group_ids_dict[samples[i]] = group_ids[i]


    for iteration in range(1, args.n_iterations+1):

        # only continue if we have less than 10k anchors with candidate scores
        if anchor_scores_topTargets.get_num_scores() >= 10000:
            break

        # init counters
        ignore_summary_scores = 0
        ignore_abundance = 0

        # if we are at more than 300k reads, don't add any more anchors
        num_reads = iteration * args.max_reads
        if num_reads >= 300000:
            read_counter_freeze = True
            status_checker.ignorelist.clear()
        else:
            read_counter_freeze = False

        # get reads from fastq files
        start_time = time.time()
        read_chunk = utils.get_read_chunk(iteration, samples)
        run_time = time.time() - start_time

        logging.info("")
        logging.info(f'i = {iteration}, {num_reads} reads')
        logging.info(f'\tFetching read chunk = {round(run_time, 2)} seconds')
        logging.info(f'\t\tchunk size = {len(read_chunk)}')

        # get summary scores
        start_time = time.time()
        phase_1, phase_2, ignore_diversity = utils.get_iteration_summary_scores(
            iteration,
            read_chunk,
            args.kmer_size,
            args.looklength,
            args.max_reads,
            group_ids_dict,
            anchor_counts,
            anchor_targets_samples,
            anchor_targets,
            anchor_scores_topTargets,
            anchor_target_distances,
            anchor_status,
            status_checker,
            read_counter_freeze,
            args.target_counts_threshold,
            args.anchor_counts_threshold,
            args.anchor_freeze_threshold,
            args.anchor_mode,
            args.window_slide
        )
        run_time = time.time() - start_time

        # ignorelist condition: ignorelist anchors that don't meet an abundance requirement
        k = math.floor(num_reads / 100000)
        if args.c_type == "num_samples":
            c = len(samples)
        elif args.c_type == "constant":
            c = 1
        anchor_min_count = k * c
        ignorelist_anchors_abundance = anchor_counts.get_ignorelist_anchors(anchor_min_count)
        for anchor in ignorelist_anchors_abundance:
            status_checker.update_ignorelist(anchor)
            ignore_abundance += 1

        # write out ignorelist
        with open(f'ignorelist_iteration_{iteration}.tsv', 'w') as outfile:
            outfile.write("\n".join(status_checker.ignorelist))

        # logging
        logging.info(f'\tIteration run time = {round(run_time, 2 )} seconds')
        logging.info(f'\t\tnumber of anchors with candidate scores = {anchor_scores_topTargets.get_num_scores()}')
        logging.info(f'\t\tnumber of phase_1 calculations = {phase_1}')
        logging.info(f'\t\tnumber of phase_2 calculations = {phase_2}')
        logging.info(f'\t\tignorelist size = {len(status_checker.ignorelist)}')
        logging.info(f'\t\t\tignore_diversity = {ignore_diversity}')
        logging.info(f'\t\t\tignore_abundance = {ignore_abundance}')
        logging.info(f'\t\t\t\tanchor_min_count = {anchor_min_count}')

    ## done with all iterations ##

    # get summary scores
    summary_scores = anchor_scores_topTargets.get_summary_scores(group_ids_dict)

    # only ignorelist if we have more than anchor_score_threshold anchors with scores
    if anchor_scores_topTargets.get_num_scores() >= args.anchor_score_threshold:

        # ignorelist condition: ignorelist the anchors with scores in [40% quantile, 60% quantile]
        lower = 0.4
        upper = 0.6
        lower_bound = summary_scores['score'].quantile([lower, upper]).loc[lower]
        upper_bound = summary_scores['score'].quantile([lower, upper]).loc[upper]
        ignorelist_anchors_percentile = (
            summary_scores[(summary_scores['score'] > lower_bound) & (summary_scores['score'] < upper_bound)]
            .index
            .to_list()
        )
        summary_scores = summary_scores[~summary_scores['score'].isin(ignorelist_anchors_percentile)]

    # output anchors with the highest scores
    summary_scores['abs_score'] = summary_scores['score'].abs()
    summary_scores = summary_scores.sort_values(
        'abs_score',
        ascending=False
    )
    keep_anchors = (
        summary_scores
        .head(args.num_keep_anchors)
        .index
        .to_list()
    )
    final_anchors = (
        pd.DataFrame(
            keep_anchors,
            columns = ['anchor']
        )
        .drop_duplicates()
    )

    ## return final anchors list
    final_anchors.to_csv(args.outfile, index=False, sep='\t')


main()
