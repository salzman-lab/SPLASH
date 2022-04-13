#!/usr/bin/env python3

import pandas as pd
import gzip
import os
import re
import sys
import argparse
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
        "--lookahead",
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
        "--read_freeze_threshold",
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
    parser.add_argument(
        "--use_std",
        type=str
    )
    parser.add_argument(
        "--compute_target_distance",
        type=str
    )
    parser.add_argument(
        "--max_ignorelist",
        type=int
    )
    parser.add_argument(
        "--distance_type",
        type=str
    )
    parser.add_argument(
        "--score_type",
        type=str
    )

    args = parser.parse_args()
    return args


def main():
    args = get_args()

    # read in samples from samplesheet, formatted as [fastq_file, optional group_id]
    sample_list = pd.read_csv(
        args.samplesheet,
        header=None
    )

    # redefine bool
    if args.use_std == "true":
        use_std = True
    elif args.use_std == "false":
        use_std = False
    # if sampelsheet only has 1 column, force use_std
    if sample_list.shape[1] == 1:
        use_std = True

    # set up logging
    logging.basicConfig(
        filename = f'get_anchors.log',
        format='%(asctime)s %(levelname)-8s %(message)s',
        level=logging.INFO,
        filemode='w',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

    # logging: print input parametrs
    utils.log_params(args, use_std)

    # get list of samples from fastq_files (if fastq_file = "file1.fastq.gz", sample = "file1")
    samples = (
        sample_list
        .iloc[:,0]
        .apply(
            lambda x:
            os.path.basename(x).split('.')[0]
        )
        .tolist()
    )

    # get sample-related dicts
    group_ids_dict, sample_index_dict = utils.get_sample_dict(sample_list, samples, use_std)

    # initialise objects
    anchor_counts = Anchors.AnchorCounts(len(samples))                          # {anchor : counts}
    anchor_targets_samples = Anchors.AnchorTargetsSamples(sample_index_dict)    # {anchor : {target : [sample_1_count, ...]}}
    anchor_targets = Anchors.AnchorTargets()                                    # {anchor : targets}
    anchor_topTargets_scores = Anchors.AnchorTopTargetsScores()                 # {anchor : scores]}
    anchor_topTargets_distances = Anchors.AnchorTopTargetsDistances()           # {anchor : {targets : distance}}
    anchor_status = Anchors.AnchorStatus()                                      # [anchor1, anchor2]
    status_checker = Anchors.StatusChecker(
        anchor_counts,
        anchor_targets_samples,
        anchor_targets,
        anchor_topTargets_scores,
        anchor_topTargets_distances,
        anchor_status,
        args.max_ignorelist
    )

    # begin iterations
    for iteration in range(1, args.n_iterations+1):

        # if we are at more than read_freeze_threshold reads, set the read_counter freeze so that we don't add any more anchors
        # if we are at more than read_freeze_threshold reads, clear the ignorelist
        num_reads = iteration * args.max_reads
        if num_reads >= args.read_freeze_threshold:
            read_counter_freeze = True
            status_checker.ignorelist.clear()
        else:
            read_counter_freeze = False

        # get reads from this iteration's fastq files
        read_chunk = utils.get_read_chunk(iteration, samples, args.n_iterations)

        # accumulate anchors for this iteration
        utils.get_iteration_summary_scores(
            args.score_type,
            iteration,
            read_chunk,
            args.kmer_size,
            args.lookahead,
            args.max_reads,
            samples,
            anchor_counts,
            anchor_targets_samples,
            anchor_targets,
            anchor_topTargets_scores,
            anchor_topTargets_distances,
            anchor_status,
            status_checker,
            read_counter_freeze,
            args.target_counts_threshold,
            args.anchor_counts_threshold,
            args.anchor_freeze_threshold,
            args.anchor_mode,
            args.window_slide,
            args.compute_target_distance,
            group_ids_dict,
            args.distance_type
        )


        # scores threshold : only continue accumulations if we have less than num_keep_anchors anchors with candidate scores
        # only do this iterative version if we are using the slow score
        ignorelist_anchors_score = []
        if args.score_type == "slow" and len(anchor_topTargets_scores) >= args.num_keep_anchors:

            # ignorelist condition : ignorelist anchors that are in the bottom 20% of abs(sum(scores))
            ignorelist_anchors_score = anchor_topTargets_scores.get_blacklist_anchors(use_std)
            # update ignorelist
            for anchor in ignorelist_anchors_score:
                status_checker.update_ignorelist(anchor, read_counter_freeze)

        # abundance threshold: ignorelist anchors that don't meet an abundance requirement
        ignore_abundance = 0
        k = math.floor(num_reads / 100000)
        if args.c_type == "num_samples":
            c = len(samples)
        elif args.c_type == "constant":
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
        logging.info(f'\t\tanchors that have failed the score requirement = {len(ignorelist_anchors_score)}')


    ## done with all iterations ##

    if args.score_type == "fast":
        for anchor in anchor_topTargets_scores:
            # get scores for each anchor
            scores, _, _ = utils.compute_phase_1_scores(
                anchor_targets_samples.get_anchor_counts_df(anchor),
                group_ids_dict,
                args.kmer_size,
                args.distance_type,
                args.score_type
            )
            # updates scores for each anchor
            anchor_topTargets_scores.initialise(anchor, scores)

    # after all iterations are done accumulating, calculate the summary score
    keep_anchors = anchor_topTargets_scores.get_final_anchors(args.num_keep_anchors, use_std)

    # filter for these anchors and drop any duplicates
    final_anchors = (
        pd.DataFrame(keep_anchors, columns = ['anchor'])
        .drop_duplicates()
    )

    ## output this final anchors list
    final_anchors.to_csv(args.outfile, index=False, sep='\t')



main()
