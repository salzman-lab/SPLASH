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
        "--samplesheet",
        type=str
    )
    parser.add_argument(
        "--target_threshold",
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
        "--read_len",
        help='up or down',
        type=int
    )
    parser.add_argument(
        "--outfile",
        help='up or down',
        type=int
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
    logging.info(f'Number of targets required to calculate phase_1 score = {args.target_threshold}')
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
        ignore_summary_scores = 0
        ignore_abundance = 0

        # if we are at more than 300k reads, don't add any more anchors
        num_reads = iteration * args.max_reads
        if num_reads >= 300000:
            read_counter_freeze = True
            # status_checker.ignorelist.clear()
        else:
            read_counter_freeze = False

        logging.info("")
        logging.info(f'i = {iteration}, {num_reads} reads')

        # get reads from fastq files
        start_time = time.time()
        read_chunk = utils.get_read_chunk(iteration, samples)
        run_time = time.time() - start_time
        logging.info(f'\tFetching read chunk = {round(run_time, 2)} seconds')
        logging.info(f'\t\tchunk size = {len(read_chunk)}')

        # get summary scores
        start_time = time.time()
        summary_scores, phase_1, phase_2, ignore_diversity, spacers_enc, before_df_len, after_df_len = utils.get_iteration_summary_scores(
            iteration,
            read_chunk,
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
            args.target_threshold,
            args.anchor_counts_threshold,
            args.anchor_freeze_threshold,
            args.anchor_mode,
            args.window_slide
        )
        run_time = time.time() - start_time
        logging.info(f'\tSummary score calculation = {round(run_time, 2 )} seconds')
        logging.info(f'\t\ttotal anchors = {len(anchor_counts.counter)}')
        logging.info(f'\t\tanchors with summary scores = {len(summary_scores)}')
        logging.info(f'\t\tnumber of phase_1 calculations = {phase_1}')
        logging.info(f'\t\tnumber of phase_2 calculations = {phase_2}')

        # get keeplist and ignorelist
        keeplist = status_checker.keeplist
        ignorelist = status_checker.ignorelist

        # keeplist condition: if we have more than 1k anchors with scores, keeplist the top and bottom 100 anhors
        if len(summary_scores) > 1000:
            summary_scores_anchors = summary_scores.index.to_list()
            keeplist_anchors = summary_scores_anchors[:100] + summary_scores_anchors[-100:]
            for anchor in keeplist_anchors:
                status_checker.update_keeplist(anchor)
                status_checker.update_ignorelist(anchor)

            keeplist_status = True
        else:
            keeplist_status = False

        if len(summary_scores) > 10000:
            # ignorelist condition: ignorelist the anchors with scores in [40% quantile, 60% quantile]
            lower_bound = summary_scores['score'].quantile([0.4, 0.6]).loc[0.4]
            upper_bound = summary_scores['score'].quantile([0.4, 0.6]).loc[0.6]
            ignorelist_anchors_percentile = (
                summary_scores[(summary_scores['score'] > lower_bound) & (summary_scores['score'] < upper_bound)]
                .index
                .to_list()
            )
            for anchor in ignorelist_anchors_percentile:
                status_checker.update_ignorelist(anchor)
                ignore_summary_scores += 1

            ignore_summary_scores_status = True
        else:
            ignore_summary_scores_status = False

        # ignorelist condition: ignorelist anchors that don't meet an abundance requirement
        k = math.floor(num_reads / 100000)
        # c = len(samples)
        c = 1
        anchor_min_count = k * c
        ignorelist_anchors_abundance = anchor_counts.get_ignorelist_anchors(anchor_min_count)
        for anchor in ignorelist_anchors_abundance:
            status_checker.update_ignorelist(anchor)
            ignore_abundance += 1

        # logging
        logging.info(f'\tKeeplisting/ignorelisting')
        logging.info(f'\t\tkeeplist activated = {keeplist_status}')
        logging.info(f'\t\t\tkeeplist size = {len(status_checker.keeplist)}')
        logging.info(f'\t\tignore_summary_scores activated = {ignore_summary_scores_status}')
        logging.info(f'\t\t\tignorelist size = {len(status_checker.ignorelist)}')
        logging.info(f'\t\t\t\tignore_diversity = {ignore_diversity}')
        logging.info(f'\t\t\t\tignore_summary_scores = {ignore_summary_scores}')
        logging.info(f'\t\t\t\tignore_abundance = {ignore_abundance}')
        logging.info(f'\t\t\t\t\tanchor_min_count = {anchor_min_count}')


        # get keeplist and ignorelist
        keeplist = status_checker.keeplist
        ignorelist = status_checker.ignorelist

        summary_scores.index.name = 'anchor'
        summary_scores.reset_index().to_csv(f'summary_scores_iteration_{iteration}.tsv', index=False, sep='\t')
        logging.info('\tDone with iteration')
    ## done with all iterations ##

    # if keeplist is not empty, use its anchors as outputs
    if keeplist:
        final_anchors_list = keeplist
        with open(args.outfile, "w") as outfile:
            outfile.write("\n".join(keeplist))

    else:
        # if keep_list is empty, then grab the anchors with the top 1000 abs_value(scores)
        final_summary_scores = summary_scores
        final_summary_scores['score'] = abs(final_summary_scores['score'])
        final_anchors = (
            final_summary_scores
            .reset_index()
            .sort_values(
                'score',
                ascending=False
            )
            .head(2000)
            [['anchor']]
        )

        ## return final anchors list
        final_anchors.to_csv(args.outfile, index=False, sep='\t')
        final_anchors_list = final_anchors['anchor'].to_list()


main()
