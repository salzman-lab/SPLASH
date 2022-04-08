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

    args = parser.parse_args()
    return args


def main():
    args = get_args()

    # read in samples from samplesheet
    # [fastq_file, optional group_id]
    sample_list = pd.read_csv(
        args.samplesheet,
        header=None
    )

    # redefine bool
    if args.use_std == "true":
        use_std = True
    elif args.use_std == "false":
        use_std = False
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

    """logging"""
    logging.info("--------------------------------------Parameters--------------------------------------")
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
    logging.info(f'anchor_score_threshold   = {args.anchor_score_threshold}')
    logging.info(f'anchor_mode              = {args.anchor_mode}')
    if args.anchor_mode == 'tile':
        logging.info(f'window_slide             = {args.window_slide}')
    logging.info(f'c_type                   = {args.c_type}')
    logging.info(f'lookahead                = {args.lookahead}')
    logging.info(f'num_keep_anchors         = {args.num_keep_anchors}')
    logging.info(f'use_std                  = {use_std}')
    logging.info(f'compute_target_distance  = {args.compute_target_distance}')
    logging.info(f'distance_type            = {args.distance_type}')

    logging.info("--------------------------------------Parameters--------------------------------------")
    logging.info('')
    """logging"""

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

    # if not using standard deviation score, make dict of {sample : group_id}
    group_ids_dict = {}
    if not use_std or sample_list.shape[1] != 1:
        group_ids = sample_list.iloc[:,1].tolist()
        for i in range(0, len(samples)):
            group_ids_dict[samples[i]] = group_ids[i]
    else:
        for i in range(0, len(samples)):
            group_ids_dict[samples[i]] = 1

    # create sample index dict
    sample_index_dict = {}
    for i in range(len(samples)):
        sample_index_dict[samples[i]] = i

    # initialise objects
    anchor_counts = Anchors.AnchorCounts(len(samples))                          # {anchor : counts}
    anchor_targets_samples = Anchors.AnchorTargetsSamples(sample_index_dict)    # {anchor : {target : [sample_1_count, ...]}}
    anchor_targets = Anchors.AnchorTargets()                                    # {anchor : targets}
    anchor_topTargets_scores = Anchors.AnchorTopTargetsScores()                 # {anchor : scores]}
    anchor_topTargets_distances = Anchors.AnchorTopTargetsDistances()                   # {anchor : {targets : distance}}
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

    # initialise dataframe of logging values
    stats = pd.DataFrame(
        columns = [
            'Run Time', 'Total Reads', 'Anchors with Scores',
            'Phase 1 Total', 'Phase 1 Kept', 'Phase 1 Ignored',
            'Phase 2 Total', 'Phase 2 topTargets', 'Phase 2 Other Targets',
            'anchor_counts', 'anchor_targets_samples', 'anchor_topTargets_scores', 'anchor_topTargets_distances',
            'ignorelist total', 'ignorelist abundance'
        ]
    )

    # begin iterations
    for iteration in range(1, args.n_iterations+1):

        # only continue accumulations if we have less than num_keep_anchors anchors with candidate scores
        if len(anchor_topTargets_scores) >= args.anchor_score_threshold:
            break

        # initialise loggging counters
        ignore_summary_scores = 0
        ignore_abundance = 0

        # if we are at more than read_freeze_threshold reads, set the read_counter freeze so that we don't add any more anchors
        # if we are at more than read_freeze_threshold reads, clear the ignorelist
        num_reads = iteration * args.max_reads
        if num_reads >= args.read_freeze_threshold:
            read_counter_freeze = True
            status_checker.ignorelist.clear()
        else:
            read_counter_freeze = False

        # get reads from this iteration's fastq files
        read_chunk = utils.get_read_chunk(iteration, samples)

        # accumulate for this iteration
        start_time = time.time() # get step start time
        phase_0, phase_1, phase_2, phase_1_compute_score, phase_1_ignore_score, phase_2_fetch_distance, phase_2_compute_distance, valid_anchor, invalid_anchor, ignorelisted_anchor, new_anchor = utils.get_iteration_summary_scores(
            iteration,
            read_chunk,
            args.kmer_size,
            args.lookahead,
            args.max_reads,
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
        run_time = time.time() - start_time # get step total run time

        # ignorelist condition: ignorelist anchors that don't meet an abundance requirement
        k = math.floor(num_reads / 100000)
        if args.c_type == "num_samples":
            c = len(samples)
        elif args.c_type == "constant":
            c = 1
        # define min number of counts required per anchor to prevent ignorelisting
        anchor_min_count = k * c
        # get the anchors that do not meet the abundance requirement of anchor_min_co
        ignorelist_anchors_abundance = anchor_counts.get_ignorelist_anchors(anchor_min_count)
        # add those anchors to the ignorelist
        for anchor in ignorelist_anchors_abundance:
            status_checker.update_ignorelist(anchor, read_counter_freeze)
            ignore_abundance += 1

        # write out ignorelist
        with open(f'ignorelist_iteration_{iteration}.tsv', 'w') as outfile:
            outfile.write("\n".join(status_checker.ignorelist))


        """logging"""
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

        logging.info(f'i = {iteration}')
        logging.info(f'\tReads procesed per file = {num_reads}')
        logging.info(f'\tReads processed total = {num_reads * len(samples)}')
        logging.info(f'\tRead_counter_freeze = {read_counter_freeze}')
        logging.info("")
        logging.info(f'\tIteration run time = {round(run_time, 2 )} seconds')
        logging.info("")
        logging.info(f'\tinvalid anchors = {invalid_anchor} ({invalid_anchor_percent}%)')
        # logging.info(f'\t\tanchor was on the ignorelist = {ignorelisted_anchor}')
        # logging.info(f'\t\tanchor is new and anchor_counter is full = {new_anchor}')
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
        logging.info(f'\t\tanchors that have failed the abundance requirement = {len(ignorelist_anchors_abundance)}')
        logging.info(f'\t\t\tabundance requirement = {anchor_min_count} minimum total anchors')
        logging.info(f'\t\tanchors that have failed the diversity requirement = {phase_1_ignore_score}')
        logging.info("")
        logging.info("-----------------------------------------------------------------------------------------------------------------")

        # add row to stats df
        stats = stats.append(
            {
                'Run Time' : run_time,
                'Total Reads' : num_reads * len(samples),
                'Anchors with Scores': len(anchor_topTargets_scores),
                'Phase 1 Total': phase_1,
                'Phase 1 Kept': phase_1_compute_score,
                'Phase 1 Ignored': phase_1_ignore_score,
                'Phase 2 Total': phase_2,
                'Phase 2 topTargets': phase_2_fetch_distance,
                'Phase 2 Other Targets': phase_2_compute_distance,
                'anchor_counts': len(anchor_counts.total_counts),
                'anchor_targets_samples': len(anchor_targets_samples.counter),
                'anchor_topTargets_scores': len(anchor_topTargets_scores),
                'anchor_topTargets_distances': len(anchor_topTargets_distances),
                'ignorelist total': len(status_checker.ignorelist),
                'ignorelist abundance': ignore_abundance
            },
            ignore_index=True
        )
        """logging"""

    ## done with all iterations ##

    # after all iterations are done accumulating, calculate the summary score
    keep_anchors = anchor_topTargets_scores.get_final_anchors(args.num_keep_anchors)

    # filter for these anchors and drop any duplicates
    final_anchors = (
        pd.DataFrame(keep_anchors, columns = ['anchor'])
        .drop_duplicates()
    )

    ## output this final anchors list
    final_anchors.to_csv(args.outfile, index=False, sep='\t')

    # output run stats
    stats.to_csv('run_stats.tsv', index=False, sep='\t')


main()
