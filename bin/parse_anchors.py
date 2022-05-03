#!/usr/bin/env python3

import pandas as pd
import gzip
import os
import re
import sys
import argparse
import logging
from Bio import SeqIO


def build_consensus(kmers, consensus_length, direction):
    """
    Takes a list (for each keyed anchor), and a sequence of bases up to consensus_length;
    Computes base composition of the next seq of bases and outputs
    total num per base, and consensus fraction and base
    """
    print(kmers)
    baseComp = ''   # the most frequent base at the current position
    baseCount = []  # count of most frequent base
    baseFrac = []   # fraction of the most frequence base
    # get length of longest kmer
    longest_kmer_len = len(max(kmers, key=len))

    # iterate over string by position
    for i in range(0, longest_kmer_len):
        # get base at each position for every sequence, if there is a base at that position
        if direction == 'up':
            pos_arr = [x[::-1][i] for x in kmers if len(x) > i]
        if direction == 'down':
            pos_arr = [x[i] for x in kmers if len(x) > i]
        # get the most frequent base, ties broken by alphabetical sort
        base = max(pos_arr, key=pos_arr.count)
        # get count of most frequent base at this position
        count = len(pos_arr)
        # get fraction of occurrence of most frequent base at this position
        frac = (pos_arr.count(base)) / len(pos_arr)

        # update lists for current position
        baseComp+=(base)            # extend string
        baseCount.append(count)     # append to list
        baseFrac.append(frac)       # append to list

    if direction == 'up':
        baseComp = baseComp[::-1]
        baseCount = baseCount[::-1]
        baseFrac = baseFrac[::-1]

    return baseComp, baseCount, baseFrac


def get_targets_consensus_seqs(fastq_file, num_parse_anchors_reads, anchor_dict, consensus_length, kmer_size, lookahead, direction):
    """
    For each read, check if there are any valid targets and/or consensus sequences,
    and append them to their respective dictionaries
    """
    # intialize read counter
    num_reads = 0

    # intialise dict
    target_dict = {}
    consensus_dict = anchor_dict.copy()

    # parse reads from fastq file
    with gzip.open(fastq_file, 'rt') as handle:
        for record in SeqIO.parse(handle, 'fastq'):

            # get read and increment counter
            read = str(record.seq)
            num_reads += 1

            # only proceed if we have not yet reached read threshold
            if num_reads < num_parse_anchors_reads:

                # tile the read by 1bp and check if tiles exist in the anchor dict
                # this is more computationally efficient than checking each read for any anchors
                for i in range(0, len(read)):

                    # define current tile as potential anchor
                    anchor = read[0+i:kmer_size+i]

                    # if potential anchor is in anchor_dict
                    if anchor in anchor_dict:

                        # get anchor start and end positions
                        anchor_start = read.index(anchor)
                        anchor_end = anchor_start + kmer_size

                        # fetch sequences, relative to position
                        if direction == 'down':
                            # get downstream target start and end positions
                            target_start = anchor_end + lookahead
                            target_end = target_start + kmer_size

                            # get downstream consensus start and end positions
                            consensus_seq_start = anchor_end
                            consensus_seq_end = consensus_seq_start + consensus_length

                        if direction == 'up':
                            # get upstream target stand and end positions
                            target_end = anchor_start - lookahead
                            target_start = max(0, anchor_end - (kmer_size + lookahead))

                            # get upstream consensus start and end positions
                            consensus_seq_end = anchor_start
                            consensus_seq_start = max(0, consensus_seq_end - consensus_length)

                        # define target and consensus sequences
                        target = read[target_start : target_end]
                        consensus_seq = read[consensus_seq_start:consensus_seq_end]

                        # if target is large enough, add it to the target_dict
                        if len(target) == kmer_size:
                            anchor_tuple = (anchor, target)

                            if anchor_tuple not in target_dict:
                                target_dict[anchor_tuple] = 1
                            else:
                                target_dict[anchor_tuple] += 1

                        # if we have less than 100 candidate consensus sequences, add it to the consensus_dict
                        if len(consensus_dict[anchor]) < 100:
                            consensus_dict[anchor].append(consensus_seq)

    return target_dict, consensus_dict


def output_consensus(consensus_dict, consensus_length, direction, out_consensus_fasta_file, out_counts_file, out_fractions_file):
    """
    Build consensus sequences and write out
    """

    # Open output file
    outfile_1 = open(out_consensus_fasta_file, "w")
    outfile_2 = open(out_counts_file, "w")
    outfile_3 = open(out_fractions_file, "w")

    # Process for each anchor that contains candidate consensus sequences
    for anchor in consensus_dict.keys():

        if len(consensus_dict.get(anchor)) > 0 :

            # build the consensus
            out = build_consensus(consensus_dict.get(anchor), consensus_length, direction)

            if len(out[1])>0:

                # write out fasta entry
                outfile_1.write(
                    f'>{anchor}\n{out[0]}\n'
                )

                # write out counts entry
                str2 = '\t'.join([str(x) for x in out[1]])
                outfile_2.write(
                    f'{anchor}\t{out[0]}\t{str2}\n'
                )

                # write out fraction entry
                str3 = '\t'.join([str(x) for x in out[2]])
                outfile_3.write(
                    f'{anchor}\t{out[0]}\t{str3}\n'
                )

    # close files
    outfile_1.close()
    outfile_2.close()
    outfile_3.close()


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--num_parse_anchors_reads",
        type=int,
        help='max number of fastq reads for input'
    )
    parser.add_argument(
        "--anchors_file",
        type=str,
        help='input list of all significant anchors'
    )
    parser.add_argument(
        "--anchors_annot",
        type=str,
        help='input anchor file from dgmfinder'
    )
    parser.add_argument("--fastq_file",
        type=str,
        help='input fasta file'
    )
    parser.add_argument("--fastq_id",
        type=str,
        help='fastq id name'
    )
    parser.add_argument("--out_consensus_fasta_file",
        type=str,
        help='output fasta file'
    )
    parser.add_argument("--out_counts_file",
        type=str,
        help='output counts file'
    )
    parser.add_argument("--out_fractions_file",
        type=str,
        help='output fractions file'
    )
    parser.add_argument("--out_target_file",
        type=str,
        help='output target file'
    )
    parser.add_argument("--consensus_length",
        type=int,
        help='lookahead length'
    )
    parser.add_argument(
        "--kmer_size",
        type=int,
        help='kmer size'
    )
    parser.add_argument(
        "--direction",
        type=str,
        help='up or down'
    )
    parser.add_argument(
        "--lookahead",
        help='up or down',
        type=int
    )
    args = parser.parse_args()
    return args


def main():
    args = get_args()

    # set up logging
    logging.basicConfig(
        filename = f'parse_anchors_{args.fastq_id}.log',
        format='%(asctime)s %(levelname)-8s %(message)s',
        level=logging.INFO,
        filemode='w',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

    """logging"""
    logging.info("--------------------------------------Parameters--------------------------------------")
    logging.info(f'consensus_length = {args.consensus_length}')
    logging.info(f'direction        = {args.direction}')
    logging.info(f'kmer_size        = {args.kmer_size}')
    logging.info(f'lookahead       = {args.lookahead}')
    logging.info("--------------------------------------Parameters--------------------------------------")
    logging.info('')
    """logging"""

    # get anchors from file
    anchors = (
        pd.read_csv(args.anchors_file, sep='\t', header=None)
        .iloc[:,0]
        .drop_duplicates()
        .head(10000)
        .tolist()
    )

    # create dict with anchors as keys
    anchor_dict = {a:[] for a in anchors}

    logging.info(f'Starting target fetching')

    # get all of the next kmers for the anchors in PREPARATION FOR BUILDING CONCENSUS
    target_dict, consensus_dict = get_targets_consensus_seqs(
        args.fastq_file,
        args.num_parse_anchors_reads,
        anchor_dict,
        args.consensus_length,
        args.kmer_size,
        args.lookahead,
        args.direction
    )

    logging.info(f'Finished target fetching')

    if not all(map(lambda x: x == [], target_dict.values())):
        # write out anchor dict for merging later
        anchor_df = (
            pd.DataFrame.from_dict(target_dict, orient='index')
            .reset_index()
            .dropna()
        )

        # reformat such that there are is a column of counts per anchor-target,
        # where the column is the fastq_id
        anchor_df.columns = ['anchor_tuple', args.fastq_id]
        anchor_df[['anchor', 'target']] = (
            pd.DataFrame(
                anchor_df['anchor_tuple'].tolist(),
                index=anchor_df.index
            )
        )

        # final output columns
        anchor_df = anchor_df[['anchor', 'target', args.fastq_id]]

    else:
        anchor_df = pd.DataFrame(columns=['anchor', 'target', args.fastq_id])

    # output anchor-target counts for this fastq
    anchor_df.to_csv(args.out_target_file, index=False, sep='\t')

    logging.info(f'Starting consensus building')
    # build and output the consensus files
    output_consensus(
        consensus_dict,
        args.consensus_length,
        args.direction,
        args.out_consensus_fasta_file,
        args.out_counts_file,
        args.out_fractions_file
    )
    logging.info(f'Finished consensus building')


main()


