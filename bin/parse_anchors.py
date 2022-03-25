#!/usr/bin/env python3

import pandas as pd
import gzip
import os
import re
import sys
import argparse
import logging


def buildConcensus(kmers, consensus_length, direction):
    """
    Takes a list (for each keyed anchor), and a sequence of bases up to consensus_length;
    Computes base composition of the next seq of bases and outputs
    total num per base, and consensus fraction and base
    """
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


def recordNextKmers(consensus_length, looklength, kmer_size, myseqs, DNAdict, signif_anchors, anchor_dict, fastq_id, direction):
    """
    anchorlist is a list -- we will loopthrough the sequence myseq and check if any of the anchorlist kmers are defined
    anchorlength is length of kmers in file
    LOOK AHEAD IN THE STRING
    """

    j = 0
    ## Consensus anchor step ##
    for myseq in myseqs:

        ## CONSENSUS_ANCHORS STEP ##
        matching_seqs = [a for a in signif_anchors if a in myseq]
        if matching_seqs:
            for anchor in matching_seqs:
                # get consensus candidate position
                if direction == 'down':
                    consensus_seq_start = myseq.index(anchor) + len(anchor)  # end of the anchor
                    consensus_seq_end = consensus_seq_start + consensus_length   # end of the anchor + consensus_length

                if direction == 'up':
                    consensus_seq_end = myseq.index(anchor) # start of the anchor
                    consensus_seq_start = max(0, consensus_seq_end - consensus_length)

                # get consensus candidate
                consensus_seq = myseq[consensus_seq_start:consensus_seq_end]

                # keep dict size small for computational efficiency
                if len(DNAdict[anchor]) < 100:
                    DNAdict[anchor].append(consensus_seq)

        ## ADJ_KMER STEP ##
        matching_anchors = [a for a in signif_anchors if a in myseq]
        if matching_anchors:

            # check for adj_kmers
            for anchor in matching_anchors:
                # get anchor position
                anchor_start = myseq.index(anchor)
                anchor_end = anchor_start + len(anchor)
                # get adjacent anchor position
                if direction == 'down':
                    adj_kmer_start = anchor_end + looklength
                    adj_kmer_end = adj_kmer_start + kmer_size

                if direction == 'up':
                    adj_kmer_end = anchor_start
                    adj_kmer_start = max(0, anchor_end - (kmer_size + looklength))

                adj_kmer = myseq[adj_kmer_start:adj_kmer_end]

                # if adj anchor exists, add adj anchor to anchor_dict
                if len(adj_kmer) == kmer_size:
                    anchor_tuple = (anchor, adj_kmer)

                    if anchor_tuple not in anchor_dict.keys():
                        anchor_dict[anchor_tuple] = 1
                    else:
                        anchor_dict[anchor_tuple] += 1
        j += 1
        if j % 100000 == 1:
            print(j)

    return DNAdict, anchor_dict


def returnSeqs(fastq_file, num_parse_anchors_reads):
    """
    GETTING REAL SEQS
    """

    myseqs = []
    tot_lines = 0

    with gzip.open(fastq_file, 'rt') as handle:
        for read_seq in handle:
            # check we're in sequence line (remainder of 2)
            tot_lines += 1
            if tot_lines%4!=2:
                continue
            # strip of new line character
            read_seq = read_seq.strip('\n')

            if len(myseqs) < num_parse_anchors_reads:
                myseqs.append(read_seq)

    return myseqs


def returnAnchors(infile, direction):
    """
    GETTING REAL SEQS
    """

    anchors=[]
    # Count reads
    tot_lines =  0

    if direction == 'up':
        ind_qval = 5
        ind_seq = 6
    elif direction == 'down':
        ind_qval = 10
        ind_seq = 11

    with gzip.open(infile, "rt") as handle:

        # parse read
        for line in handle:

            # check we're in sequence line (remainder of 2)
            tot_lines += 1
            if tot_lines>1 :
                # strip of new line character
                qval = line.strip().split("\t")[ind_qval]
                seq = line.strip().split("\t")[ind_seq]

                if 'QVAL' in qval: # exclude header lies
                    next
                else:
                    ## edit so that either qup or qdown is less than 0.01
                    if float(qval) < .01 : # lots of clusters
                        anchors.append(seq)
    return list(set(anchors))


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
    parser.add_argument("--out_adj_kmer_file",
        type=str,
        help='output adj_kmer file'
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
        "--looklength",
        help='up or down',
        type=int
    )
    args = parser.parse_args()
    return args


def write_out(nextseqs, consensus_length, direction, out_consensus_fasta_file, out_counts_file, out_fractions_file):
    """
    write out files
    """

    # filse for writing
    outfile_1 = open(out_consensus_fasta_file, "w")
    outfile_2 = open(out_counts_file, "w")
    outfile_3 = open(out_fractions_file, "w")

    for kk in nextseqs.keys():
    # gets the value as an array
    # syntax for getting the values of a key
        if len(nextseqs.get(kk)) > 0 :  # build concensus
            out = buildConcensus(nextseqs.get(kk), consensus_length, direction)
            if len(out[1])>0:
                outfile_1.write(
                    f'>{kk}\n{out[0]}\n'
                )

                str2 = '\t'.join([str(x) for x in out[1]])
                outfile_2.write(
                    f'{kk}\t{out[0]}\t{str2}\n'
                )

                str3 = '\t'.join([str(x) for x in out[2]])
                outfile_3.write(
                    f'{kk}\t{out[0]}\t{str3}\n'
                )

    outfile_1.close()
    outfile_2.close()
    outfile_3.close()


def main():
    args = get_args()

    # get reads from fastq
    myseqs = returnSeqs(
        args.fastq_file,
        args.num_parse_anchors_reads
    )

    # dict for significant anchors and their downstream kmers
    anchor_dict = {}

    # create dict from signif anchors
    signif_anchors_df = (
        pd.read_csv(
            args.anchors_file,
            sep='\t'
        )
    )
    signif_anchors = (
        signif_anchors_df['anchor']
        .drop_duplicates()
        .tolist()
    )

    # DNA dictionary stores the set of reads after each anchor in the anchorlist
    DNAdict = {}
    for an in signif_anchors:
        DNAdict[an] = []

    # get all of the next kmers for the anchors in PREPARATION FOR BUILDING CONCENSUS
    nextseqs, anchor_dict = recordNextKmers(
        args.consensus_length,
        args.looklength,
        args.kmer_size,
        myseqs,
        DNAdict,
        signif_anchors,
        anchor_dict,
        args.fastq_id,
        args.direction
    )

    if anchor_dict:
        # write out anchor dict for merging later
        anchor_df = (
            pd.DataFrame.from_dict(
                anchor_dict,
                orient='index'
            )
            .reset_index()
            .dropna()
        )

        anchor_df.columns = ['anchor_tuple', args.fastq_id]
        anchor_df[['anchor', 'target']] = (
            pd.DataFrame(
                anchor_df['anchor_tuple'].tolist(),
                index=anchor_df.index
            )
        )

        anchor_df = anchor_df[['anchor', 'target', args.fastq_id]]
    else:
        anchor_df = pd.DataFrame.from_dict(anchor_dict)

    anchor_df.to_csv(
        args.out_adj_kmer_file,
        index=False,
        sep='\t'
    )

    write_out(
        nextseqs,
        args.consensus_length,
        args.direction,
        args.out_consensus_fasta_file,
        args.out_counts_file,
        args.out_fractions_file
    )



main()
