#!/usr/bin/env python3

import argparse
import os
import gzip

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--infile",
        type=str
    )
    parser.add_argument(
        "--fastq_id",
        type=str
    )
    parser.add_argument(
        "--num_lines",
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
        "--anchor_mode",
        type=str
    )
    parser.add_argument(
        "--window_slide",
        type=int
    )
    parser.add_argument(
        "--outfile",
        type=str
    )
    parser.add_argument(
        "--is_10X",
        action='store_true'
    )
    args = parser.parse_args()
    return args


def get_anchor_target(read, lookahead, kmer_size, step_size):
    last_base = len(read) - (lookahead + 2 * kmer_size)

    for i in range(0, last_base, step_size):
        # get anchor
        anchor = read[0+i : kmer_size+i]

        # get target start and end positions, as a function of anchor end
        target_start = (kmer_size+i) + lookahead
        target_end = target_start + kmer_size

        # get target
        target = read[target_start : target_end]
    # return
    return anchor, target


def main():
    args = get_args()
    print(args.is_10X)

    if args.anchor_mode == 'chunk':
        step_size = args.kmer_size
    elif args.anchor_mode == 'tile':
        step_size = args.window_slide

    x = 0

    file = gzip.open(args.outfile, 'wb')
    with gzip.open(args.infile, 'rt') as infile:
        for line in infile:

            ## if this 10X, parse each line as [cbc+umi, read]
            if args.is_10X:
                (cbc_umi, read) = line.split(" ")
                cbc = cbc_umi[0:16]
                anchor, target = get_anchor_target(read, args.lookahead, args.kmer_size, step_size)

                # write out anchor and target, with the cbc as sample
                if "N" not in anchor and "N" not in target:
                    file.write(str.encode(f'{anchor+target} {cbc}\n'))

            ## if this is not 10X, parse as fastq file
            else:

                if args.num_lines != 0:
                    if x > args.num_lines*4:
                        break

                x += 1
                if x % 4 == 2:
                    read = line.strip()

                    anchor, target = get_anchor_target(read, args.lookahead, args.kmer_size, step_size)

                    if "N" not in anchor and "N" not in target:
                        file.write(str.encode(f'{anchor+target} {args.fastq_id}\n'))
    file.close()


main()
