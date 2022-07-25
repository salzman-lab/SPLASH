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
        "--run_type",
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
    args = parser.parse_args()
    return args


def main():
    args = get_args()

    x = 0

    file = gzip.open(args.outfile, 'wb')
    with gzip.open(args.infile, 'rt') as infile:
        for line in infile:
            if args.num_lines != 0:
                if x > args.num_lines:
                    break

            x += 1

            cbc = line.split(' ')[0].split('_')[1]
            umi = line.split(' ')[0].split('_')[2]
            sample = cbc + umi

            if x % 4 == 2:
                read = line.strip()
                file.write(str.encode(f'{sample} {read}\n'))

    file.close()


main()
