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
        "--id",
        type=str
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

            x += 1
            if x % 4 == 1:
                cbc = line.split(' ')[0].split('_')[1]
                umi = line.split(' ')[0].split('_')[2]
                sample = f'{cbc}____{umi}____{args.id}'

            if x % 4 == 2:
                read = line.strip()
                file.write(str.encode(f'{sample} {read}\n'))

    file.close()


main()
