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
    args = parser.parse_args()
    return args


def main():
    args = get_args()

    sample = os.path.basename(args.infile).split('.')[0]

    current_substr = None

    with open(args.infile, 'rt') as infile:
        for line in infile:
            try:
                count, kmer, sample = tuple(line.strip().split(" "))
                substr = kmer[0:3]

                # this is the first 3mer, open the file and write out
                if not current_substr:
                    current_substr = substr
                    outfile = open(f"stratified_{substr}.txt", 'a')
                    outfile.write(f"{count} {kmer} {sample}\n")
                else:
                    # we are still on the same 3mer, write out
                    if substr == current_substr:
                        outfile.write(f"{count} {kmer} {sample}\n")
                    # we have encountered a new 3mer, close the old one, open a new one, and write out
                    else:
                        outfile.close()
                        current_substr = substr
                        outfile = open(f"stratified_{substr}.txt", 'a')
                        outfile.write(f"{count} {kmer} {sample}\n")
            except:
                pass


main()

