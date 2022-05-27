#!/usr/bin/env python3

import argparse
import gzip
import argparse
import pandas as pd
import os
import glob


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--outfile_ann_anchors",
        type=str
    )
    parser.add_argument(
        "--outfile_ann_targets",
        type=str
    )
    args = parser.parse_args()
    return args


def main():
    args = get_args()

    # read in anchor_hits
    anchor_hits_paths = glob.glob("anchor_hits*tsv")

    anchor_anns = pd.read_csv(anchor_hits_paths[0].strip(), sep='\t')
    # iteratively merge in hits on anchor column
    for hits_path in anchor_hits_paths[1:]:
        hits = pd.read_csv(hits_path.strip(), sep='\t')
        anchor_anns = anchor_anns.merge(hits, on='anchor')

    anchor_anns.to_csv(
        args.outfile_ann_anchors,
        sep='\t',
        index=False
    )

    # read in target_hits
    target_hits_paths = glob.glob("target_hits*tsv")

    target_anns = pd.read_csv(target_hits_paths[0].strip(), sep='\t')
    # iteratively merge in hits on taret column
    for hits_path in target_hits_paths[1:]:
        hits = pd.read_csv(hits_path.strip(), sep='\t')
        target_anns = target_anns.merge(hits, on='target')

    target_anns.to_csv(
        args.outfile_ann_targets,
        sep='\t',
        index=False
    )



main()
