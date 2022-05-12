#!/usr/bin/env python3

import argparse
from collections import defaultdict
import pandas as pd
import pickle


def get_args():
    parser = argparse.ArgumentParser(description="add annotation columns for splicing")
    parser.add_argument("-i","--in_file",help="the file to add columns to. Must be human data, be tab separated, and have columns chrR1A, chrR1B, juncPosR1A, and juncPosR1B. Will create columns exon_annR1A, exon_annR1B, both_ann, splice_ann, and sort_junc (the last is just an artifact of computation)")
    parser.add_argument("-o","--out_file",help="file to save the output to. If you just want to add the columns to the original file you can pass in the same path as in_file")
    parser.add_argument("-e","--exon_pickle",help="the pickle file for exon annotation")
    parser.add_argument("-s","--splice_pickle",help="the pickle file for splice junction annotation")
    args = parser.parse_args()
    return args


def add_exon_columns(temp_df, exon_bounds):
    for suffix in ["A","B"]:
        temp_df["exon_annR1" + suffix] = False
    for name2, group in temp_df.groupby("chrR1A"):
        temp_df.loc[group.index,"exon_annR1" + suffix] = group["juncPosR1" + suffix].isin(exon_bounds[name2])

    temp_df["both_ann"] = (temp_df["exon_annR1B"] & temp_df["exon_annR1A"]).astype("bool")
    return temp_df


def add_splice_ann_column(temp_df, splices):
    temp_df["sort_junc"] = [tuple(sorted([x,y])) for x, y in zip(temp_df.juncPosR1A,temp_df.juncPosR1B)]
    temp_df["splice_ann"] = False

    for name2, group in temp_df.groupby("chrR1A"):
        sub_group = group[group["chrR1A"].astype(str) == group["chrR1A"].astype(str)]
        if name2 in splices:
            temp_df.loc[sub_group.index,"splice_ann"] = sub_group["sort_junc"].isin(splices[name2])
    return temp_df

def main():
    args = get_args()

    exon_bounds = pickle.load(open(args.exon_pickle, "rb"))
    splices = pickle.load(open(args.splice_pickle, "rb"))

    exon_bounds = defaultdict(set,exon_bounds)
    splices = defaultdict(set,splices)

    df = pd.read_csv(args.in_file, sep = "\t", header=None)
    df.columns = ['chrR1A', 'juncPosR1A', 'juncPosR1B', 'col1',  'col2', 'col3', 'col4', 'col5', 'col6']
    df = add_exon_columns(df, exon_bounds)
    df = add_splice_ann_column(df, splices)
    df.to_csv(args.out_file, sep="\t", index = False)

main()
