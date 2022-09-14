#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np
import os
from Bio import SeqIO

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--fasta_name",
        type=str
    )
    args = parser.parse_args()
    return args

def main():
    args = get_args()
    ## use biopython here
    reads = []

    for record in SeqIO.parse(f"anchor_{args.fasta_name}.fasta", 'fasta'):
        reads.append(str(record.id))
    df = pd.DataFrame(reads, columns=['anchor'])


    ann_tuples = [
        ["strand", f"{args.fasta_name}_genome_strand.txt"],
        ["gene", f"{args.fasta_name}_genome_genes.txt"],
        ["gene_MAPQ", f"{args.fasta_name}_genome_mapq.txt"],
        ["transcript", f"{args.fasta_name}_transcriptome_hit.txt"],
        ["transcriptome_MAPQ", f"{args.fasta_name}_transcriptome_mapq.txt"],
    ]
    print(df.head())
    print("_____")
    for ann_name, ann_file in ann_tuples:
        try:
            ann = pd.read_csv(ann_file, sep='\t', header=None, names=[args.fasta_name, ann_name])
            print(ann.head())
            # clean up
            if "MAPQ" in ann_name:
                ann = ann.replace(0.0, np.nan)

            df = pd.merge(
                df,
                ann,
                left_on='anchor',
                right_on=args.fasta_name,
                how='outer'
            )
        except:
            df[ann_name] = np.nan

    df = df.replace('*', np.nan)

    df.to_csv(f"genome_annotations_anchors.tsv", index=False, sep='\t')


main()
