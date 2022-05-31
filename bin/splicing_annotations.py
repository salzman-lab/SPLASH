#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--infile",
        type=str
    )
    parser.add_argument(
        "--fasta",
        type=str
    )
    parser.add_argument(
        "--genome_annotations_anchors",
        type=str
    )
    parser.add_argument(
        "--consensus_genes",
        type=str
    )
    parser.add_argument(
        "--reported_alignments",
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

    df = pd.read_csv(
        args.infile,
        sep='\t',
        names=['chr', 'start', 'end', 'sample_anchor', 'id', 'gtf_strand', 'gtf_gene', 'gtf_AS_start', 'gtf_AS_end']
    )


    seqs = {}
    with open(args.fasta, 'r') as fasta:
        for line in fasta:
            if line.startswith('>'):
                seqs[line.strip().strip('>')] = next(fasta).strip()

    consensus = (
        pd.DataFrame.from_dict(seqs, orient='index')
        .reset_index()
    )
    consensus.columns = ['anchor', 'consensus']

    df = pd.merge(df, consensus, left_on='sample_anchor', right_on='anchor', how='left')

    consensus_genes = pd.read_csv(args.consensus_genes, sep='\t', names=['sample_anchor', 'consensus_gene'])
    df = pd.merge(df, consensus_genes, on='sample_anchor', how='left')

    df[['sample', 'anchor']] = df['sample_anchor'].str.split('____', expand=True)

    anchor_ann = pd.read_csv(args.genome_annotations_anchors, sep='\t', usecols=['anchor', 'local_gene', 'end_to_end_gene'])
    df = pd.merge(df, anchor_ann, on='anchor', how='left')

    reported_alignments = pd.read_csv(args.reported_alignments, sep='\t', names=['sample_anchor', 'consensus_reported_alignment'])
    df = pd.merge(df, reported_alignments, on='sample_anchor', how='left')

    df['position'] = df['start'] - 1
    df = df.drop(['start', 'end', 'sample_anchor'], axis=1)

    df.rename(
        columns = {
            "chr" : "called_exon_chr",
            "id" : "called_exon_id",
            "gtf_strand" : "ann_exon_strand",
            "gtf_gene" : "ann_exon_gene",
            "gtf_AS_start" : "ann_exon_AS_start",
            "gtf_AS_end" : "ann_exon_AS_end",
            "local_gene" : "anchor_local_gene",
            "end_to_end_gene" : "anchor_end_to_end_gene",
            "position" : "called_exon_id_position"
        },
        inplace=True
    )

    df = df[[
        'called_exon_chr', 'called_exon_id', 'called_exon_id_position', 'ann_exon_gene', 'ann_exon_strand', 'ann_exon_AS_start', 'ann_exon_AS_end',
        'sample', 'anchor', 'consensus', 'anchor_local_gene', 'anchor_end_to_end_gene', 'consensus_gene', 'consensus_reported_alignment']]

    df = df.replace(".", np.nan)

    df.to_csv(args.outfile, sep='\t', index=False)


main()
