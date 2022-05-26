#!/usr/bin/env python3

import argparse
import pandas as pd

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

    df[['sample', 'anchor']] = df['sample_anchor'].str.split('____', expand=True)

    df = df.drop('sample_anchor', axis=1)

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

    df = pd.merge(df, consensus, on='anchor', how='left')

    anchor_ann = pd.read_csv(args.genome_annotations_anchors, sep='\t', usecols=['anchor', 'local_gene', 'end_to_end_gene'])

    df = pd.merge(df, anchor_ann, on='anchor', how='left')

    df['position'] = df['start'] - 1

    df = df.drop(['start', 'end'], axis=1)

    df.columns = ['called_exon_chr', 'called_exon_id', 'ann_exon_strand', 'ann_exon_gene', 'ann_exon_AS_start', 'ann_exon_AS_end', 'sample', 'anchor', 'consensus', 'anchor_local_gene', 'anchor_end_to_end_gene', 'called_exon_id_position']
    df = df[['called_exon_chr', 'called_exon_id', 'called_exon_id_position', 'ann_exon_strand', 'ann_exon_gene', 'ann_exon_AS_start', 'ann_exon_AS_end', 'sample', 'anchor', 'consensus', 'anchor_local_gene', 'anchor_end_to_end_gene']]

    df.to_csv(args.outfile, sep='\t', index=False)


main()
