#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--unmapped_fasta",
        type=str
    )
    parser.add_argument(
        "--ann_called_exons",
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
        "--outfile_unmapped",
        type=str
    )
    parser.add_argument(
        "--outfile_annotations",
        type=str
    )
    args = parser.parse_args()
    return args


def main():
    args = get_args()

    ## read in unmapped fasta reads into df
    unmapped = {}
    with open(args.unmapped_fasta, 'r') as unmapped_fasta:
        for line in unmapped_fasta:
            if line.startswith('>'):
                ## strip away unnecessary information
                unmapped[line.strip().strip('>').split(" ")[0]] = next(unmapped_fasta).strip()
    ## conert to df
    unmapped = (
        pd.DataFrame.from_dict(unmapped, orient='index')
        .reset_index()
    )
    unmapped.columns = ['sample_anchor', 'consensus']
    ## create columns
    unmapped[['sample', 'anchor']] = unmapped['sample_anchor'].str.split('____', expand=True)
    unmapped = unmapped.drop('sample_anchor', axis=1)
    ## output
    unmapped.to_csv(args.outfile_unmapped, index=False, sep='\t')

    ## read in annotated called exons
    df = pd.read_csv(
        args.ann_called_exons,
        sep='\t',
        names=['chr', 'start', 'end', 'sample_anchor', 'id', 'gtf_strand', 'gtf_gene', 'gtf_AS_start', 'gtf_AS_end']
    )

    ## proceed with merging if annotated called exons not empty
    if not df.empty:
        ## read in consensus fasta sequences into df
        seqs = {}
        with open(args.fasta, 'r') as fasta:
            for line in fasta:
                if line.startswith('>'):
                    seqs[line.strip().strip('>')] = next(fasta).strip()
        ## convert to df
        consensus = (
            pd.DataFrame.from_dict(seqs, orient='index')
            .reset_index()
        )
        consensus.columns = ['anchor', 'consensus']

        ## merge in consensus sequences
        df = pd.merge(df, consensus, left_on='sample_anchor', right_on='anchor', how='left')

        ## read in consensus annotated with genes and merge in
        consensus_genes = pd.read_csv(args.consensus_genes, sep='\t', names=['sample_anchor', 'consensus_gene'])
        df = pd.merge(df, consensus_genes, on='sample_anchor', how='left')

        ## create columns
        df[['sample', 'anchor']] = df['sample_anchor'].str.split('____', expand=True)

        ## read in anchor annotations and merge in
        anchor_ann = pd.read_csv(args.genome_annotations_anchors, sep='\t', usecols=['anchor', 'local_gene', 'end_to_end_gene'])
        df = pd.merge(df, anchor_ann, on='anchor', how='left')

        ## try reading in reported alignments, if they exist
        try:
            reported_alignments = pd.read_csv(args.reported_alignments, sep='\t', names=['sample_anchor', 'consensus_reported_alignment'])
        except:
            reported_alignments = pd.DataFrame(columns=['sample_anchor', 'consensus_reported_alignment'])

        ## merge in
        df = pd.merge(df, reported_alignments, on='sample_anchor', how='left')

        ## create actual position coordinate
        df['position'] = df['start'] - 1
        ## remove unnecessary columns
        df = df.drop(['start', 'end', 'sample_anchor'], axis=1)

        ## rename columns
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

        ## reorder columns
        df = df[[
            'called_exon_chr', 'called_exon_id', 'called_exon_id_position', 'ann_exon_gene', 'ann_exon_strand', 'ann_exon_AS_start', 'ann_exon_AS_end',
            'sample', 'anchor', 'consensus', 'anchor_local_gene', 'anchor_end_to_end_gene', 'consensus_gene', 'consensus_reported_alignment']]

        ## clean up from bedtools
        df = df.replace(".", np.nan)

        ## output
        df.to_csv(args.outfile_annotations, sep='\t', index=False)


main()
