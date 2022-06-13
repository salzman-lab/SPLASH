#!/usr/bin/env python3

import argparse
import pandas as pd
from Bio.Seq import Seq
import numpy as np
from Bio.Blast.Applications import NcbiblastnCommandline
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
pd.options.mode.chained_assignment = None  # default='warn'


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--anchor_scores",
        type=str
    )
    parser.add_argument(
        "--anchor_targets_counts",
        type=str
    )
    parser.add_argument(
        "--annotated_anchors",
        type=str
    )
    parser.add_argument(
        "--annotated_targets",
        type=str
    )
    parser.add_argument(
        "--run_blast",
        action='store_true'
    )
    parser.add_argument(
        "--outfile",
        type=str
    )
    args = parser.parse_args()
    return args


def get_blast_df(seqs, seq_type):
    blast_fasta_file = f"blast_{seq_type}.fasta"
    outfile = open(blast_fasta_file, 'w')
    for i in range(len(seqs)):
        outfile.write(f'>seq{i}_{seqs[i]}_NN\n{seqs[i]}\n')
    outfile.close()

    cline = NcbiblastnCommandline(
        outfmt="6 qseqid evalue stitle",
        query=blast_fasta_file,
        remote=True,
        db="nt",
        out="-",
        evalue=0.2,
        task="blastn",
        dust="no",
        word_size=24,
        reward=1,
        penalty=-3,
        max_target_seqs=200
    )

    # reformat blastn output to dataframe
    output = cline()[0].strip()
    rows = [line.split('\t') for line in output.splitlines()]

    blast_df = pd.DataFrame(rows, columns=['qseqid', 'evalue', 'stitle'])

    blast_df[seq_type] = blast_df['qseqid'].str.split('_').str[1]

    blast_df = (
        blast_df
        .set_index(seq_type)
        .drop("qseqid", axis=1)
    )

    # get top annotation with evalue
    blast_df[f"{seq_type}_top_ann"] = (
        blast_df
        .groupby(seq_type)['evalue', 'stitle']
        .apply(max)['stitle']
    )

    # get number of annotations
    blast_df[f"{seq_type}_num_ann"] = (
        blast_df
        .groupby(seq_type)['stitle']
        .apply(
            lambda x:
            x.value_counts()[0]
        )
    )

    # clean up
    blast_df = (
        blast_df
        .drop(['evalue', 'stitle'], axis=1)
        .reset_index()
        .drop_duplicates()
    )

    blast_df[f"{seq_type}_annotation_source"] = 'blastn'

    return blast_df


def get_ann(row, seq_type):
    priority_list = [
        f'{seq_type}_hits_UniVec',
        f'{seq_type}_hits_illumina_adapters',
        f'{seq_type}_hits_escherichia_phage_phiX174',
        f'{seq_type}_hits_RF_all',
        f'{seq_type}_hits_dfam_te_eukaryota',
        f'{seq_type}_hits_tncentral_te_prokaryotes_final',
        f'{seq_type}_hits_mge_aclame_genes_all_0',
        f'{seq_type}_hits_ice_iceberg',
        f'{seq_type}_hits_direct_repeats',
        f'{seq_type}_hits_spacers',
        f'{seq_type}_hits_eukaryota_its1_itstonedb',
        f'{seq_type}_hits_eukaryota_its2_biozentrum',
        f'{seq_type}_hits_WBcel235',
        f'{seq_type}_hits_TAIR10',
        f'{seq_type}_hits_grch38_1kgmaj'
    ]

    for ref in priority_list:
        if row.loc[ref] != '*':
            row.loc[f"{seq_type}_top_ann"] = ref
            row.loc[f"{seq_type}_top_ann_hit"] = row.loc[ref]
            row.loc[f"{seq_type}_top_ann_hit_pos"] = row.loc[ref.replace("hits", "hits_pos")]
            break
    return row


def add_summary(df, ann_table, seq_type, run_blast, top_anchors):
    c_list = [c for c in ann_table.columns if 'hits_pos' not in c]


    ann_table['all_unannotated'] = (
        ann_table[c_list]
        .drop(seq_type, axis=1)
        .replace('*', np.NaN)
        .isnull()
        .all(axis=1)
    )

    # df with bowtie2 ann
    ann_df = ann_table[ann_table['all_unannotated'] == False]

    ann_df[f"{seq_type}_top_ann"] = None
    ann_df[f"{seq_type}_top_ann_hit"] = None
    ann_df[f"{seq_type}_top_ann_hit_pos"] = None

    ann_df = ann_df.apply(
        lambda row:
        get_ann(row, seq_type), axis=1
    )

    ann_df[f"{seq_type}_num_ann"] = (
            ann_df[c_list]
            .drop(seq_type, axis=1)
            .replace('*', np.nan)
            .notnull()
            .astype(int)
            .sum(axis=1)
        )

    ann_df[f"{seq_type}_annotation_source"] = "bowtie2"
    bowtie_df = ann_df[[seq_type, f"{seq_type}_top_ann", f"{seq_type}_top_ann_hit", f"{seq_type}_top_ann_hit_pos", f"{seq_type}_num_ann", f"{seq_type}_annotation_source"]]

    if run_blast:
        # df without any bowtie2 ann
        unann_df = ann_table[ann_table['all_unannotated'] == True]

        if seq_type == 'anchor':
            seqs = (
                pd.merge(
                    unann_df[['anchor']].drop_duplicates(),
                    df[['anchor', 'pv_hash']].drop_duplicates()
                )
                .sort_values('pv_hash', ascending=False)
                .head(2)['anchor']
                .to_list()
            )

        elif seq_type == 'target':
            seqs = (
                df[df['anchor'].isin(top_anchors)]
                .groupby(['anchor'])
                .apply(
                    lambda x:
                    x.sort_values(['total_anchor_target_counts'], ascending = False)
                )
                .reset_index(drop=True)
                .groupby('anchor')
                .head(5)['target']
                .to_list()
            )

        blast_df = get_blast_df(seqs, seq_type)

        addition = pd.concat([blast_df, bowtie_df], axis=0)

    else:
        addition = bowtie_df


    if run_blast:
        if seq_type == 'anchor':
            return addition, seqs
        elif seq_type == 'target':
            return addition, None
    else:
        return addition, None


def main():
    args = get_args()

    ## read in anchors and their scores
    scores = (
        pd.read_csv(args.anchor_scores, sep='\t')
        .drop_duplicates()
    )

    ## account for decoy scores
    if len(scores.columns) == 2:
        scores.columns = ['anchor', 'decoy_pvalue']

    scores_cols = [c for c in scores.columns if "cj_" not in c]
    scores = scores[scores_cols]

    print(scores.head())

    ## read in anchor target counts file
    anchor_targets_counts = pd.read_csv(args.anchor_targets_counts, sep='\t')

    make_target_columns = True
    if len(anchor_targets_counts.columns) == 2:
        anchor_targets_counts.columns = ['anchor', 'target']
        make_target_columns = False

    ## read in annotated anchors and targets files
    ann_anchors = pd.read_csv(args.annotated_anchors, sep='\t')
    ann_targets = pd.read_csv(args.annotated_targets, sep='\t')

    if make_target_columns:
        # get top 5 most abundant targets per anchor by their indices
        top_targets = (
            anchor_targets_counts
            .groupby('anchor')
            .apply(
                lambda x:
                x.nlargest(5, ['total_anchor_target_counts'])
            )
            .drop('anchor', axis=1)
            ['target']
            .tolist()
        )

        ## make a new column of number of samples that have this specific anchor
        anchor_targets_counts['n_samples_anchor'] = (
            anchor_targets_counts
            .drop(['target', 'total_anchor_target_counts'], axis=1)
            .groupby('anchor')
            .transform(sum)
            .astype(bool)
            .sum(axis=1)
        )

        ## make a new column of number of samples that have this specific anchor-target
        anchor_targets_counts['n_samples_anchor_target'] = (
            anchor_targets_counts
            .drop(['anchor', 'target', 'total_anchor_target_counts'], axis=1)
            .astype(bool)
            .sum(axis=1)
            .reset_index(drop=True)
        )

        anchor_targets_counts = anchor_targets_counts[['anchor', 'target', 'total_anchor_target_counts', 'n_samples_anchor', 'n_samples_anchor_target']]

    df = pd.merge(anchor_targets_counts, scores, on='anchor')

    if make_target_columns:
        df['rank_target_counts'] = (
            df
            .groupby('anchor')['total_anchor_target_counts']
            .rank('average', ascending=False)
        )

    df['rcAnchor'] = (
        df['anchor']
        .apply(
            lambda x:
            str(Seq(x).reverse_complement())
        )
    )
    df['rcTarget'] = (
        df['target']
        .apply(
            lambda x:
            str(Seq(x).reverse_complement())
        )
    )

    df['anchor_matches_rcTarget'] = df['rcTarget'].isin(df['anchor'])
    df['rcAnchor_matches_target'] = df['rcAnchor'].isin(df['target'])

    summarized_anchors, top_anchors = add_summary(df, ann_anchors, 'anchor', args.run_blast, None)
    summarized_targets, _ = add_summary(df, ann_targets, 'target', args.run_blast, top_anchors)

    df = pd.merge(df, summarized_anchors, on='anchor', how='left')
    df = pd.merge(df, summarized_targets, on='target', how='left')

    df.to_csv(args.outfile, sep='\t', index=False)


main()
