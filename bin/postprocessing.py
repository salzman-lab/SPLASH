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
        "--anchor_target_counts",
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

def get_blast_df(unann_df, seq_type):
    """
    Function to run blast on the top n sequences without bowtie2 annotations
    """
    # only run blast on top n sequences, as ranked by pvalues
    seqs = (
        unann_df
        .sort_values('bf.cor.p', ascending=False)
        .head(2)[seq_type]
        .to_list()
    )

    # write out anchors without annotations to a fasta file
    blast_fasta_file = f"blast_{seq_type}.fasta"
    outfile = open(blast_fasta_file, 'w')
    for i in range(len(seqs)):
        outfile.write(f'>seq{i}_{seqs[i]}_NN\n{seqs[i]}\n')
    outfile.close()

    # run blastn
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

    # reformat columns if there are blast results
    if not blast_df.empty:
        # get seqs from seq_ids
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
        )

        blast_df[f"{seq_type}_top_ann_hit"] = np.nan
        blast_df[f"{seq_type}_annotation_source"] = 'blastn'

        blast_df = blast_df.drop_duplicates()

        return blast_df
    else:
        return unann_df


def get_ann(row, seq_type):
    priority_list = [
        f'{seq_type}_hits_UniVec',
        f'{seq_type}_hits_illumina_adapters',
        f'{seq_type}_hits_dfam_te_eukaryota',
        f'{seq_type}_hits_tncentral_te_prokaryotes_final',
        f'{seq_type}_hits_mge_aclame_genes_all_0',
        f'{seq_type}_hits_ice_iceberg',
        f'{seq_type}_hits_direct_repeats',
        f'{seq_type}_hits_spacers',
        f'{seq_type}_hits_eukaryota_its1_itstonedb',
        f'{seq_type}_hits_eukaryota_its2_biozentrum',
        f'{seq_type}_hits_RF_all',
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


def summarize(ann_table, df, seq_type, run_blast):

    c_list = [c for c in ann_table.columns if 'hits_pos' not in c]

    # ann_table = ann_table[c_list]

    # try:
    #     ann_table = ann_table.drop([f'{seq_type}_hits_concat_sherlock_fastas', axis=1)
    # except:
    #     pass
    ann_table['all_unannotated'] = (
        ann_table[c_list]
        .drop(seq_type, axis=1)
        .replace('*', np.NaN)
        .isnull()
        .all(axis=1)
    )

    ann_df = ann_table[ann_table['all_unannotated'] == False]
    unann_df = ann_table[ann_table['all_unannotated'] == True]

    ann_df[f"{seq_type}_top_ann"] = None
    ann_df[f"{seq_type}_top_ann_hit"] = None
    ann_df[f"{seq_type}_top_ann_hit_pos"] = None

    ann_df = ann_df.apply(lambda row: get_ann(row, seq_type), axis=1)

    # get the number of annotations
    ann_df[f"{seq_type}_num_ann"] = (
        ann_df[c_list]
        .drop(seq_type, axis=1)
        .replace('*', np.nan)
        .notnull()
        .astype(int)
        .sum(axis=1)
    )
    ann_df[f"{seq_type}_annotation_source"] = "bowtie2"

    # subset
    bowtie_df = ann_df[[seq_type, f"{seq_type}_top_ann", f"{seq_type}_top_ann_hit", f"{seq_type}_num_ann", f"{seq_type}_annotation_source"]]

    ## unannotated df
    unann_df = (
        pd.merge(
            df,
            unann_df,
            on=seq_type
        )
        [[seq_type, 'bf.cor.p']]
    )

    unann_df[f"{seq_type}_top_ann"] = None
    unann_df[f"{seq_type}_top_ann_hit"] = None
    unann_df[f"{seq_type}_num_ann"] = None
    unann_df[f"{seq_type}_annotation_source"] = None
    unann_df[f"{seq_type}_top_ann_hit_pos"] = None

    if run_blast:
        blast_df = get_blast_df(unann_df, seq_type)
    else:
        blast_df = unann_df

    blast_df = blast_df.drop_duplicates()

    bowtie_df.to_csv(f'bowtie_{seq_type}.csv', index=False, sep='\t')
    blast_df.to_csv(f'blast_{seq_type}.csv', index=False, sep='\t')

    summarized_df = pd.concat([bowtie_df, blast_df], axis=0)

    return summarized_df


def main():
    args = get_args()

    scores = pd.read_csv(args.anchor_scores, sep='\t')
    anchor_targets_counts = pd.read_csv(args.anchor_target_counts, sep='\t')
    ann_anchors = pd.read_csv(args.annotated_anchors, sep='\t')
    ann_targets = pd.read_csv(args.annotated_targets, sep='\t')


    anchor_targets_counts['total_target_counts'] = (
        anchor_targets_counts
        .drop(['anchor', 'target'], axis=1)
        .sum(axis=1)
    )

    df = pd.merge(
        anchor_targets_counts,
        scores,
        on='anchor'
    )

    df['rank_target_counts'] = (
        df
        .groupby('anchor')['total_target_counts']
        .rank('average')
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

    summarized_anchors = summarize(ann_anchors, df, 'anchor', args.run_blast)
    summarized_targets = summarize(ann_targets, df, 'target', args.run_blast)

    # df.to_csv('df.csv', index=False, sep='\t')

    df = pd.merge(df, summarized_anchors, on='anchor', how='left')
    # df.to_csv('summarized_anchors.csv', index=False, sep='\t')

    df = pd.merge(df, summarized_targets, on='target', how='left')
    # df.to_csv('summarized_targets.csv', index=False, sep='\t')

    df.drop_duplicates().to_csv(args.outfile, sep='\t', index=False)


main()
