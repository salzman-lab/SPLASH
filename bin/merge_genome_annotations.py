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

    for record in SeqIO.parse(f"{args.fasta_name}.fasta", 'fasta'):
        reads.append(str(record.id))
    df = pd.DataFrame(reads, columns=[args.fasta_name])


    # ann_tuples = [
    #     ["strand", f"{args.fasta_name}_genome_strand.txt"],
    #     ["gene", f"{args.fasta_name}_genome_genes.txt"],
    #     ["gene_MAPQ", f"{args.fasta_name}_genome_mapq.txt"],
    #     ["distance_upstream_exon_start", f"{args.fasta_name}_genome_upstream_exon_starts.txt"],
    #     ["distance_upstream_exon_end", f"{args.fasta_name}_genome_upstream_exon_ends.txt"],
    #     ["distance_downstream_exon_start", f"{args.fasta_name}_genome_downstream_exon_starts.txt"],
    #     ["distance_downstream_exon_end", f"{args.fasta_name}_genome_downstream_exon_ends.txt"],
    #     ["transcript", f"{args.fasta_name}_transcriptome_hit.txt"],
    #     ["transcriptome_MAPQ", f"{args.fasta_name}_transcriptome_mapq.txt"],
    # ]

    #######
    ### Super messy now for diagnostic purposes
    #######

    ann_tuples = [
        ["local_strand", f"local/{args.fasta_name}_genome_strand.txt"],
        ["local_gene", f"local/{args.fasta_name}_genome_genes.txt"],
        ["local_gene_MAPQ", f"local/{args.fasta_name}_genome_mapq.txt"],
        ["local_distance_upstream_exon_start", f"local/{args.fasta_name}_genome_upstream_exon_starts.txt"],
        ["local_distance_upstream_exon_end", f"local/{args.fasta_name}_genome_upstream_exon_ends.txt"],
        ["local_distance_downstream_exon_start", f"local/{args.fasta_name}_genome_downstream_exon_starts.txt"],
        ["local_distance_downstream_exon_end", f"local/{args.fasta_name}_genome_downstream_exon_ends.txt"],
        ["local_transcript", f"local/{args.fasta_name}_transcriptome_hit.txt"],
        ["local_transcriptome_MAPQ", f"local/{args.fasta_name}_transcriptome_mapq.txt"],
        ["end_to_end_strand", f"end_to_end/{args.fasta_name}_genome_strand.txt"],
        ["end_to_end_gene", f"end_to_end/{args.fasta_name}_genome_genes.txt"],
        ["end_to_end_gene_MAPQ", f"end_to_end/{args.fasta_name}_genome_mapq.txt"],
        ["end_to_end_distance_upstream_exon_start", f"end_to_end/{args.fasta_name}_genome_upstream_exon_starts.txt"],
        ["end_to_end_distance_upstream_exon_end", f"end_to_end/{args.fasta_name}_genome_upstream_exon_ends.txt"],
        ["end_to_end_distance_downstream_exon_start", f"end_to_end/{args.fasta_name}_genome_downstream_exon_starts.txt"],
        ["end_to_end_distance_downstream_exon_end", f"end_to_end/{args.fasta_name}_genome_downstream_exon_ends.txt"],
        ["end_to_end_transcript", f"end_to_end/{args.fasta_name}_transcriptome_hit.txt"],
        ["end_to_end_transcriptome_MAPQ", f"end_to_end/{args.fasta_name}_transcriptome_mapq.txt"],
    ]

    for ann_name, ann_file in ann_tuples:
        if os.path.exists(ann_file):
            ann = pd.read_csv(ann_file, sep='\t', header=None, names=[args.fasta_name, ann_name])

            # clean up
            if "MAPQ" in ann_name:
                ann = ann.replace(0.0, np.nan)

            df = pd.merge(df, ann, on=args.fasta_name, how='outer')
        else:
            df[ann_name] = np.nan

    df = df.replace('*', np.nan)

    df.to_csv(f"genome_annotations_{args.fasta_name}.tsv", index=False, sep='\t')


main()
