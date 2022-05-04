#!/usr/bin/env python3

import pandas as pd
import numpy as np
import os

for seq_type in ["anchors", "targets"]:
    df = pd.read_csv(f"{seq_type}.txt", header=None, names=[seq_type])

    ann_tuples = [
        ["strand", f"{seq_type}_genome_strand.txt"],
        ["gene", f"{seq_type}_genome_genes.txt"],
        ["gene_MAPQ", f"{seq_type}_genome_mapq.txt"],
        ["distance_exon_start", f"{seq_type}_genome_exon_starts_distances.txt"],
        ["distance_exon_end", f"{seq_type}_genome_exon_ends_distances.txt"],
        ["transcript", f"{seq_type}_transcriptome_hit.txt"],
        ["transcriptome_MAPQ", f"{seq_type}_transcriptome_mapq.txt"],
    ]

    for ann_name, ann_file in ann_tuples:
        if os.path.exists(ann_file):
            ann = pd.read_csv(ann_file, sep='\t', header=None, names=[seq_type, ann_name])
            df = pd.merge(df, ann, on=seq_type, how='outer')
        else:
            df[ann_name] = np.nan

    df = df.replace('*', np.nan)

    df.to_csv(f"genome_annotations_{seq_type}.tsv", index=False, sep='\t')
