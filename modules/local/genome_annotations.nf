
process GENOME_ANNOTATIONS {

    tag "${index_name}"
    label 'process_low'

    input:
    path anchor_genome_bam
    path target_genome_bam
    path anchor_trans_bam
    path target_trans_bam
    path anchor
    path target
    path gene_bed
    path exon_starts_bed
    path exon_ends_bed

    output:
    path anchor_hits
    path target_hits

    script:
    anchor_hits         = "genome_annotations_anchor.tsv"
    target_hits         = "genome_annotations_target.tsv"
    """
    seq_types="anchor target"

    for seq_type in \${seq_types}
    do
        samtools view \${seq_type}_transcriptome.bam \\
            | cut -f1,3 \\
            > \${seq_type}_transcriptome_hit.txt

        samtools view \${seq_type}_transcriptome.bam \\
            | cut -f1,5 \\
            > \${seq_type}_transcriptome_mapq.txt

        bedtools bamtobed -i \${seq_type}_genome.bam \\
            | sort -k1,1 -k2,2n \\
            > \${seq_type}_genome.bed

        if [[ \$(wc -l \${seq_type}_genome.bed | awk '{print \$1}') -gt 0 ]]
        then
            cut -f4,6 \${seq_type}_genome.bed \\
                | sort \\
                > \${seq_type}_genome_strand.txt
            cut -f4,5 \${seq_type}_genome.bed \\
                | sort \\
                > \${seq_type}_genome_mapq.txt
            bedtools intersect -a \${seq_type}_genome.bed -b ${gene_bed} -wb \\
                | cut -f1-4,10,6 \\
                | bedtools groupby -i - -g 1,2,3,4,5 -c 6 -o collapse \\
                | cut -f4,6 \\
                | sort \\
                > \${seq_type}_genome_genes.txt
            bedtools closest -a \${seq_type}_genome.bed -b ${exon_starts_bed} -D ref -id -t first \\
                | cut -f4,13 \\
                | sort \\
                > \${seq_type}_genome_upstream_exon_starts.txt
            bedtools closest -a \${seq_type}_genome.bed -b ${exon_ends_bed} -D ref -id -t first \\
                | cut -f4,13 \\
                | sort \\
                > \${seq_type}_genome_upstream_exon_ends.txt
            bedtools closest -a \${seq_type}_genome.bed -b ${exon_starts_bed} -D ref -iu -t first \\
                | cut -f4,13 \\
                | sort \\
                > \${seq_type}_genome_downstream_exon_starts.txt
            bedtools closest -a \${seq_type}_genome.bed -b ${exon_ends_bed} -D ref -iu -t first \\
                | cut -f4,13 \\
                | sort \\
                > \${seq_type}_genome_downstream_exon_ends.txt
        fi
    done

    merge_genome_annotations.py
    """
}
